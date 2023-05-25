"""
KD树近似查询库
一个合理的使用方法：
1. 调用 loadDataset 加载数据集
2. 调用 buildKDTrees 构建KD树
3. 调用 loadModels 加载模型
4. 使用 query 函数进行查询
"""
import ctypes
from ctypes import CDLL, POINTER, Structure, c_int, c_float, c_double
import numpy as np
import pandas as pd
import os.path as osp
import os
from timeit import default_timer as timer
import atexit
from itertools import combinations, product
from tqdm import tqdm
from scipy.special import comb
import shutil

WORKING_DIR = osp.dirname(osp.dirname(osp.abspath(__file__)))
# WORKING_DIR = '*/2021201626/'
HANDOUTS_DIR = osp.dirname(WORKING_DIR)
DATA_DIR = osp.join(WORKING_DIR, "tmp")  # 存放临时数据的目录，很重要，要保证至少有1GB的空间（10GB以上为佳）
MODEL_DIR = osp.join(DATA_DIR, "models")

shutil.rmtree(DATA_DIR, ignore_errors=True)
os.mkdir(DATA_DIR)
# if not osp.exists(DATA_DIR):
#     os.mkdir(DATA_DIR)

CODE_DIR = osp.join(WORKING_DIR, "codes")
DATASET_DIR = osp.join(HANDOUTS_DIR, "data")
WORKLOAD_DIR = osp.join(HANDOUTS_DIR, "workload")

TINY_DATASET_PATH = osp.join(DATASET_DIR, "idebench-dataset-tiny.csv")
DATASET_PATH = osp.join(DATASET_DIR, "idebench-dataset.csv")

COLUMNS = [
    "YEAR_DATE",
    "DEP_DELAY",
    "TAXI_OUT",
    "TAXI_IN",
    "ARR_DELAY",
    "AIR_TIME",
    "DISTANCE",
    "UNIQUE_CARRIER",
    "ORIGIN",
    "ORIGIN_STATE_ABR",
    "DEST",
    "DEST_STATE_ABR",
]
DISCRETE_COLUMNS = [
    "UNIQUE_CARRIER",
    "ORIGIN",
    "ORIGIN_STATE_ABR",
    "DEST",
    "DEST_STATE_ABR",
]
CONTINUOUS_COLUMNS = [
    "YEAR_DATE",
    "DEP_DELAY",
    "TAXI_OUT",
    "TAXI_IN",
    "ARR_DELAY",
    "AIR_TIME",
    "DISTANCE",
]
Q_ID_LIST = ["Q1", "Q2", "Q3", "Q4"]

COLUMN2INDEX = {c: i for i, c in enumerate(COLUMNS)}
COLUMN2INDEX.update({"*": -1})
COLUMN2INDEX.update({"_None_": -1})
INDEX2COLUMN = COLUMNS

DATASET, VALUE2ID, ID2VALUE = None, None, None

global_mode = "memory"  # 'memory' or 'performance'


def set_mode(_mode):
    global global_mode
    global_mode = _mode


def get_dataset():
    df = pd.read_csv(DATASET_PATH)
    df = df[COLUMNS]
    return df


def get_workload_path(q_id="tiny"):
    """offline, online"""
    offline_workload_path = osp.join(
        WORKLOAD_DIR, "{}-workload-{}.json".format("offline", q_id)
    )
    online_workload_path = osp.join(
        WORKLOAD_DIR, "{}-workload-{}.json".format("online", q_id)
    )
    return offline_workload_path, online_workload_path


def get_workloads(online=True, offline=True, q_id_list=Q_ID_LIST):
    assert online or offline
    workloads = []
    for q_id in q_id_list:
        offline_workload_path, online_workload_path = get_workload_path(q_id)
        if online:
            online_workloads = pd.read_json(
                online_workload_path, orient="records", lines=True
            )
            workloads.append(online_workloads)
        if offline:
            offline_workloads = pd.read_json(
                offline_workload_path, orient="records", lines=True
            )
            workloads.append(offline_workloads)
    return pd.concat(workloads, ignore_index=True)


def factorize(df, columns=DISCRETE_COLUMNS):
    id2value = {}
    value2id = {}
    df = df.copy()
    for col in columns:
        values, index = pd.factorize(df[col])
        df[col] = values
        col_id = COLUMN2INDEX[col]
        id2value.update({col_id: index.values})
        value2id.update({col_id: {v: i for i, v in enumerate(index)}})
        id2value.update({col: index.values})
        value2id.update({col: {v: i for i, v in enumerate(index)}})

    return df, id2value, value2id


def get_modelName(predicate, groupBy) -> str:
    col = [int(p[0]) for p in predicate]
    if groupBy not in col and groupBy != -1:
        col.append(groupBy)
    col = sorted(col)
    modelName = "_".join([str(c) for c in col])
    return modelName


def standlize(workloads: pd.DataFrame):
    """
    数值化并排序，排序是为了让使用相同模型的查询相邻，避免反复加载模型。
    """

    # result_col
    op_map = {"count": 0, "sum": 1, "avg": 2}
    standlized_result_col = []
    for result_col in workloads["result_col"]:
        standlized_result_col.append(
            [(op_map[i[0]], COLUMN2INDEX[i[1]]) for i in result_col if len(i) >= 2]
        )
    workloads["result_col"] = standlized_result_col

    # predicate
    standlized_predicate = []
    for predicate in workloads["predicate"]:
        predicate_tmp = []
        for i in predicate:
            col = COLUMN2INDEX[i["col"]]
            if i["col"] in DISCRETE_COLUMNS:
                lb = ub = VALUE2ID[i["col"]][i["lb"]]
            else:
                if i["lb"] == "_None_":
                    lb = 0
                else:
                    lb = float(i["lb"])
                if i["ub"] == "_None_":
                    ub = 1e9
                else:
                    ub = float(i["ub"])
            predicate_tmp.append((col, lb, ub))
        standlized_predicate.append(predicate_tmp)
    workloads["predicate"] = standlized_predicate

    # groupby
    standlized_groupby = workloads["groupby"].apply(COLUMN2INDEX.get)
    workloads["groupby"] = standlized_groupby

    # model
    workloads["model"] = workloads.apply(
        lambda x: get_modelName(x["predicate"], x["groupby"]), axis=1
    )

    # sort
    arg = workloads["model"].argsort()
    workloads = workloads.iloc[arg, :]

    return workloads, arg


def _get_valueMaps_from_dataset(dataset):
    dataset, id2value, value2id = factorize(dataset)
    maps = [value2id, id2value]
    np.save(osp.join(DATA_DIR, "valueMaps.npy"), maps)
    dataset.to_csv(osp.join(DATA_DIR, "factorizedDataset.csv"), index=False)
    return dataset, value2id, id2value


def map_factorized_init(dataset=None):
    global DATASET, VALUE2ID, ID2VALUE
    if osp.exists(osp.join(DATA_DIR, "valueMaps.npy")) and osp.exists(
        osp.join(DATA_DIR, "factorizedDataset.csv")
    ):
        DATASET = pd.read_csv(osp.join(DATA_DIR, "factorizedDataset.csv"))
        VALUE2ID, ID2VALUE = np.load(
            osp.join(DATA_DIR, "valueMaps.npy"), allow_pickle=True
        )
    elif dataset is not None:
        DATASET, VALUE2ID, ID2VALUE = _get_valueMaps_from_dataset(dataset)


lib = CDLL(osp.join(CODE_DIR, "libaqp.so"))


class Operation(Structure):
    _fields_ = [("col", c_int), ("op", c_int)]


class Predication(Structure):
    _fields_ = [("col", c_int), ("lb", c_float), ("ub", c_float)]


class GroupAnswer(Structure):
    _fields_ = [("id", c_int), ("value", c_float)]


class Answer(Structure):
    _fields_ = [("group_ans", POINTER(GroupAnswer)), ("size", c_int)]


class AnswerBatch(Structure):
    _fields_ = [("ans", POINTER(Answer)), ("size", c_int)]


def loadDataset(dataset):
    if not osp.exists(osp.join(DATA_DIR, "factorizedDataset.csv")):
        dataset = dataset[COLUMNS]
        map_factorized_init(dataset)
        dataset = DATASET.copy()
        values = dataset.values.astype(np.float32)
        values = values.flatten()
        lib.loadData(values.ctypes.data_as(POINTER(c_float)), dataset.shape[0])


def buildKDTrees(force=True, deltaDepth=None, buildK=None):
    if not osp.exists(MODEL_DIR):
        os.mkdir(MODEL_DIR)
    mode = global_mode
    if not force and osp.exists(osp.join(MODEL_DIR, "model_list.txt")):
        return
    print(mode, deltaDepth, buildK)
    if osp.exists(osp.join(MODEL_DIR, "model_list.txt")):
        os.remove(osp.join(MODEL_DIR, "model_list.txt"))
    if mode == "performance":
        if deltaDepth is None:
            deltaDepth = -3
        if buildK is None:
            buildK = 0.1
        for pred_num in [1, 2, 3]:
            for col in tqdm(
                combinations(range(0, 12), pred_num), total=int(comb(12, pred_num))
            ):
                col_np = np.array(col, dtype=np.int32)
                lib.build(
                    col_np.ctypes.data_as(POINTER(c_int)),
                    len(col_np),
                    deltaDepth,
                    buildK,
                )
    else:  # mode == 'memory'
        if deltaDepth is None:
            deltaDepth = 1
        if buildK is None:
            buildK = 1
        for di_pred_num in [0, 1, 2, 3]:
            for di_col in tqdm(
                combinations(range(7, 12), di_pred_num), total=int(comb(5, di_pred_num))
            ):
                if di_pred_num < 3:
                    col_np = np.r_[
                        np.array(range(7), dtype=np.int32),
                        np.array(di_col, dtype=np.int32),
                    ]
                else:
                    col_np = np.array(di_col, dtype=np.int32)
                lib.build(
                    col_np.ctypes.data_as(POINTER(c_int)),
                    len(col_np),
                    deltaDepth,
                    buildK,
                )


def query(workload):
    if global_mode == "performance":
        mode = 0
    else:  # mode == 'memory'
        mode = 1

    ops = np.array(workload["result_col"], dtype=Operation)
    preds = np.array(workload["predicate"], dtype=Predication)
    groupBy_col = workload["groupby"]
    # print(ops, preds, groupBy_col, workload['ground_truth'])
    ans = lib.aqpQuery(
        ops.ctypes.data_as(POINTER(Operation)),
        len(ops),
        preds.ctypes.data_as(POINTER(Predication)),
        len(preds),
        groupBy_col,
        mode,
    )
    ans = ans.contents
    size = ans.size
    ret = []
    for i in range(size):
        g_ans = ans.group_ans[i]
        if g_ans.id < 0:
            ret.append([g_ans.value])
        else:
            id = g_ans.id
            ret.append([ID2VALUE[groupBy_col][id], g_ans.value])
    return ret


@atexit.register
def clear():
    lib.clear()
    shutil.rmtree(DATA_DIR)
    print("clear")


def lib_init():
    lib.aqpQuery.argtypes = [
        POINTER(Operation),
        c_int,
        POINTER(Predication),
        c_int,
        c_int,
        c_int,
    ]
    lib.aqpQuery.restype = POINTER(Answer)

    lib.loadData.argtypes = [POINTER(c_float), c_int]
    lib.loadData.restype = None

    lib.clear.argtypes = []
    lib.clear.restype = None

    lib.build.argtypes = [POINTER(c_int), c_int, c_int, c_float]
    lib.build.restype = None

    lib.init.argtypes = [ctypes.c_char_p]
    lib.init.restype = None
    dir = MODEL_DIR
    lib.init(dir.encode("utf-8"))

    map_factorized_init(None)


def loadModels():
    lib.load_models()


lib_init()
