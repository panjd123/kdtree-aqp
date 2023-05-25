import pandas as pd
import json
from tqdm import tqdm
import sys
import os.path as osp
import numpy as np
lib_dir = osp.join(osp.dirname(osp.abspath(__file__)), 'codes')
sys.path.append(lib_dir)
import kdtree_aqp as aqplib

# MODE = 'performance' # 'memory' or 'performance'
# MODE = 'memory' # 'memory' or 'performance'
# aqplib.set_mode(MODE)

def aqp_online(data: pd.DataFrame, Q: list) -> list:
    workloads = pd.json_normalize([json.loads(i) for i in Q]) #把Q转换为DataFrame，方便处理
    
    workloads, arg = aqplib.standlize(workloads) # 数值化+排序
    
    results = []
    
    for _, workload in tqdm(workloads.iterrows()):
        ans = aqplib.query(workload)
        results.append(ans)
    
    arg = np.argsort(arg)
    results = [results[i] for i in arg]
    
    results = [json.dumps(i, ensure_ascii=False) for i in results]
    
    return results

def aqp_offline(data: pd.DataFrame, Q: list) -> None: 
    '''无需返回任何值
    必须编写的aqp_offline函数，用data和Offline-workload，构建采样、索引、机器学习相关的结构或模型
    data是一个DataFrame，Q是一个包含多个json格式字符串的list
    如果你的算法无需构建结构或模型，该函数可以为空
    '''
    aqplib.loadDataset(data)
    aqplib.buildKDTrees(force=False) # 构造KD树
    aqplib.loadModels() # 加载到内存