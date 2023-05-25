#include "libaqp.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

static bool is_init = false;

/**** dataset ****/

static FLOAT_T* dataset = nullptr;
static int dataset_size;
static DATA_T* data = nullptr;
static int& data_size = dataset_size;

static size_t working_memory = 0;
static size_t max_working_memory = 0;
static size_t total_memory = 0;

/**** Build KD-Tree ****/

static COL_T split_axises[12];
static int split_axis_num;
static int max_depth = 20;
static FILE* model_file;

// k larger, performance better
// k smaller, accuracy better
static float build_k = 0.1;

double data_cross_ratio(const BOUND_T& bound_in, const BOUND_T& bound_out) {
    double ratio = 1;
    for (int i = 0; i < DATA_DIM; i++) {
        if (bound_in[i][0] == bound_in[i][1]) {
            ratio *= (bound_out[i][0] <= bound_in[i][0] && bound_in[i][0] <= bound_out[i][1]);
        } else {
            double l, r;
            l = std::max(bound_in[i][0], bound_out[i][0]);
            r = std::min(bound_in[i][1], bound_out[i][1]);
            ratio *= (r - l) / (bound_in[i][1] - bound_in[i][0]);
        }
    }
    return ratio;
}

int kd_contain(const BOUND_T& bound_in, const BOUND_T& bound_out) {
    for (int i = 0; i < split_axis_num; i++) {
        int split_axis = split_axises[i];
        if (bound_in[split_axis][0] < bound_out[split_axis][0] ||
            bound_in[split_axis][1] > bound_out[split_axis][1]) {
            return 0;
        }
    }
    return 1;
}

int kd_cross(const BOUND_T& bound_in, const BOUND_T& bound_out) {
    for (int i = 0; i < split_axis_num; i++) {
        int split_axis = split_axises[i];
        if (bound_in[split_axis][0] > bound_out[split_axis][1] ||
            bound_in[split_axis][1] < bound_out[split_axis][0]) {
            return 0;
        }
    }
    return 1;
}

Node* buildKDTree(DATA_T* data, int l, int r, int depth) {
    if (l > r) {
        return nullptr;
    }
    Node* u = new Node;

    if (l == r || depth >= max_depth || split_axis_num == 0) {
        u->lchild = u->rchild = nullptr;
        // u->split_value = 0;
        u->count = r - l + 1;
        for (int i = 0; i < DATA_DIM; i++) {
            u->sum[i] = 0;
            u->bound[i][0] = 1e9;
            u->bound[i][1] = -1e9;
            for (int j = l; j <= r; j++) {
                u->sum[i] += data[j][i];
                u->bound[i][0] = std::min(u->bound[i][0], data[j][i]);
                u->bound[i][1] = std::max(u->bound[i][1], data[j][i]);
            }
        }
        return u;
    }

    COL_T split_dim = split_axises[depth % split_axis_num];
    // std::nth_element(data + l, data + (l + r) / 2, data + r + 1,
    //                  [split_dim](const DATA_T& a, const DATA_T& b) {
    //                      return a[split_dim] < b[split_dim];
    //                  });
    int perfomance_median = (l + r) / 2;

    /* get the position of the mid value */
    int accuracy_median = l;

    if (build_k != 1) {
        FLOAT_T min(1e9), max(-1e9), mid_value;
        for (int i = l; i <= r; i++) {
            if (data[i][split_dim] < min) {
                min = data[i][split_dim];
            }
            if (data[i][split_dim] > max) {
                max = data[i][split_dim];
            }
        }
        mid_value = (min + max) / 2;
        for (int i = l; i <= r; i++) {
            if (data[i][split_dim] < mid_value) {
                accuracy_median++;
            }
        }
    }

    // consider both performance and accuracy
    int median = perfomance_median * build_k + accuracy_median * (1 - build_k);

#ifdef INFO
    printf("%.2f %.2f %.2f\n", min, mid_value, max);
    printf("perfomance_median: %d, accuracy_median: %d, median: %d\n", perfomance_median, accuracy_median, median);
#endif

    if (median == r) {
        u->lchild = u->rchild = nullptr;
        // u->split_value = 0;
        u->count = r - l + 1;
        for (int i = 0; i < DATA_DIM; i++) {
            u->sum[i] = 0;
            u->bound[i][0] = 1e9;
            u->bound[i][1] = -1e9;
            for (int j = l; j <= r; j++) {
                u->sum[i] += data[j][i];
                u->bound[i][0] = std::min(u->bound[i][0], data[j][i]);
                u->bound[i][1] = std::max(u->bound[i][1], data[j][i]);
            }
        }
        return u;
    }

    std::nth_element(data + l, data + median, data + r + 1,
                     [split_dim](const DATA_T& a, const DATA_T& b) {
                         return a[split_dim] < b[split_dim];
                     });

    // u->split_value = data[median][split_dim];
    u->count = r - l + 1;

    // [l, median]
    u->lchild = buildKDTree(data, l, median, depth + 1);
    // [median + 1, r]
    u->rchild = buildKDTree(data, median + 1, r, depth + 1);

    for (int i = 0; i < DATA_DIM; i++) {
        u->sum[i] = 0;
        u->bound[i][0] = 1e9;
        u->bound[i][1] = -1e9;
        if (u->lchild != nullptr) {
            u->sum[i] += u->lchild->sum[i];
            u->bound[i][0] = std::min(u->bound[i][0], u->lchild->bound[i][0]);
            u->bound[i][1] = std::max(u->bound[i][1], u->lchild->bound[i][1]);
        }
        if (u->rchild != nullptr) {
            u->sum[i] += u->rchild->sum[i];
            u->bound[i][0] = std::min(u->bound[i][0], u->rchild->bound[i][0]);
            u->bound[i][1] = std::max(u->bound[i][1], u->rchild->bound[i][1]);
        }
    }

    return u;
}

void _queryRange(Node* u, const BOUND_T& bound, FLOAT_T* sum, double& count) {
    if (u == nullptr) {
        return;
    }
    if (IS_LEAF(u) || kd_contain(u->bound, bound)) {
        double ratio = data_cross_ratio(u->bound, bound);
        count += u->count * ratio;
        for (int i = 0; i < DATA_DIM; i++) {
            sum[i] += u->sum[i] * ratio;
        }
        return;
    }
    if (u->lchild && kd_cross(u->lchild->bound, bound)) {
        _queryRange(u->lchild, bound, sum, count);
    }
    if (u->rchild && kd_cross(u->rchild->bound, bound)) {
        _queryRange(u->rchild, bound, sum, count);
    }
}

void saveKDTree(Node* u) {
    if (!model_file || !u) {
        return;
    }
    fwrite(u, sizeof(Node), 1, model_file);
    saveKDTree(u->lchild);
    saveKDTree(u->rchild);
}

void saveIDKDTree(Node* u, int id) {
    if (!model_file || !u) {
        return;
    }
    fwrite(&id, sizeof(int), 1, model_file);
    saveKDTree(u);
}

Node* loadKDTree() {
    Node* u = new Node;
    working_memory += sizeof(Node);
    size_t _ = fread(u, sizeof(Node), 1, model_file);
    if (_ < 1) {
        printf("loadKDTree error\n");
        return nullptr;
    }
    if (u->lchild != nullptr) {
        u->lchild = loadKDTree();
    }
    if (u->rchild != nullptr) {
        u->rchild = loadKDTree();
    }
    return u;
}

void clearKDTree(Node* u) {
    if (u == nullptr) {
        return;
    }
    working_memory += sizeof(Node);
    clearKDTree(u->lchild);
    clearKDTree(u->rchild);
    delete u;
}

/**** Query Function ****/

// >=0 for continuous, -1 for discrete
static int col_map[COL_NUM];
static int value_num[COL_NUM];
static std::unordered_map<std::string, std::unordered_map<int, Node*>> model_map;
static std::vector<std::string> model_list;

std::string get_model_name(COL_VALUE_T& col_value) {
    std::string model_name = "";
    model_name.reserve(100);
    for (int i = 0; i < col_value.size(); i++) {
        model_name += std::to_string(col_value[i].first);
        // if (col_value[i].second >= 0) {
        //     model_name += "=" + std::to_string(col_value[i].second);
        // }
        if (i != col_value.size() - 1) {
            model_name += "_";
        }
    }
    return model_name;
}

Node* get_root(COL_VALUE_T col_value) {
    int size = col_value.size();
    int root_idx = 0;
    int tmp = 1;
    for (int i = 0; i < size; i++) {
        int col = col_value[i].first;
        int value = col_value[i].second;
        if (value >= 0) {
            root_idx += tmp * value;
            tmp *= value_num[col];
        }
    }
    std::string model_name = get_model_name(col_value);
    if (model_map.find(model_name) == model_map.end()) {
// printf("model %s not found\n", model_name.c_str());
#if 0
        printf("load %s\n", model_name.c_str());
#endif
        load_model(model_name);
        model_list.push_back(model_name);
    }
    Node* root = model_map[model_name][root_idx];
    return root;
}

void load_col_type() {
    // FILE* col_file = fopen("col.txt", "rb");
    // fread(col_map, sizeof(int), COL_NUM, col_file);
    // fread(value_num, sizeof(int), COL_NUM, col_file);
    for (int i = 0; i < 7; i++)
        col_map[i] = i;
    for (int i = 7; i < COL_NUM; i++)
        col_map[i] = -1;
    for (int i = 0; i < 7; i++)
        value_num[i] = 1;
    value_num[7] = 26;
    value_num[8] = 363;
    value_num[9] = 53;
    value_num[10] = 366;
    value_num[11] = 53;
}

static std::string MODEL_DIR;

std::string get_model_path(std::string model_name) {
    return MODEL_DIR + "/model_" + model_name + ".bin";
}

void load_model(const std::string model_name) {
    if (model_map.find(model_name) != model_map.end()) {
        return;
    }
    while (total_memory > MEM_LIMIT && model_list.size() > 0) {
        int rd = rand() % model_list.size();
        clear_model(model_list[rd]);
        model_list.erase(model_list.begin() + rd);
    }
    std::string model_path = get_model_path(model_name);
    model_file = fopen(model_path.c_str(), "rb");
    model_map[model_name] = std::unordered_map<int, Node*>();
    int idx;
    working_memory = 0;
    while (fread(&idx, sizeof(int), 1, model_file) == 1) {
        model_map[model_name][idx] = loadKDTree();
    }
    max_working_memory = std::max(max_working_memory, working_memory);
    total_memory += working_memory;
}

void clear_model(const std::string model_name) {
    if (model_map.find(model_name) == model_map.end()) {
        return;
    }
    working_memory = 0;
    for (auto& it : model_map[model_name]) {
        clearKDTree(it.second);
    }
    model_map.erase(model_name);
    total_memory -= working_memory;
}

extern "C" void load_models() {
    FILE* model_list_file = fopen((MODEL_DIR + "/model_list.txt").c_str(), "r");
    char model_name[100];
    std::vector<std::string> models;
    while (fscanf(model_list_file, "%s", model_name) != EOF) {
        models.push_back(model_name);
    }
    if (model_map.size() == models.size()) {
        return;
    }
    clear_models();
    printf("After clear: %zu\n", model_map.size());
    for (int i = 0; i < models.size() && total_memory < MEM_LIMIT; i++) {
        load_model(models[i]);
        model_list.push_back(models[i]);
    }
    printf("After load: %zu\n", model_map.size());
}

extern "C" void init(const char* dir) {
    if (!is_init) {
        MODEL_DIR = dir;
        load_col_type();
        is_init = true;
    }
}

void queryRange(Node* root, BOUND_T& bound, FLOAT_T* sum, double& count) {
#if 0
    printf("queryRange\n");
#endif
    memset(sum, 0, sizeof(FLOAT_T) * DATA_DIM);
    count = 0;
    _queryRange(root, bound, sum, count);
}

// according to the predication, extract the bound and col_value for query
void extract_pred(Predication* pred,
                  int pred_num,
                  BOUND_T& bound,
                  COL_VALUE_T& col_value,
                  MODE mode = MODE::PERFORMANCE) {
#if 0
    printf("extract pred\n");
#endif
    for (int i = 0; i < DATA_DIM; i++) {
        bound[i][0] = 1e9;
        bound[i][1] = -1e9;
    }
    if (mode == MODE::PERFORMANCE) {
        for (int i = 0; i < pred_num; i++) {
            if (col_map[pred[i].col] >= 0) {  // continuous
                bound[col_map[pred[i].col]][0] = pred[i].lb;
                bound[col_map[pred[i].col]][1] = pred[i].ub;
                if (split_axis_num <= 2) {  // largest model's split axis num is 3
                    // model_name += std::to_string(pred[i].col) + "_";
                    col_value.push_back(std::make_pair(pred[i].col, -1));
                    split_axises[split_axis_num] = col_map[pred[i].col];
                    split_axis_num++;
                }
            } else {
                col_value.push_back(std::make_pair(pred[i].col, int(pred[i].lb)));
            }
        }
    } else {
        int count_c = 0;
        for (int i = 0; i < 7; i++) {
            split_axises[i] = i;
            col_value.push_back(std::make_pair(i, -1));
        }
        split_axis_num = 7;
        for (int i = 0; i < pred_num; i++) {
            if (col_map[pred[i].col] >= 0) {  // continuous
                bound[col_map[pred[i].col]][0] = pred[i].lb;
                bound[col_map[pred[i].col]][1] = pred[i].ub;
                count_c++;
            } else {
                col_value.push_back(std::make_pair(pred[i].col, int(pred[i].lb)));
            }
        }
        if (count_c == 0 && pred_num == 3) {
            // 取第7位之后
            COL_VALUE_T col_value_tmp(col_value.begin() + 7, col_value.end());
            col_value = col_value_tmp;
            split_axis_num = 0;
        }
    }
    for (int i = 0; i < DATA_DIM; i++) {
        if (bound[i][0] == 1e9)
            bound[i][0] = -1e9;
        if (bound[i][1] == -1e9)
            bound[i][1] = 1e9;
    }
}

static Answer* _lastAns;

void clearans() {
    if (_lastAns != nullptr) {
        delete[] _lastAns->group_ans;
        delete _lastAns;
        _lastAns = nullptr;
    }
}

Answer* aqp_group_query(Predication* pred,
                        int pred_num,
                        Operation* ops,
                        int op_num,
                        COL_T groupBy_col,
                        MODE mode = MODE::PERFORMANCE) {
    clearans();

    BOUND_T bound;
    COL_VALUE_T col_value;

    FLOAT_T* sum = new FLOAT_T[DATA_DIM];
    double count = 0;
    split_axis_num = 0;

    Answer* ans = new Answer();
    _lastAns = ans;

    extract_pred(pred, pred_num, bound, col_value, mode);

#if 0
    printf("query init\n");
#endif
    if (groupBy_col != -1) {
        int groupBy_col_vectorIndex = -1;
        bool in_pred = !std::none_of(col_value.begin(), col_value.end(),
                                     [groupBy_col](std::pair<int, int> p) { return p.first == groupBy_col; });
        if (!in_pred) {
            ans->size = value_num[groupBy_col] * op_num;
            col_value.push_back(std::make_pair(groupBy_col, -1));
        } else {
            ans->size = op_num;
        }
        ans->group_ans = new GroupAnswer[ans->size];
        std::sort(col_value.begin(), col_value.end());
        groupBy_col_vectorIndex = std::find_if(col_value.begin(), col_value.end(),
                                               [groupBy_col](std::pair<int, int> p) { return p.first == groupBy_col; }) -
                                  col_value.begin();
        int bg = 0, ed = value_num[groupBy_col];
        if (in_pred) {
            bg = col_value[groupBy_col_vectorIndex].second;
            ed = bg + 1;
        }
        for (int i = bg; i < ed; i++) {
            col_value[groupBy_col_vectorIndex].second = i;
            Node* root = get_root(col_value);
            queryRange(root, bound, sum, count);
            for (int j = 0; j < op_num; j++) {
                int idx = in_pred ? j : i * op_num + j;
                ans->group_ans[idx].id = i;
                switch (ops[j].op) {
                    case OP::SUM:
                        ans->group_ans[idx].value = round(sum[ops[j].col] * 10) / 10;
                        break;
                    case OP::AVG:
                        if (count == 0) {
                            ans->group_ans[idx].value = 1;
                        } else {
                            ans->group_ans[idx].value = sum[ops[j].col] / count;
                        }
                        break;
                    case OP::COUNT:
                        ans->group_ans[idx].value = round(count);
                        break;
                    default:
                        break;
                }
            }
        }
    } else {
        ans->size = op_num;
        ans->group_ans = new GroupAnswer[ans->size];
        std::sort(col_value.begin(), col_value.end());
        Node* root = get_root(col_value);
        queryRange(root, bound, sum, count);
        for (int j = 0; j < op_num; j++) {
            ans->group_ans[j].id = -1;
            switch (ops[j].op) {
                case OP::SUM:
                    ans->group_ans[j].value = round(sum[ops[j].col] * 10) / 10;
                    break;
                case OP::AVG:
                    if (count == 0) {
                        ans->group_ans[j].value = 1;
                    } else {
                        ans->group_ans[j].value = sum[ops[j].col] / count;
                    }
                    break;
                case OP::COUNT:
                    ans->group_ans[j].value = round(count);
                    break;
                default:
                    break;
            }
        }
    }
    delete[] sum;
    return ans;
}

/**** KDTree Test Function ****/

void printKDTree(Node* u, int depth) {
    if (u->rchild) {
        printKDTree(u->lchild, depth + 1);
    }
    for (int i = 0; i < depth; i++) {
        printf("  ");
    }
    printf("[%d %.0lf]\n", u->count, u->sum[0]);
    if (u->lchild) {
        printKDTree(u->rchild, depth + 1);
    }
}

void testKDTree(Node* u, int depth) {
    if (u->lchild) {
        testKDTree(u->lchild, depth + 1);
    }
    if (u->rchild) {
        testKDTree(u->rchild, depth + 1);
    }
    if (IS_LEAF(u) && depth != 0) {
        if (u->count != 1) {
            printf("error\n");
        }
    }
}

extern "C" void loadData(FLOAT_T* _data, int n) {
    clearData();

    dataset_size = n;

    dataset = new FLOAT_T[n * COL_NUM];
    for (int i = 0; i < n * COL_NUM; i++) {
        dataset[i] = _data[i];
    }

    data = new DATA_T[n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < DATA_DIM; j++) {
            data[i][j] = _data[i * COL_NUM + j];
        }
    }
}

void clearData() {
    if (dataset)
        delete[] dataset;
    if (data)
        delete[] data;
}

extern "C" void build(INT_T* col, int size, int delta_depth, float _build_k) {
    build_k = _build_k;
    FILE* model_list_file = fopen((MODEL_DIR + "/model_list.txt").c_str(), "a");
    DATA_T* tmp_data = new DATA_T[dataset_size];
    int tmp_data_size = dataset_size;

    split_axis_num = 0;

    int discrete_axises[12], discrete_axis_num = 0;

    for (int i = 0; i < size; i++) {
        int c = col[i];
        if (IS_CONTINUED(c)) {
            split_axises[split_axis_num] = c;
            split_axis_num++;
        } else {
            discrete_axises[discrete_axis_num] = c;
            discrete_axis_num++;
        }
    }

    // std::copy(data, data + dataset_size, tmp_data);
    int* d = new int[dataset_size];
    for (int i = 0; i < dataset_size; i++)
        d[i] = i;

    std::sort(d, d + dataset_size, [&](int a, int b) {
        for (int i = 0; i < discrete_axis_num; i++) {
            if (dataset[a * COL_NUM + discrete_axises[i]] != dataset[b * COL_NUM + discrete_axises[i]]) {
                return dataset[a * COL_NUM + discrete_axises[i]] < dataset[b * COL_NUM + discrete_axises[i]];
            }
        }
        return false;
    });

    for (int i = 0; i < dataset_size; i++) {
        for (int j = 0; j < DATA_DIM; j++) {
            tmp_data[i][j] = dataset[d[i] * COL_NUM + j];
        }
    }
    COL_VALUE_T col_value;
    for (int i = 0; i < size; i++) {
        col_value.push_back(std::make_pair(col[i], -1));
    }
    std::string model_name = get_model_name(col_value);
    std::string model_path = get_model_path(model_name);

    model_file = fopen(model_path.c_str(), "wb");

    int l = 0, r = 0;

    while (l != tmp_data_size) {
        /* get idx */
        int idx = 0;
        int tmp = 1;
        for (int i = 0; i < discrete_axis_num; i++) {
            idx = idx + tmp * dataset[d[l] * COL_NUM + discrete_axises[i]];
            tmp *= value_num[discrete_axises[i]];
        }

        /* get r */
        r = l;
        while (r < tmp_data_size - 1) {
            for (int i = 0; i < discrete_axis_num; i++) {
                if (dataset[d[l] * COL_NUM + discrete_axises[i]] !=
                    dataset[d[r + 1] * COL_NUM + discrete_axises[i]]) {
                    goto get_r_end;
                }
            }
            r++;
        }
    get_r_end:

        /* build */
        /*
        -1:     8.1 GB  5.50 s  104 s   9.9e-11
        -2:     5.0 GB  3.46 s  82.4 s  1.07e-10
        -3:     3.0 GB  2.10 s  70.1 s  1.17e-10
        -4:     1.8 GB  1.32 s  62.4 s  2e-10
        -5:     1.1 GB  0.84 s  57.1 s  6e-10
        -6:     675 MB  0.57 s  54.0 s  5e-9
        -6(Of)          0.60 s  52.7 s  5e-9
        -8:     347 MB  0.34 s  50.2 s  1e-8
        -10:    263 MB  0.28 s  47.2 s  1e-7
        -12:    244 MB  0.26 s  45.5 s  5e-6
        -15:    240 MB  0.27 s  42.3 s  3e-5
        */
        max_depth = std::max(1, int(log2(r - l + 1) + delta_depth));
        // printf("%s %d %d\n", model_name.c_str(), l, r);
        Node* root = buildKDTree(tmp_data, l, r, 0);
        saveIDKDTree(root, idx);
        clearKDTree(root);
        l = r + 1;
    }

    fprintf(model_list_file, "%s\n", model_name.c_str());

#ifdef INFO
    printf("model_name=%s\nmodel_path=%s\n", model_name.c_str(), model_path.c_str());
    printKDTree(root, 0);
#endif

    fclose(model_file);
    fclose(model_list_file);
    delete[] tmp_data;
    delete[] d;
}

void clear_models() {
    for (auto& model_name : model_list) {
        clear_model(model_name);
    }
    model_list.clear();
}

extern "C" void clear() {
    clearData();
    clearans();
    clear_models();
}

extern "C" Answer* aqpQuery(Operation* ops,
                            int ops_size,
                            Predication* preds,
                            int preds_size,
                            COL_T groupBy_col,
                            MODE mode) {
    return aqp_group_query(preds, preds_size, ops, ops_size, groupBy_col, mode);
}