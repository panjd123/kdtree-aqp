#pragma once
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <utility>
#include <vector>

const int DATA_DIM = 7;
const int COL_NUM = 12;
using COL_T = int;
enum OP {
    COUNT,
    SUM,
    AVG,
};
using OP_T = OP;
using INT_T = int;
using FLOAT_T = float;
// using DATA_T = FLOAT_T[DATA_DIM];
using DATA_T = std::array<FLOAT_T, DATA_DIM>;
using BOUND_T = float[DATA_DIM][2];
using COL_VALUE_T = std::vector<std::pair<int, int>>;

#define IS_LEAF(u) ((u)->lchild == nullptr && (u)->rchild == nullptr)
#define IS_CONTINUED(c) ((c <= 6))
#define IS_DISCRETE(c) ((c >= 7))
#define GB (1024ull * 1024 * 1024)
#define MEM_LIMIT (10 * GB)

enum MODE {
    PERFORMANCE,
    MEMORY,
};

struct Operation {
    OP_T op;
    COL_T col;
};

struct Predication {
    COL_T col;
    FLOAT_T lb;
    FLOAT_T ub;
};

struct GroupAnswer {
    INT_T id;
    FLOAT_T value;
};

struct Answer {
    GroupAnswer* group_ans;
    int size;
};

struct AnswerBatch {
    Answer* ans;
    int size;
};

struct Node {
    struct Node* lchild;
    struct Node* rchild;
    int count;
    FLOAT_T sum[DATA_DIM];
    BOUND_T bound;
};

/* KD tree module */

// Building a tree based on the basic parameters of the KD tree
extern "C" void build(INT_T* col, int size, int delta_depth, float _build_k);

// Calculate the cross ratio for approximate calculations (consider all dimensions)
double data_cross_ratio(const BOUND_T& bound_in, const BOUND_T& bound_out);

// Whether it includes (only considering the KD tree segmentation dimension)
int kd_contain(const BOUND_T& bound_in, const BOUND_T& bound_out);

// Whether to intersect (considering only the KD tree segmentation dimension)
int kd_cross(const BOUND_T& bound_in, const BOUND_T& bound_out);

// Use the data in the range of [L, R) to establish a KD tree
Node* buildKDTree(DATA_T* data, int l, int r, int depth);

// Write tree to hard disk
void saveKDTree(Node* u);

// Write a ID before writing a tree
void saveIDKDTree(Node* u, int id);

// Load a tree to memory
Node* loadKDTree();

// Release tree memory
void clearKDTree(Node* u);

void testKDTree(Node* u, int depth);

// Print tree to the screen
void printKDTree(Node* u, int depth);

/* Query module */

// Get the model name corresponding to the column
std::string get_model_name(COL_VALUE_T& col_value);

// Get the model path corresponding to model_name
std::string get_model_path(std::string model_name);

// Get the root node of the KD tree corresponding to the column
Node* get_root(COL_VALUE_T col_value);

// Recursive query
void _queryRange(Node* u, const BOUND_T& bound, FLOAT_T* sum, double* count);

// Initialization query and recursive query
void queryRange(Node* root, BOUND_T& bound, FLOAT_T* sum, double& count);

// Get the answer to the query
extern "C" Answer* aqpQuery(Operation* ops, int, Predication* preds, int, COL_T, MODE);

// Release the memory of the last answer
void clearans();

/* Initialization module */

void load_col_type();

extern "C" void init(const char* dir);

/* Load dataset module */

extern "C" void loadData(FLOAT_T* _data, int n);

void clearData();

/* Loading model */

void load_model(const std::string model_name);

void clear_model(const std::string model_name);

extern "C" void load_models();

/* Destructive module */

void clear_models();

extern "C" void clear();