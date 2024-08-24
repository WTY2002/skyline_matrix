/**
* @author: WTY
* @date: 2024/8/22
* @description: Definition of constants, operations, and header files
*/

#ifndef RANGE_SEARCH_H
#define RANGE_SEARCH_H
#include "Matrix_encryption.h"
#include <queue>
#include <fstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <sys/mman.h>
#include <thread>
#include <mutex>
#include <atomic>
#include <cmath>

// 定义kd树节点结构体
struct KDNode {
    int dimension; // 切割维度
    double value;  // 切割值
    KDNode* left;  // 左子树
    KDNode* right; // 右子树
    vector<double> point; // 存储k维向量（仅在叶子节点中有效）
    vector<VectorXd> encrypted_point; // 存储加密后的k维向量
    VectorXd encrypted_value; // 存储加密后的中间节点

    // 构造函数，初始化为叶子节点
    KDNode(const vector<double>& pt) : dimension(-1), value(0.0), point(pt), left(nullptr), right(nullptr) {}

    // 构造函数，初始化为内部节点
    KDNode(int dim, double val) : dimension(dim), value(val), left(nullptr), right(nullptr) {}

    // 析构函数，递归释放子节点
    ~KDNode() {
        delete left;
        delete right;
    }
};

// 存储kd树所有叶子节点
extern vector<KDNode*> nodes;

// 加密矩阵，用于加密数据点
extern MatrixXd encryptMatrix1;

// 加密矩阵，用于查询边界
extern MatrixXd encryptMatrix2;

// 保存随机数r11,r12,r13
extern double r11, r12, r13;

// 保存随机数r21,r22,r23
extern double r21, r22, r23;

// 保存kd树根节点
// extern KDNode* root;

// 设置创建kd树的个数
const int N = 10;

// 保存每棵kd树的树根
extern vector<KDNode*> kdTrees;

/**
 * @Method: readDataFromFile
 * @Description: 读取文件中的doubles，并返回一个vector<vector<double>>类型的数据
 * @param char* filename 文件名
 * @return vector<vector<double>> doubles数据
 */
vector<vector<double>> readDataFromFile(const char* filename);

/**
 * @Method: buildKDTree
 * @Description: 递归构建kd树
 * @param vector<vector<double>>& points 待构建的doubles数据
 * @param int depth 当前深度
 * @return KDNode* kd树根节点
 */
KDNode* buildKDTree(vector<vector<double>>& points, int depth);

/**
 * @Method: dealData
 * @Description: 对数据集进行预计算
 * @param char* fileString 读取数据集的地址
 * @return 状态码，1：成功；0：失败
 */
int dealData(char* fileString);

/**
 * @Method: kd_search
 * @Description: 查询kd树
 * @param KDNode* node kd树的根节点
 * @param vector<VectorXd>& encrypted_query_data 查询请求的加密数据
 * @param VectorXd& lower_bound_vector 查询请求的下界
 * @param VectorXd& upper_bound_vector 查询请求的上界
 * @param vector<vector<VectorXd>>& res 查询结果
 */
void kd_search(KDNode* node, vector<VectorXd>& encrypted_query_data, VectorXd& lower_bound_vector, VectorXd& upper_bound_vector, vector<vector<VectorXd>>& res);

/**
 * @Method: range_search
 * @Description: 发起查询请求，并返回查询结果
 * @param char* fileString 读取数据的地址
 * @param char* resultFilePath 输出数据的地址
 * @return 状态码，1：成功；0：失败
 */
int skyline_search(char *fileString, char *resultFilePath);


#endif //RANGE_SEARCH_H
