/**
* @author: WTY
* @date: 2024/8/22
* @description: Definition of constants, operations, and header files
*/

#include "Skyline_search.h"
#include "Matrix_encryption.h"

// 存储kd树所有叶子节点
vector<KDNode *> nodes;

// 加密矩阵，用于加密数据点
MatrixXd encryptMatrix1;

// 加密矩阵，用于查询边界
MatrixXd encryptMatrix2;

// 保存随机数r11,r12,r13
double r11, r12, r13;

// 保存随机数r21,r22,r23
double r21, r22, r23;

// 保存kd树根节点
// KDNode* root;

// 保存每棵kd树的树根
vector<KDNode *> kdTrees;

// 定义一个全局互斥锁
mutex nodes_mutex;

/**
 * @Method: readDataFromFile
 * @Description: 读取文件中的doubles，并返回一个vector<vector<double>>类型的数据
 * @param char* filename 文件名
 * @return vector<vector<double>> doubles数据
 */
vector<vector<double> > readDataFromFile(const char *filename) {
    vector<vector<double> > data_list;
    ifstream infile(filename);

    if (!infile.is_open()) {
        cerr << "Error opening file" << endl;
        return data_list;
    }

    string line;
    while (getline(infile, line)) {
        vector<double> row;
        stringstream ss(line);
        double number;

        while (ss >> number) {
            row.push_back(number);
        }

        data_list.push_back(row);
    }

    infile.close();
    return data_list;
}


/**
 * @Method: buildKDTree
 * @Description: 递归构建kd树
 * @param vector<vector<double>>& points 待构建的doubles数据
 * @param int depth 当前深度
 * @return KDNode* kd树根节点
 */
KDNode *buildKDTree(vector<vector<double> > &points, int depth) {
    // 如果数据集为空，返回空指针
    if (points.empty()) {
        return nullptr;
    }

    // 如果数据集只有一个点，返回叶子节点
    if (points.size() == 1) {
        KDNode *node = new KDNode(points[0]);

        // 使用互斥锁保护对 nodes 数组的写操作
        std::lock_guard<std::mutex> lock(nodes_mutex);
        nodes.push_back(node); // 将叶子节点存入 nodes 数组

        return node;
    }

    // 计算切割维度
    int dim = depth % points[0].size();

    // 按切割维度对数据集进行排序
    sort(points.begin(), points.end(), [dim](const vector<double> &a, const vector<double> &b) {
        return a[dim] < b[dim];
    });

    // 创建内部节点
    KDNode *node;

    // 如果只有两个点，取平均值作为切割点
    if (points.size() == 2) {
        double median_data = (points[0][dim] + points[1][dim]) / 2;
        // 创建内部节点
        node = new KDNode(dim, median_data);

        // 将中间节点加密
        vector<double> t(points[0].size() + 4, 0);
        t[dim] = r11;
        t[points[0].size()] = -median_data * r11;
        t[points[0].size() + 1] = r12 * r11;
        t[points[0].size() + 2] = -r12 * r11;
        t[points[0].size() + 3] = r13;
        VectorXd v = Eigen::Map<VectorXd>(t.data(), t.size()); // 将vector<double>转换为Eigen::VectorXd
        // 加密
        v = encryptMatrix2.transpose() * v;
        node->encrypted_value = v;

        // 将数据分为两部分并递归构建左右子树
        vector<vector<double> > leftPoints(points.begin(), points.begin() + 1);
        vector<vector<double> > rightPoints(points.begin() + 1, points.end());

        // 递归构建左右子树
        node->left = buildKDTree(leftPoints, depth + 1);
        node->right = buildKDTree(rightPoints, depth + 1);
    } else {
        // 取中位数作为切割点
        int median = points.size() / 2;

        // 创建内部节点
        node = new KDNode(dim, points[median][dim]);

        // 将中间节点加密
        vector<double> t(points[0].size() + 4, 0);
        t[dim] = r11;
        t[points[0].size()] = -points[median][dim] * r11;
        t[points[0].size() + 1] = r12 * r11;
        t[points[0].size() + 2] = -r12 * r11;
        t[points[0].size() + 3] = r13;
        VectorXd v = Eigen::Map<VectorXd>(t.data(), t.size()); // 将vector<double>转换为Eigen::VectorXd
        // 加密
        v = encryptMatrix2.transpose() * v;
        node->encrypted_value = v;

        // 将数据分为两部分并递归构建左右子树
        vector<vector<double> > leftPoints(points.begin(), points.begin() + median + 1);
        vector<vector<double> > rightPoints(points.begin() + median + 1, points.end());

        // 递归构建左右子树
        node->left = buildKDTree(leftPoints, depth + 1);
        node->right = buildKDTree(rightPoints, depth + 1);
    }

    return node;
}

/**
 * @Method: dealData
 * @Description: 对数据集进行预计算
 * @param char* fileString 读取数据集的地址
 * @return 状态码，1：成功；0：失败
 */
int dealData(char *fileString) {
    auto start_time = chrono::high_resolution_clock::now();
    // 读取数据集
    vector<vector<double> > points = readDataFromFile(fileString);

    // 获取结束时间点
    auto end_time = chrono::high_resolution_clock::now();
    // 计算时间间隔
    chrono::duration<double, milli> total_duration = end_time - start_time;
    // 输出时间间隔
    printf("数据读取的时间是：%f 毫秒\n", total_duration.count());
    fflush(stdout);

    // 检查数据集是否为空
    if (points.empty()) {
        return 0; // 数据集为空，返回错误码
    }

    nodes.reserve(points.size()); // 预先分配空间

    // 生成加密矩阵
    encryptMatrix1 = generateInvertibleMatrix(6);

    // 生成三个随机值r11,r12,r13，其中r11 > r13 > 0
    r11 = generateRandomDouble();
    r12 = generateRandomDouble();
    r13 = generateRandomDouble();
    while (r11 <= r13) {
        r11 = generateRandomDouble();
        r13 = generateRandomDouble();
    }

    // 生成加密矩阵
    encryptMatrix2 = generateInvertibleMatrix(points[0].size() + 4);

    // 生成三个随机值r21,r22,r23，其中r21 > r23 > 0
    r21 = generateRandomDouble();
    r22 = generateRandomDouble();
    r23 = generateRandomDouble();
    while (r21 <= r23) {
        r21 = generateRandomDouble();
        r23 = generateRandomDouble();
    }

    start_time = chrono::high_resolution_clock::now();

    // 将数据集均分为N份
    int numPoints = points.size();
    int pointsPerGroup = numPoints / N;
    int remainder = numPoints % N;

    vector<vector<vector<double> > > dividedPoints(N); // 创建一个大小为N的二维向量，用于存储N个数据子集

    int startIndex = 0;
    for (int i = 0; i < N; i++) {
        int endIndex = startIndex + pointsPerGroup;
        if (i < remainder) {
            endIndex++;
        }
        for (int j = startIndex; j < endIndex; j++) {
            dividedPoints[i].push_back(points[j]);
        }
        startIndex = endIndex;
    }

    // 保存每棵kd树的树根
    kdTrees.resize(N);

    // 创建N个线程，每个线程处理一个数据子集并构建一棵kd树
    vector<thread> threads(N);
    for (int i = 0; i < N; i++) {
        threads[i] = thread([i, &dividedPoints]() {
            kdTrees[i] = buildKDTree(dividedPoints[i], 0);
        });
    }

    // 等待所有线程完成
    for (int i = 0; i < N; i++) {
        if (threads[i].joinable()) {
            threads[i].join();
        }
    }

    // 输出所有叶子结点的个数
    cout << "所有叶子结点的个数：" << nodes.size() << endl;

    // 构建kd树
    // root = buildKDTree(points, 0);

    // 获取结束时间点
    end_time = chrono::high_resolution_clock::now();
    // 计算时间间隔
    total_duration = end_time - start_time;
    // 输出时间间隔
    printf("构建kd树的时间是：%f 毫秒\n", total_duration.count());
    fflush(stdout);


    start_time = chrono::high_resolution_clock::now();

    // // 对kd树的所有叶子节点进行加密，并保存到encrypted_point中
    // for (int i = 0; i < nodes.size(); i++) { // 对每个叶子节点进行加密
    //    for (int j = 0; j < nodes[i]->point.size(); j++) { // 对每个叶子节点的每个维度进行加密
    //        vector<double> t(6); // 创建一个大小为6的向量t
    //        t[0] = nodes[i]->point[j] * nodes[i]->point[j] * r11;
    //        t[1] = nodes[i]->point[j] * r11;
    //        t[2] = r11;
    //        t[3] = r12 * r11;
    //        t[4] = -r12 * r11;
    //        t[5] = r13;
    //        VectorXd v = Eigen::Map<VectorXd>(t.data(), t.size()); // 将vector<double>转换为Eigen::VectorXd
    //        // 加密
    //        v = encryptMatrix1.transpose() * v;
    //        nodes[i]->encrypted_point.push_back(v); // 将加密后的向量存入encrypted_point
    //    }
    // }

    // 对kd树的所有叶子节点进行加密，并保存到encrypted_point中
    vector<thread> threads2(N); // 创建N个线程
    int nodesPerThread = nodes.size() / N; // 每个线程需要处理的节点数
    int remainder2 = nodes.size() % N; // 处理不均匀分配的剩余节点数

    // 为每个线程分配工作
    for (int t = 0; t < N; t++) {
        int startIdx = t * nodesPerThread + min(t, remainder2);
        int endIdx = startIdx + nodesPerThread + (t < remainder2 ? 1 : 0);

        threads2[t] = thread([startIdx, endIdx]() {
            for (int i = startIdx; i < endIdx; i++) {
                // 对分配的节点进行加密
                for (int j = 0; j < nodes[i]->point.size(); j++) {
                    // 对每个叶子节点的每个维度进行加密
                    vector<double> t(6); // 创建一个大小为6的向量t
                    t[0] = nodes[i]->point[j] * nodes[i]->point[j] * r11;
                    t[1] = nodes[i]->point[j] * r11;
                    t[2] = r11;
                    t[3] = r12 * r11;
                    t[4] = -r12 * r11;
                    t[5] = r13;
                    VectorXd v = Eigen::Map<VectorXd>(t.data(), t.size()); // 将vector<double>转换为Eigen::VectorXd
                    // 加密
                    v = encryptMatrix1.transpose() * v;
                    nodes[i]->encrypted_point.push_back(v); // 将加密后的向量存入encrypted_point
                }
            }
        });
    }

    // 等待所有线程完成
    for (int t = 0; t < N; t++) {
        if (threads2[t].joinable()) {
            threads2[t].join();
        }
    }


    // 获取结束时间点
    end_time = chrono::high_resolution_clock::now();
    // 计算时间间隔
    total_duration = end_time - start_time;
    // 输出时间间隔
    printf("加密kd树的时间是：%f 毫秒\n", total_duration.count());
    fflush(stdout);
    printf("--------------------------------------------------------------\n");
    return 1;
}

/**
 * @Method: kd_search
 * @Description: 查询kd树
 * @param KDNode* node kd树的根节点
 * @param vector<VectorXd>& encrypted_query_data 查询请求的加密数据
 * @param VectorXd& lower_bound_vector 查询请求的下界
 * @param VectorXd& upper_bound_vector 查询请求的上界
 * @param vector<vector<VectorXd>>& res 查询结果
 */
void kd_search(KDNode *node, vector<VectorXd> &encrypted_query_data, VectorXd &lower_bound_vector,
               VectorXd &upper_bound_vector, vector<vector<VectorXd> > &res) {
    // 判断当前节点是否为空
    if (node == nullptr) {
        return;
    }
    // 判断当前节点是否为叶子节点
    if (node->left == nullptr && node->right == nullptr) {
        // 判断每一维的加密数据是否在当前节点的加密数据范围内
        for (int i = 0; i < encrypted_query_data.size(); i++) {
            if (node->encrypted_point[i].dot(encrypted_query_data[i]) > 0) {
                return;
            }
        }
        res.push_back(node->encrypted_point); // 将当前节点的加密数据存入结果
        return;
    }

    if (node->encrypted_value.dot(lower_bound_vector) <= 0) {
        kd_search(node->left, encrypted_query_data, lower_bound_vector, upper_bound_vector, res); // 递归查询左子树
    }

    if (node->encrypted_value.dot(upper_bound_vector) > 0) {
        kd_search(node->right, encrypted_query_data, lower_bound_vector, upper_bound_vector, res); // 递归查询右子树
    }
}

/**
 * @Method: range_search
 * @Description: 发起查询请求，并返回查询结果
 * @param char* fileString 读取数据的地址
 * @param char* resultFilePath 输出数据的地址
 * @return 状态码，1：成功；0：失败
 */
int skyline_search(char *fileString, char *resultFilePath) {
    // 读取查询请求(假设数据维度是d维，前d行是范围查询请求，每行2个数据；最后一行是skyline查询请求，含d个数据，空格分隔)
    vector<vector<double> > query_data = readDataFromFile(fileString);

    // 检查查询请求是否为空
    if (query_data.empty()) {
        return 0; // 查询请求为空，返回错误码
    }

    // 逆矩阵，用于解密
    MatrixXd encryptMatrixInverse2 = calculateInverseMatrix(encryptMatrix2);

    // 保存下边界
    vector<double> lower_bound(query_data.size() - 1 + 4);
    // 保存上边界
    vector<double> upper_bound(query_data.size() - 1 + 4);
    for (int i = 0; i < query_data.size() - 1; i++) {
        lower_bound[i] = query_data[i][0] * r21; // 下边界
        upper_bound[i] = query_data[i][1] * r21; // 上边界
    }
    lower_bound[query_data.size() - 1] = r21;
    upper_bound[query_data.size() - 1] = r21;
    lower_bound[query_data.size() - 1 + 1] = r22 * r21;
    upper_bound[query_data.size() - 1 + 1] = r22 * r21;
    lower_bound[query_data.size() - 1 + 2] = r22 * r21;
    upper_bound[query_data.size() - 1 + 2] = r22 * r21;
    lower_bound[query_data.size() - 1 + 3] = -r23;
    upper_bound[query_data.size() - 1 + 3] = -r23;

    VectorXd lower_bound_vector = Eigen::Map<VectorXd>(lower_bound.data(), lower_bound.size());
    VectorXd upper_bound_vector = Eigen::Map<VectorXd>(upper_bound.data(), upper_bound.size());
    lower_bound_vector = encryptMatrixInverse2 * lower_bound_vector;
    upper_bound_vector = encryptMatrixInverse2 * upper_bound_vector;


    // 逆矩阵，用于解密
    MatrixXd encryptMatrixInverse1 = calculateInverseMatrix(encryptMatrix1);

    // 加密的查询范围
    vector<VectorXd> encrypted_query_data;

    for (int i = 0; i < query_data.size() - 1; i++) {
        // 对每个查询范围进行加密
        vector<double> t(6); // 创建一个大小为6的向量t
        t[0] = r21;
        t[1] = -(query_data[i][0] + query_data[i][1]) * r21;
        t[2] = query_data[i][0] * query_data[i][1] * r21;
        t[3] = r22 * r21;
        t[4] = r22 * r21;
        t[5] = -r23;
        VectorXd v = Eigen::Map<VectorXd>(t.data(), t.size()); // 将vector<double>转换为Eigen::VectorXd
        // 加密
        v = encryptMatrixInverse1 * v;
        encrypted_query_data.push_back(v);
    }

    auto start_time = chrono::high_resolution_clock::now();

    // 先进行范围查询,res用于存储范围查询的结果
    vector<vector<VectorXd> > res;
    // kd_search(root, encrypted_query_data, lower_bound_vector, upper_bound_vector, res);
    // 对N颗kd树进行范围查询
    for (int i = 0; i < N; i++) {
        kd_search(kdTrees[i], encrypted_query_data, lower_bound_vector, upper_bound_vector, res);
    }

    cout << "范围查询的结果数量为：" << res.size() << endl;

    // 获取结束时间点
    auto end_time = chrono::high_resolution_clock::now();
    // 计算时间间隔
    chrono::duration<double, milli> total_duration = end_time - start_time;
    // 输出时间间隔
    printf("范围查询的时间是：%f 毫秒\n", total_duration.count());
    fflush(stdout);


    // 对skyline查询请求进行加密
    vector<VectorXd> encrypted_skyline_query;
    for (int i = 0; i < query_data[query_data.size() - 1].size(); i++) {
        // 对每一维的数据加密
        vector<double> t(6); // 创建一个大小为6的向量t
        t[0] = 1;
        t[1] = -2 * query_data[query_data.size() - 1][i];
        t[2] = pow(query_data[query_data.size() - 1][i], 2);
        t[3] = 0;
        t[4] = 0;
        t[5] = 0;
        VectorXd v = Eigen::Map<VectorXd>(t.data(), t.size()); // 将vector<double>转换为Eigen::VectorXd
        // 加密
        v = encryptMatrixInverse1 * v;
        encrypted_skyline_query.push_back(v);
    }

    start_time = chrono::high_resolution_clock::now();

    // 后进行范围查询,ans用于存储skyline查询的结果
    vector<bool> flag(res.size(), true); // 用于标记res中哪些点满足skyline查询

    for (int i = 0; i < res.size(); i++) {
        bool flag1 = true; // 判断是否每一维的距离都大于等于
        bool flag2 = false; // 判断是否至少有一维的距离大于
        for (int j = 0; j < res.size(); j++) {
            if (i != j) {
                for (int k = 0; k < res[i].size(); k++) {
                    if (res[i][k].dot(encrypted_skyline_query[k]) < res[j][k].dot(encrypted_skyline_query[k])) {
                        flag1 = false;
                        break;
                    } else if (res[i][k].dot(encrypted_skyline_query[k]) > res[j][k].dot(encrypted_skyline_query[k])) {
                        flag2 = true;
                    }
                }
            }
            if (flag1 && flag2) {
                flag[i] = false; // 不满足skyline查询
            }
        }
    }

    // 获取结束时间点
    end_time = chrono::high_resolution_clock::now();
    // 计算时间间隔
    total_duration = end_time - start_time;
    // 输出时间间隔
    printf("skyline查询的时间是：%f 毫秒\n", total_duration.count());
    fflush(stdout);

    printf("--------------------------------------------------------------\n");

    // 将ans内的数据写入文件
    ofstream resultFile(resultFilePath);
    if (resultFile.is_open()) {
        for (int i = 0; i < res.size(); i++) {
            if (flag[i]) {
                for (int j = 0; j < res[i].size(); j++) {
                    // 将ans[i][j]解密
                    VectorXd decrypted_data = res[i][j].transpose() * encryptMatrixInverse1;
                    resultFile << decrypted_data[1] / r11 << " "; // 写入文件
                }
                resultFile << endl;
            }
        }
        resultFile.close(); // 关闭文件
    } else {
        cerr << "Unable to open file " << resultFilePath << endl;
        return 0;
    }

    // // 将res内的数据写入文件
    // ofstream resultFile(resultFilePath);
    // if (resultFile.is_open()) {
    //     for (int i = 0; i < res.size(); i++) {
    //         for (int j = 0; j < res[i].size(); j++) {
    //             // 将res[i][j]解密
    //             VectorXd decrypted_data = res[i][j].transpose() * encryptMatrixInverse1;
    //             resultFile << decrypted_data[1] / r11 << " "; // 写入文件
    //         }
    //         resultFile << endl;
    //     }
    //     resultFile.close(); // 关闭文件
    // } else {
    //     cerr << "Unable to open file " << resultFilePath << endl;
    //     return 0;
    // }
    return 1;
}
