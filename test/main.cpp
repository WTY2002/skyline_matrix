#include "Matrix_encryption.h"
#include "Skyline_search.h"

// 构建kd树
KDNode* buildKDTree_test(vector<vector<double>>& points, int depth) {
    // 如果数据集为空，返回空指针
    if (points.empty()) {
        return nullptr;
    }

    // 如果数据集只有一个点，返回叶子节点
    if (points.size() == 1) {
        KDNode* node = new KDNode(points[0]);
        return node;
    }

    // 计算切割维度
    int dim = depth % points[0].size();

    // 按切割维度对数据集进行排序
    sort(points.begin(), points.end(), [dim](const vector<double>& a, const vector<double>& b) {
        return a[dim] < b[dim];
    });

    // 创建内部节点
    KDNode* node;

    // 如果只有两个点，取平均值作为切割点
    if (points.size() == 2) {
        double median_data = (points[0][dim] + points[1][dim]) / 2;
        // 创建内部节点
        node = new KDNode(dim, median_data);

        // 将数据分为两部分并递归构建左右子树
        vector<vector<double>> leftPoints(points.begin(), points.begin() + 1);
        vector<vector<double>> rightPoints(points.begin() + 1, points.end());

        // 递归构建左右子树
        node->left = buildKDTree_test(leftPoints, depth + 1);
        node->right = buildKDTree_test(rightPoints, depth + 1);
    } else {
        // 取中位数作为切割点
        int median = points.size() / 2;

        // 创建内部节点
        node = new KDNode(dim, points[median][dim]);

        // 将数据分为两部分并递归构建左右子树
        vector<vector<double>> leftPoints(points.begin(), points.begin() + median + 1);
        vector<vector<double>> rightPoints(points.begin() + median + 1, points.end());

        // 递归构建左右子树
        node->left = buildKDTree_test(leftPoints, depth + 1);
        node->right = buildKDTree_test(rightPoints, depth + 1);
    }

    return node;
}

// 打印kd树的递归函数
void printKDTree(const KDNode* node, int depth) {
    if (node == nullptr) {
        return;
    }

    // 打印当前节点的缩进
    for (int i = 0; i < depth; ++i) {
        cout << "    ";
    }

    // 判断是叶子节点还是中间节点
    if (node->left == nullptr && node->right == nullptr) {
        // 叶子节点，打印向量
        cout << "Leaf: (";
        for (size_t i = 0; i < node->point.size(); ++i) {
            cout << node->point[i];
            if (i < node->point.size() - 1) {
                cout << ", ";
            }
        }
        cout << ")" << endl;
    } else {
        // 中间节点，打印切割维度和切割值
        cout << "Node: dimension = " << node->dimension << ", value = " << node->value << endl;
    }

    // 递归打印左右子树
    printKDTree(node->left, depth + 1);
    printKDTree(node->right, depth + 1);
}

// 求kd树的高度
int getHeight(const KDNode* node) {
    if (node == nullptr) {
        return 0;
    }

    int leftHeight = getHeight(node->left);
    int rightHeight = getHeight(node->right);

    return max(leftHeight, rightHeight) + 1;
}

void test() {
    // 示例数据集，k=3维
    vector<vector<double>> points = {
        {2.1, 3.2, 1.5},
        {5.4, 2.8, 9.1},
        {1.5, 4.5, 2.8},
        {3.9, 3.6, 7.5},
        {7.2, 1.1, 6.3}
    };

    // 构建kd树
    KDNode* root = buildKDTree_test(points, 0);

    // 打印kd树
    cout << "KD-Tree structure:" << endl;
    printKDTree(root, 0);

    // 获取kd树的高度
    int height = getHeight(root);
    cout << "Height of KD-Tree: " << height << endl;


    // 释放kd树
    delete root;
}


// 主函数示例
int main() {

    // test();

    char* fileString = "/root/wty/data.txt";
    char* query = "/root/wty/query.txt";
    char* res = "/root/wty/result.txt";


    auto start_time = chrono::high_resolution_clock::now();

    dealData(fileString); // 预处理数据

    // 获取结束时间点
    auto end_time = chrono::high_resolution_clock::now();
    // 计算时间间隔
    chrono::duration<double, milli> total_duration = end_time - start_time;
    // 输出时间间隔
    printf("======>>>数据加密外包的时间是：%f 毫秒\n", total_duration.count());
    fflush(stdout);

    printf("--------------------------------------------------------------\n");

    // int height = getHeight(root); // 获取kd树的高度
    // cout << "Height of KD-Tree: " << height << endl;

    auto start_time2 = chrono::high_resolution_clock::now();

    skyline_search(query, res); // 查询

    // 获取结束时间点
    auto end_time2 = chrono::high_resolution_clock::now();
    // 计算时间间隔
    chrono::duration<double, milli> total_duration2 = end_time2 - start_time2;
    // 输出时间间隔
    printf("======>>>查询的总时间是：%f 毫秒\n", total_duration2.count());
    fflush(stdout);

    printf("--------------------------------------------------------------\n");


    // 释放kd树
    // delete root;
    for (int i = 0; i < kdTrees.size(); i++) {
        delete kdTrees[i];
    }

    return 0;
}

