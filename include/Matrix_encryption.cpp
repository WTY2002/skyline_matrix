/**
* @author: WTY
* @date: 2024/8/22
* @description: Definition of constants, operations, and header files
*/

#include "Matrix_encryption.h"



/**
 * @Method: generateInvertibleMatrix
 * @Description: 生成一个N*N的随机可逆矩阵
 * @param int N: 矩阵的维度
 * @return MatrixXd: 返回生成的随机可逆矩阵
 */
MatrixXd generateInvertibleMatrix(int N) {
    MatrixXd matrix;
    do {
        // 生成随机矩阵
        matrix = MatrixXd::Random(N, N);
    } while (matrix.determinant() == 0);  // 检查矩阵是否可逆
    return matrix;
}

/**
 * @Method: calculateInverseMatrix
 * @Description: 计算矩阵的逆矩阵
 * @param const MatrixXd& matrix: 输入矩阵
 * @return MatrixXd: 返回矩阵的逆矩阵
 */
MatrixXd calculateInverseMatrix(const MatrixXd& matrix) {
    return matrix.inverse();
}

/**
 * @Method: generateRandomDouble
 * @Description: 生成一个1到100之间的随机浮点数
 * @return double: 返回生成的随机浮点数
 */
double generateRandomDouble() {
    // 创建随机数引擎，使用随机设备作为种子
    random_device rd;
    mt19937 generator(rd());
    
    // 定义均匀分布的范围
    uniform_real_distribution<double> distribution(1, 100);
    
    // 生成并返回一个随机数
    return distribution(generator);
}