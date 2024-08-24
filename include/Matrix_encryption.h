/**
* @author: WTY
* @date: 2024/8/22
* @description: Definition of constants, operations, and header files
*/

#ifndef MATRIX_ENCRYPTION_H
#define MATRIX_ENCRYPTION_H

#include <iostream>
#include <ctime>
#include <chrono>
#include<vector>
#include <Eigen/Dense>
#include <Eigen/Dense>
#include<random>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;


/**
 * @Method: generateInvertibleMatrix
 * @Description: 生成一个N*N的随机可逆矩阵
 * @param int N: 矩阵的维度
 * @return MatrixXd: 返回生成的随机可逆矩阵
 */
MatrixXd generateInvertibleMatrix(int N);

/**
 * @Method: calculateInverseMatrix
 * @Description: 计算矩阵的逆矩阵
 * @param const MatrixXd& matrix: 输入矩阵
 * @return MatrixXd: 返回矩阵的逆矩阵
 */
MatrixXd calculateInverseMatrix(const MatrixXd& matrix);

/**
 * @Method: generateRandomDouble
 * @Description: 生成一个1到100之间的随机浮点数
 * @return double: 返回生成的随机浮点数
 */
double generateRandomDouble();



#endif //MATRIX_ENCRYPTION_H
