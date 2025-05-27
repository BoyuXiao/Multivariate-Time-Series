#include <vector>
#include <Eigen/Dense> // 使用 Eigen 库进行 PCA 计算
#include "funcs.h"

// 对 Q 进行降维并返回投影矩阵
Eigen::MatrixXd pca_reduce_with_projection(const std::vector<std::vector<double>>& data, int target_dim, std::vector<std::vector<double>>& reduced_data) {
    int m = data.size();    // 行数
    int n = data[0].size(); // 列数（原始维度）

    // 将数据转换为 Eigen 矩阵
    Eigen::MatrixXd mat(m, n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            mat(i, j) = data[i][j];
        }
    }

    // 中心化数据（减去均值）
    Eigen::VectorXd mean = mat.colwise().mean();
    mat.rowwise() -= mean.transpose();

    // 计算协方差矩阵
    Eigen::MatrixXd cov = (mat.adjoint() * mat) / double(m - 1);

    // 计算特征值和特征向量
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
    Eigen::MatrixXd eigenvectors = eig.eigenvectors(); // 特征向量
    Eigen::VectorXd eigenvalues = eig.eigenvalues();  // 特征值

    // 选择前 target_dim 个主成分
    Eigen::MatrixXd projection_matrix = eigenvectors.rightCols(target_dim);

    // 投影数据到低维空间
    Eigen::MatrixXd reduced_mat = mat * projection_matrix;

    // 将结果填充到已初始化的 reduced_data
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < target_dim; ++j) {
            reduced_data[i][j] = reduced_mat(i, j);
        }
    }
   
    return projection_matrix;
}
// 使用 Q 的投影矩阵对 A 进行降维
void reduce_A_with_Q_projection(const std::vector<std::vector<double>>& A,
                                const Eigen::MatrixXd& projection_matrix,
                                std::vector<std::vector<double>>& A_reduced) {
    int m = A.size();    // 行数
    int n = A[0].size(); // 列数（原始维度）

    // 将 A 转换为 Eigen 矩阵
    Eigen::MatrixXd mat_A(m, n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            mat_A(i, j) = A[i][j];
        }
    }

    // 中心化 A（减去均值）
    Eigen::VectorXd mean_A = mat_A.colwise().mean();
    mat_A.rowwise() -= mean_A.transpose();

    // 使用 Q 的投影矩阵对 A 进行降维
    Eigen::MatrixXd reduced_A = mat_A * projection_matrix;

    // 将结果填充到已初始化的 A_reduced
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < projection_matrix.cols(); ++j) {
            A_reduced[i][j] = reduced_A(i, j);
        }
    }
}