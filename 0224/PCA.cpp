
#include <vector>
#include <Eigen/Dense>
#include <algorithm>
#include <stdexcept>
#include <numeric> // for std::iota
using namespace std;

// 对 Q 进行降维并返回投影矩阵
Eigen::MatrixXd pca_reduce_with_projection(const vector<vector<double>>& data, int target_dim, vector<vector<double>>& reduced_data, Eigen::VectorXd& saved_mean) {
    int m = data.size();
    int n = data[0].size();

    // 初始化 reduced_data
    reduced_data.resize(m, vector<double>(target_dim));

    // 将数据转换为 Eigen 矩阵
    Eigen::MatrixXd mat(m, n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            mat(i, j) = data[i][j];
        }
    }

    // 中心化数据
    Eigen::VectorXd mean = mat.colwise().mean();
    saved_mean = mean;
    mat.rowwise() -= mean.transpose();

    // 计算协方差矩阵
    Eigen::MatrixXd cov = (mat.adjoint() * mat) / (m - 1);

    // 计算特征值和特征向量
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
    Eigen::MatrixXd eigenvectors = eig.eigenvectors();
    Eigen::VectorXd eigenvalues = eig.eigenvalues();

    // 对特征值降序排序
    vector<int> indices(n);
    iota(indices.begin(), indices.end(), 0);
    sort(indices.begin(), indices.end(), [&eigenvalues](int a, int b) {
        return eigenvalues[a] > eigenvalues[b];
    });

    // 构建投影矩阵
    Eigen::MatrixXd projection_matrix(n, target_dim);
    for (int i = 0; i < target_dim; ++i) {
        projection_matrix.col(i) = eigenvectors.col(indices[i]);
    }

    // 投影并填充结果
    Eigen::MatrixXd reduced_mat = mat * projection_matrix;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < target_dim; ++j) {
            reduced_data[i][j] = reduced_mat(i, j);
        }
    }

    return projection_matrix;
}

// 使用 Q 的投影矩阵对 A 进行降维
void reduce_A_with_Q_projection(const vector<vector<double>>& A,
                                const Eigen::MatrixXd& projection_matrix,
                                const Eigen::VectorXd& mean_Q,
                                vector<vector<double>>& A_reduced) {

    int m = A.size();
    int n = A[0].size();
    
    // 初始化 A_reduced
    A_reduced.resize(m, vector<double>(projection_matrix.cols()));

    // 转换到 Eigen 矩阵
    Eigen::MatrixXd mat_A(m, n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            mat_A(i, j) = A[i][j];
        }
    }

    // 中心化并使用投影矩阵
    mat_A.rowwise() -= mean_Q.transpose();
    Eigen::MatrixXd reduced_A = mat_A * projection_matrix;

    // 填充结果
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < projection_matrix.cols(); ++j) {
            A_reduced[i][j] = reduced_A(i, j);
        }
    }
}