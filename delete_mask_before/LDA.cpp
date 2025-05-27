#include <vector>
#include <Eigen/Dense>
#include <algorithm>
#include <stdexcept>
#include <numeric> // for std::iota
using namespace std;

// 对 Q 进行 LDA 降维并返回投影矩阵
Eigen::MatrixXd lda_reduce_with_projection(const vector<vector<double>>& data,
                                           const vector<int>& labels,
                                           int target_dim,
                                           vector<vector<double>>& reduced_data) {
    int m = data.size();
    int n = data[0].size();

    // 初始化 reduced_data
    reduced_data.clear();
    reduced_data.resize(m, vector<double>(target_dim));

    // 将数据转换为 Eigen 矩阵
    Eigen::MatrixXd mat(m, n);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            mat(i, j) = data[i][j];
        }
    }

    // 计算类内散度矩阵 Sw 和类间散度矩阵 Sb
    Eigen::MatrixXd Sw = Eigen::MatrixXd::Zero(n, n);
    Eigen::MatrixXd Sb = Eigen::MatrixXd::Zero(n, n);

    // 计算总体均值
    Eigen::VectorXd mean_total = mat.colwise().mean();

    // 类别数量
    int num_classes = *max_element(labels.begin(), labels.end()) + 1;

    vector<int> class_counts(num_classes, 0);
    vector<Eigen::VectorXd> class_means(num_classes, Eigen::VectorXd::Zero(n));

    // 计算每个类别的均值和样本数量
    for (int i = 0; i < m; ++i) {
        int label = labels[i];
        class_means[label] += mat.row(i);
        class_counts[label]++;
    }

    for (int c = 0; c < num_classes; ++c) {
        class_means[c] /= class_counts[c];
    }

    // 计算 Sw 和 Sb
    for (int i = 0; i < m; ++i) {
        int label = labels[i];
        Eigen::VectorXd diff_within = mat.row(i) - class_means[label];
        Sw += diff_within.transpose() * diff_within;

        Eigen::VectorXd diff_between = class_means[label] - mean_total;
        Sb += class_counts[label] * diff_between.transpose() * diff_between;
    }

    // 计算广义特征值和特征向量
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eig(Sb, Sw);
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

// 使用 Q 的 LDA 投影矩阵对 A 进行降维
void reduce_A_with_Q_lda_projection(const vector<vector<double>>& A,
                                    const Eigen::MatrixXd& projection_matrix,
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

    // 使用投影矩阵
    Eigen::MatrixXd reduced_A = mat_A * projection_matrix;

    // 填充结果
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < projection_matrix.cols(); ++j) {
            A_reduced[i][j] = reduced_A(i, j);
        }
    }
}
