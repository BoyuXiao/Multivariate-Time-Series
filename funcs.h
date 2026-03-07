
#ifndef MULTI_DTW_FUNCS_H
#define MULTI_DTW_FUNCS_H
#include "def.h"
#include <vector>
#include <string>
#include <cmath>
#include "bucket.h"
#include "cell.h"
#include <Eigen/Dense>
using namespace std;
void LB_MV_T_hot_map(vector<vector<double>>& lowerBound,vector<vector<double>>& upperBound, vector<vector<double>>& A, int m, int K);
// double LB_MV_double_T(vector<vector<double>>& A,vector<vector<double>>& Q, int r);
double LB_MV_double_T(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r);
double LB_MV_double_IDEA_T_V2(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& A_t,vector<vector<double>>& Q_t,vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r);
double LB_MV_double_KIM_T(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r);
double LB_MV_double_T_V2(vector<vector<double>>& A,vector<vector<double>>& Q,vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r);
double LB_MV_double_IDEA_T(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r);
void lower_upper_lemire_index(vector<double>& Q, int m, int r, vector<pair<int,double>>& lowerBound, vector<pair<int,double>>& upperBound);
double LB_MV_double(vector<vector<double>>& A,vector<vector<double>>& Q,vector<double>&lb_record,vector<char>&mask,int r);
// double LB_GRAPH(vector<vector<double>>& A,vector<vector<double>>& Q,
//     vector<double>& U_mu_Q,vector<double>& L_mu_Q,
//     vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
//     vector<double>& Dbl, vector<double>& Dbr,
//     vector<double>& Dtl, vector<double>& Dtr,
//     vector<vector<vector<Cell>>>& T,
//     vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,vector<double>&lb_record,int r);
vector<double> getLB_oneQ_qbox(const vector<vector<double>>& X,
                             const vector<vector<vector<double>>>& others,
                             const vector<pair<vector<double>, vector<double>>>& qbounds);
void LB_P_MV_early(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r,double threshold,int &flag,double &lb);
void LB_P_MV_early(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r,double threshold,int &flag,double &lb,vector<double>&cb);
tuple<double, int, int, double, int> DTWDistanceWindowLB_Ordered_LBPC_(
    int K,
    double Q,
    const vector<vector<double>>& query,
    const vector<vector<vector<double>>>& references,
    int W);
double DTW_a(const vector<vector<double>>& s1,
             const vector<vector<double>>& s2,
             int windowSize,
             double bestdist);
tuple< int, int> DTWDistanceWindowLB_Ordered_LBTI(
    const vector<vector<double>>& query,
    const vector<vector<vector<double>>>& references,
    int W,double threshold, int P);
double LB_GRAPH(vector<vector<double>>& A,vector<vector<double>>& Q,vector<vector<double>>& A_T,vector<vector<double>>& Q_T,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<double>& mu_A,vector<double>& sigma_A,
    vector<double>& Dbl, vector<double>& Dbr,
    vector<double>& Dtl, vector<double>& Dtr,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,vector<double>&lb_record,int r) ;


double LB_MV_graph(vector<double>& Q, vector<vector<double>>& A, int K);
void upper_lemire(vector<double>& Q, int m, int r, vector<double>& upperBound);
void lower_lemire(vector<double>& Q, int m, int r, vector<double>& lowerBound);
void only_mu(vector<vector<double>>& A, vector<double>& mu_A);
void mu_sigma(vector<vector<double>>& A, vector<double>& mu_A, vector<double>& sigma_A);
void lower_upper_lemire(vector<double>& Q, int m, int r, vector<double>& lowerBound, vector<double>& upperBound);
void TLB(vector<vector<vector<double>>>& Test, vector<vector<vector<double>>>& Train,string name);
double euclideanDistance_2(vector<double>& point1, vector<double>& point2);
double dtw(vector<vector<double>>& A, vector<vector<double>>& B, int K, int m, int r);
double dtw(vector<vector<double>>& A, vector<vector<double>>& B, int K, int m, int r,double threshold);
double dtw(vector<vector<double>>& A, vector<vector<double>>& B, int K, int m, int r,double threshold,vector<double>&cb);
double LB_Keogh_2(vector<double>& A, vector<double>& Q, int m, int r);
double LB_Keogh_old(vector<double>& Q, vector<double>& A, int m, int r);
double LB_MV(vector<vector<double>>& Q, vector<vector<double>>& A, int K, int m, int r);
// double LB_MV_T(vector<vector<double>>& Q, vector<vector<double>>& A, int K, int m, int r);
double LB_MV_T(vector<vector<double>>& lowerBound,vector<vector<double>>& upperBound, vector<vector<double>>& A, int K, int m);
double LB_MC1(vector<vector<double>>& Q, vector<vector<double>>& A,vector<double>& mu_Q, int K, int m, int r);
double LB_MC2(vector<vector<double>>& Q, vector<vector<double>>& A,vector<double>& mu_Q, vector<double>& sigma_Q, int K, int m, int r);
vector<vector<double>> Preprocessing_MC3(vector<vector<double>>& Q,
    vector<double>& mu_Q,vector<double>& sigma_Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<vector<Bucket>>& T_mu,vector<vector<Bucket>>& T_sigma,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,int r);
vector<vector<double>> Preprocessing_MC3E(vector<vector<double>>& Q,
    vector<double>& mu_Q,vector<double>& sigma_Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,int r);
vector<vector<double>> Preprocessing_MC3E_DIVIDE(vector<vector<double>>& Q,
    vector<double>& mu_Q,vector<double>& sigma_Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,int r);
double LB_MC3(vector<vector<double>>& Q, vector<vector<double>>& A,
    vector<double>& mu_Q,vector<double>& sigma_Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<double>& Dbl, vector<double>& Dbr,
    vector<double>& Dtl, vector<double>& Dtr,
    vector<vector<Bucket>>& T_mu,vector<vector<Bucket>>& T_sigma,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,
    int K);
double LB_MC3E(vector<vector<double>>& A,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<double>& mu_A,vector<double>& sigma_A,
    vector<double>& Dbl, vector<double>& Dbr,
    vector<double>& Dtl, vector<double>& Dtr,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num);
double LB_MC_NAIVE1(vector<vector<double>>& Q, vector<vector<double>>& A, vector<double>& mu_Q, int K, int m, int r);
double LB_MC_NAIVE2(vector<vector<double>>& Q, vector<vector<double>>& A,vector<double>& mu_Q, vector<double>& sigma_Q, int K, int m, int r);
vector<vector<vector<double>>> readCSV(const string& filePath);
vector<vector<double>> transpose(const std::vector<std::vector<double>>& matrix);
void splitA_sort(const std::vector<std::vector<double>>& A,
                 const std::vector<int>& split_indices,
                 std::vector<std::vector<double>>& A1,
                 std::vector<std::vector<double>>& A2);
void splitQ_sort(const std::vector<std::vector<double>>& Q,
            std::vector<std::vector<double>>& Q1,
            std::vector<std::vector<double>>& Q2,
            vector<int>& split_indices,
            int r);
void split_half(const std::vector<std::vector<double>>& Q,
            std::vector<std::vector<double>>& Q1,
            std::vector<std::vector<double>>& Q2);
void split_third(const std::vector<std::vector<double>>& Q,
                 std::vector<std::vector<double>>& Q1,
                 std::vector<std::vector<double>>& Q2,
                 std::vector<std::vector<double>>& Q3);
void splitQ_correlation(const std::vector<std::vector<double>>& Q,
                        std::vector<std::vector<double>>& Q1,
                        std::vector<std::vector<double>>& Q2,
                        std::vector<int>& split_indices,
                        int r);
Eigen::MatrixXd pca_reduce_with_projection(const std::vector<std::vector<double>>& data, int target_dim, std::vector<std::vector<double>>& reduced_data,Eigen::VectorXd& saved_mean);
void reduce_A_with_Q_projection(const std::vector<std::vector<double>>& A,
                                const Eigen::MatrixXd& projection_matrix,
                                const Eigen::VectorXd& mean_Q,
                                std::vector<std::vector<double>>& A_reduced);
void splitQ_min_metric_heuristic(const std::vector<std::vector<double>>& Q,
                                 std::vector<std::vector<double>>& Q1,
                                 std::vector<std::vector<double>>& Q2,
                                 std::vector<int>& split_indices);
void splitQ_min_metric(const std::vector<std::vector<double>>& Q,
                       std::vector<std::vector<double>>& Q1,
                       std::vector<std::vector<double>>& Q2,
                       std::vector<int>& split_indices);
void split_two_dim(const std::vector<std::vector<double>>& Q,
            vector<double>& Q1,
            vector<double>& Q2);
void splitQ_min_metric_4_to_2(const std::vector<std::vector<double>>& Q,
                       std::vector<std::vector<double>>& Q1,
                       std::vector<std::vector<double>>& Q2,
                       std::vector<int>& split_indices);
vector<vector<double>> Preprocessing_MC3E_PCA(vector<vector<double>>& Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,int r,
    vector<Eigen::MatrixXd>& projection_matrix_vector,vector<Eigen::VectorXd>& mean_Q_vector);
void lower_upper_lemire_pca(int m, vector<double>& mu_A, vector<double>& sigma_A,
                  double& U_mu_Q, double& L_mu_Q, double& U_sigma_Q, double& L_sigma_Q);
void reduce_single_A_with_Q_projection(const vector<double>& A,
                                       const Eigen::MatrixXd& projection_matrix,
                                       const Eigen::VectorXd& mean_Q,
                                       vector<double>& A_reduced);
double LB_MC3E_2dim_pca(vector<vector<double>>& A,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<double>& Dbl, vector<double>& Dbr,
    vector<double>& Dtl, vector<double>& Dtr,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,
    vector<Eigen::MatrixXd>& projection_matrix_vector,vector<Eigen::VectorXd>& mean_Q_vector);
vector<vector<vector<double>>> Preprocessing_MC3E_PCA_4dim(vector<vector<double>>& Q,
    vector<double>& U_mu_Q1,vector<double>& L_mu_Q1,
    vector<double>& U_sigma_Q1,vector<double>& L_sigma_Q1,
    vector<double>& U_mu_Q2,vector<double>& L_mu_Q2,
    vector<double>& U_sigma_Q2,vector<double>& L_sigma_Q2,
    vector<vector<vector<Cell>>>& T1,
    vector <double>& intervalLength_mu1,vector <double>& intervalLength_sigma1 ,
    vector<vector<vector<Cell>>>& T2,
    vector <double>& intervalLength_mu2,vector <double>& intervalLength_sigma2 ,vector <int>& num,int r,
    vector<Eigen::MatrixXd>& projection_matrix_vector,vector<Eigen::VectorXd>& mean_Q_vector);
double LB_MC3E_4dim_pca(vector<vector<double>>& A,
    vector<double>& U_mu_Q1,vector<double>& L_mu_Q1,
    vector<double>& U_sigma_Q1,vector<double>& L_sigma_Q1,
    vector<double>& U_mu_Q2,vector<double>& L_mu_Q2,
    vector<double>& U_sigma_Q2,vector<double>& L_sigma_Q2,
    vector<double>& Dbl1, vector<double>& Dbr1,
    vector<double>& Dtl1, vector<double>& Dtr1,
    vector<double>& Dbl2, vector<double>& Dbr2,
    vector<double>& Dtl2, vector<double>& Dtr2,
    vector<vector<vector<Cell>>>& T1,
    vector <double>& intervalLength_mu1,vector <double>& intervalLength_sigma1 ,
    vector<vector<vector<Cell>>>& T2,
    vector <double>& intervalLength_mu2,vector <double>& intervalLength_sigma2 ,vector <int>& num,
    vector<Eigen::MatrixXd>& projection_matrix_vector,vector<Eigen::VectorXd>& mean_Q_vector);
void split_first_and_second(const std::vector<std::vector<double>>& Q,
            std::vector<std::vector<double>>& Q1);
double LB_MC3E_sort(
    vector<double>& mu_A,vector<double>& sigma_A,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<double>& Dbl, vector<double>& Dbr,
    vector<double>& Dtl, vector<double>& Dtr,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,
    int m,int K);
void splitA_sort_single(const vector<double>& A,
                 const std::vector<int>& split_indices,
                 vector<double>& A1,
                 vector<double>& A2);
void mu_sigma_single(vector<double>& A, double& mu_A, double& sigma_A);
double LB_MC3E_sort_strength(vector<vector<double>>& A,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<double>& Dbl, vector<double>& Dbr,
    vector<double>& Dtl, vector<double>& Dtr,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,
    vector<vector<int>> split_indices_sort);
vector<vector<double>> Preprocessing_MC3E_SORT_strength(vector<vector<double>>& Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,int r,
    vector<vector<int>> &split_indices_sort);
void splitQ_sort_strength(const std::vector<std::vector<double>>& Q,
            std::vector<std::vector<double>>& Q1,
            std::vector<std::vector<double>>& Q2,
            vector<int>& split_indices);
vector<vector<vector<double>>> Preprocessing_MC3E_sort_4dim_strength(vector<vector<double>>& Q,
    vector<double>& U_mu_Q1,vector<double>& L_mu_Q1,
    vector<double>& U_sigma_Q1,vector<double>& L_sigma_Q1,
    vector<double>& U_mu_Q2,vector<double>& L_mu_Q2,
    vector<double>& U_sigma_Q2,vector<double>& L_sigma_Q2,
    vector<vector<vector<Cell>>>& T1,
    vector <double>& intervalLength_mu1,vector <double>& intervalLength_sigma1 ,
    vector<vector<vector<Cell>>>& T2,
    vector <double>& intervalLength_mu2,vector <double>& intervalLength_sigma2 ,vector <int>& num,int r,
    vector<vector<int>> &split_indices_sort);
double LB_MC3E_4dim_sort_strength(vector<vector<double>>& A,
    vector<double>& U_mu_Q1,vector<double>& L_mu_Q1,
    vector<double>& U_sigma_Q1,vector<double>& L_sigma_Q1,
    vector<double>& U_mu_Q2,vector<double>& L_mu_Q2,
    vector<double>& U_sigma_Q2,vector<double>& L_sigma_Q2,
    vector<double>& Dbl1, vector<double>& Dbr1,
    vector<double>& Dtl1, vector<double>& Dtr1,
    vector<double>& Dbl2, vector<double>& Dbr2,
    vector<double>& Dtl2, vector<double>& Dtr2,
    vector<vector<vector<Cell>>>& T1,
    vector <double>& intervalLength_mu1,vector <double>& intervalLength_sigma1 ,
    vector<vector<vector<Cell>>>& T2,
    vector <double>& intervalLength_mu2,vector <double>& intervalLength_sigma2 ,vector <int>& num,
     vector<vector<int>> split_indices_sort);
vector<vector<vector<double>>> Preprocessing_MC3E_sigma_4dim_strength(vector<vector<double>>& Q,
    vector<double>& U_mu_Q1,vector<double>& L_mu_Q1,
    vector<double>& U_sigma_Q1,vector<double>& L_sigma_Q1,
    vector<double>& U_mu_Q2,vector<double>& L_mu_Q2,
    vector<double>& U_sigma_Q2,vector<double>& L_sigma_Q2,
    vector<vector<vector<Cell>>>& T1,
    vector <double>& intervalLength_mu1,vector <double>& intervalLength_sigma1 ,
    vector<vector<vector<Cell>>>& T2,
    vector <double>& intervalLength_mu2,vector <double>& intervalLength_sigma2 ,vector <int>& num,int r,
    vector<vector<int>> &split_indices_sort);
Eigen::MatrixXd lda_reduce_with_projection(const vector<vector<double>>& data,
                                           const vector<int>& labels,
                                           int target_dim,
                                           vector<vector<double>>& reduced_data);
void reduce_A_with_Q_lda_projection(const vector<vector<double>>& A,
                                    const Eigen::MatrixXd& projection_matrix,
                                    vector<vector<double>>& A_reduced);
double LB_MC3E_4dim_sort_strength_mumu_sigmasigma(vector<vector<double>>& A,
    vector<double>& U_mu_Q1,vector<double>& L_mu_Q1,
    vector<double>& U_sigma_Q1,vector<double>& L_sigma_Q1,
    vector<double>& U_mu_Q2,vector<double>& L_mu_Q2,
    vector<double>& U_sigma_Q2,vector<double>& L_sigma_Q2,
    vector<double>& Dbl1, vector<double>& Dbr1,
    vector<double>& Dtl1, vector<double>& Dtr1,
    vector<double>& Dbl2, vector<double>& Dbr2,
    vector<double>& Dtl2, vector<double>& Dtr2,
    vector<vector<vector<Cell>>>& T1,
    vector <double>& intervalLength_mu1,vector <double>& intervalLength_sigma1 ,
    vector<vector<vector<Cell>>>& T2,
    vector <double>& intervalLength_mu2,vector <double>& intervalLength_sigma2 ,vector <int>& num,
     vector<vector<int>> split_indices_sort);
vector<vector<vector<double>>> Preprocessing_MC3E_sort_4dim_strength_mumu_sigmasigma(vector<vector<double>>& Q,
    vector<double>& U_mu_Q1,vector<double>& L_mu_Q1,
    vector<double>& U_sigma_Q1,vector<double>& L_sigma_Q1,
    vector<double>& U_mu_Q2,vector<double>& L_mu_Q2,
    vector<double>& U_sigma_Q2,vector<double>& L_sigma_Q2,
    vector<vector<vector<Cell>>>& T1,
    vector <double>& intervalLength_mu1,vector <double>& intervalLength_sigma1 ,
    vector<vector<vector<Cell>>>& T2,
    vector <double>& intervalLength_mu2,vector <double>& intervalLength_sigma2 ,vector <int>& num,int r,
    vector<vector<int>> &split_indices_sort);
double LB_TI(const vector<vector<double>>& X,const vector<vector<double>>& Y,int W,int P);
double LB_PC(const vector<vector<double>>& Q,
                  const vector<vector<double>>& R,vector<vector<vector<vector<double>>>>& bboxes);
vector<vector<vector<vector<double>>>> LB_PC_Preprocess(
    const vector<vector<double>>& s1,
    int windowSize,
    int K,
    double Q);
double LB_MV_test(vector<vector<double>>& Q, vector<vector<double>>& A, int K, int m, int r);
double distance(const std::vector<double>& a, const std::vector<double>& b);
double LB_P_MV(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r);
double LB_WEB_MV(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,vector<vector<double>>& L_UQ,vector<vector<double>>& U_LQ,int r);
#endif
