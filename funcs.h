
#ifndef MULTI_DTW_FUNCS_H
#define MULTI_DTW_FUNCS_H
#include "def.h"
#include <vector>
#include <string>
#include <cmath>
#include "cell.h"
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
vector<double> getLB_oneQ_qbox(const vector<vector<double>>& X,
                             const vector<vector<vector<double>>>& others,
                             const vector<pair<vector<double>, vector<double>>>& qbounds);
double LB_MV_graph(vector<double>& Q, vector<vector<double>>& A, int K);
void upper_lemire(vector<double>& Q, int m, int r, vector<double>& upperBound);
void lower_lemire(vector<double>& Q, int m, int r, vector<double>& lowerBound);
void only_mu(vector<vector<double>>& A, vector<double>& mu_A);
void mu_sigma(vector<vector<double>>& A, vector<double>& mu_A, vector<double>& sigma_A);
void lower_upper_lemire(vector<double>& Q, int m, int r, vector<double>& lowerBound, vector<double>& upperBound);
void TLB(vector<vector<vector<double>>>& Test, vector<vector<vector<double>>>& Train, const string& name,
         const struct Nn1Config& config);
double euclideanDistance_2(vector<double>& point1, vector<double>& point2);
double dtw(vector<vector<double>>& A, vector<vector<double>>& B, int K, int m, int r);
double dtw(vector<vector<double>>& A, vector<vector<double>>& B, int K, int m, int r,double threshold);
double dtw(vector<vector<double>>& A, vector<vector<double>>& B, int K, int m, int r,double threshold,vector<double>&cb);
double LB_Keogh_2(vector<double>& A, vector<double>& Q, int m, int r);
double LB_Keogh_2_early(vector<double>& A, vector<double>& Q, int m, int r, double best_dist);
double LB_Keogh_2_env_early(vector<double>& A, const vector<double>& q_lower, const vector<double>& q_upper, int m,
                            double best_dist);
double LB_Keogh_old(vector<double>& Q, vector<double>& A, int m, int r);
double LB_MV(vector<vector<double>>& Q, vector<vector<double>>& A, int K, int m, int r);
double LB_MV_early(vector<vector<double>>& Q, vector<vector<double>>& A, int K, int m, int r, double best_dist);
double LB_MV_early_env(vector<vector<double>>& A_t, int K, int m,
                       const vector<vector<double>>& q_lower, const vector<vector<double>>& q_upper, double best_dist);
// double LB_MV_T(vector<vector<double>>& Q, vector<vector<double>>& A, int K, int m, int r);
double LB_MV_T(vector<vector<double>>& lowerBound,vector<vector<double>>& upperBound, vector<vector<double>>& A, int K, int m);
double LB_MC1(vector<vector<double>>& Q, vector<vector<double>>& A,vector<double>& mu_Q, int K, int m, int r);
double LB_MC1_early(vector<vector<double>>& Q, vector<vector<double>>& A,vector<double>& mu_Q, int K, int m, int r, double best_dist);
double LB_MC1_early_cached(const vector<double>& mu_A, const vector<double>& mu_Q,
                           const vector<double>& mu_q_lower, const vector<double>& mu_q_upper,
                           int K, int m, double best_dist);
double LB_MC2(vector<vector<double>>& Q, vector<vector<double>>& A,vector<double>& mu_Q, vector<double>& sigma_Q, int K, int m, int r);
double LB_MC2_early(vector<vector<double>>& Q, vector<vector<double>>& A,vector<double>& mu_Q, vector<double>& sigma_Q, int K, int m, int r, double best_dist);
double LB_MC2_early_cached(const vector<double>& mu_A, const vector<double>& sigma_A,
                           const vector<double>& mu_Q, const vector<double>& sigma_Q,
                           const vector<double>& mu_q_lower, const vector<double>& mu_q_upper,
                           const vector<double>& sigma_q_lower, const vector<double>& sigma_q_upper,
                           int K, int m, double best_dist);
vector<vector<double>> Preprocessing_MC3(vector<vector<double>>& Q,
    vector<double>& mu_Q,vector<double>& sigma_Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,int r);
double LB_MC3(vector<vector<double>>& A,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<double>& mu_A,vector<double>& sigma_A,
    vector<double>& Dbl, vector<double>& Dbr,
    vector<double>& Dtl, vector<double>& Dtr,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num);
double LB_MC3_early(vector<vector<double>>& A,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<double>& mu_A,vector<double>& sigma_A,
    vector<double>& Dbl, vector<double>& Dbr,
    vector<double>& Dtl, vector<double>& Dtr,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,
    double best_dist);
double LB_MC3E(vector<vector<double>>& A,vector<vector<double>>& Q,vector<vector<double>>& A_T,vector<vector<double>>& Q_T,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<double>& mu_A,vector<double>& sigma_A,
    vector<double>& Dbl, vector<double>& Dbr,
    vector<double>& Dtl, vector<double>& Dtr,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,vector<double>&lb_record,int r);
double LB_MC3E_early(vector<vector<double>>& A,vector<vector<double>>& Q,vector<vector<double>>& A_T,vector<vector<double>>& Q_T,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<double>& mu_A,vector<double>& sigma_A,
    vector<double>& Dbl, vector<double>& Dbr,
    vector<double>& Dtl, vector<double>& Dtr,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,vector<double>&lb_record,int r,
    double best_dist);
double LB_MC_NAIVE1(vector<vector<double>>& Q, vector<vector<double>>& A, vector<double>& mu_Q, int K, int m, int r);
double LB_MC_NAIVE2(vector<vector<double>>& Q, vector<vector<double>>& A,vector<double>& mu_Q, vector<double>& sigma_Q, int K, int m, int r);
vector<vector<vector<double>>> readCSV(const string& filePath);
vector<vector<double>> transpose(const std::vector<std::vector<double>>& matrix);
void split_two_dim(const std::vector<std::vector<double>>& Q,
            vector<double>& Q1,
            vector<double>& Q2);
double DTW_a(const vector<vector<double>>& s1,
             const vector<vector<double>>& s2,
             int windowSize,
             double bestdist);
double LB_TI(const vector<vector<double>>& X,const vector<vector<double>>& Y,int W,int P);
double LB_TI_early(const vector<vector<double>>& X,const vector<vector<double>>& Y,int W,int P,double best_dist);
double LB_PC(const vector<vector<double>>& Q,
                  const vector<vector<double>>& R,vector<vector<vector<vector<double>>>>& bboxes);
double LB_PC_early(const vector<vector<double>>& Q,
                  const vector<vector<double>>& R,vector<vector<vector<vector<double>>>>& bboxes,double best_dist);
vector<vector<vector<vector<double>>>> LB_PC_Preprocess(
    const vector<vector<double>>& s1,
    int windowSize,
    int K,
    double Q);
double LB_MV_test(vector<vector<double>>& Q, vector<vector<double>>& A, int K, int m, int r);
double distance(const std::vector<double>& a, const std::vector<double>& b);
double LB_P_MV(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r);
double LB_P_MV_early(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,int r,double best_dist);
double LB_WEB_MV(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,vector<vector<double>>& L_UQ,vector<vector<double>>& U_LQ,int r);
double LB_WEB_MV_early(vector<vector<double>>& A,vector<vector<double>>& Q, vector<vector<double>>& Q_L,vector<vector<double>>& Q_U,vector<vector<double>>& L_UQ,vector<vector<double>>& U_LQ,int r,double best_dist);
#endif
