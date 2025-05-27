#include <iostream>
#include <fstream>
#include "funcs.h"
#include <ctime>

//TODO 2dim可以只算一次
void TLB(vector<vector<vector<double>>>& Test, vector<vector<vector<double>>>& Train,string name) {
    int n1 = Test.size();
    int n2 = Train.size();
    int m = Test[0].size();
    int K = Test[0][0].size();
    vector<double> r_num = {0.01,0.02,0.05,0.1,0.2};

    clock_t start;
    clock_t end;
    clock_t lb_mv_time;

    std::ofstream outputFile;
    std::ofstream outputFile111;
    // outputFile111.open("..//data.txt",std::ios::app);

    double TLB_MV = 0;
    double TLB_MC1 = 0;
    double TLB_MC2 = 0;
    double TLB_MC3 = 0;
    double TLB_MC3E = 0;
    double TLB_MC3_PLUS = 0;
    double TLB_MC3_PLUS_SORT = 0;
    double TLB_MC3_PLUS_COV = 0;
    double TLB_MC3_PLUS_SIGMA = 0;
    double TLB_PCA = 0;
    double TLB_PCA4 = 0;
    double TLB_MC_NAIVE1 = 0;
    double TLB_MC_NAIVE2 = 0;

    double Time_DTW = 0;
    double Time_MV = 0;
    double Time_MC1 = 0;
    double Time_MC2 = 0;
    double Time_MC3 = 0;
    double Time_MC3E = 0;
    double Time_MC3_PLUS = 0;
    double Time_MC3_PLUS_SORT = 0;
    double Time_MC3_PLUS_COV = 0;
    double Time_MC3_PLUS_SIGMA = 0;
    double Time_PCA = 0;
    double Time_PCA4 = 0;
    double Time_MC_NAIVE1 = 0;
    double Time_MC_NAIVE2 = 0;

    double pruning_rate_MV = 0;
    double pruning_rate_MC1 = 0;
    double pruning_rate_MC2 = 0;
    double pruning_rate_MC3 = 0;
    double pruning_rate_MC3E = 0;
    double pruning_rate_MC3_PLUS = 0;
    double pruning_rate_MC3_PLUS_SORT = 0;
    double pruning_rate_MC3_PLUS_COV = 0;
    double pruning_rate_MC3_PLUS_SIGMA = 0;
    double pruning_rate_PCA = 0;
    double pruning_rate_PCA4 = 0;
    double pruning_rate_MC_NAIVE1 = 0;
    double pruning_rate_MC_NAIVE2 = 0;

    double epsilon = MAX_VAL;

    double count_mv = 0;
    double count_mc1 = 0;
    double count_mc2 = 0;
    double count_mc3 = 0;
    double count_mc3e = 0;
    double count_mc3_plus = 0;
    double count_mc3_plus_sort = 0;
    double count_mc3_plus_cov = 0;
    double count_mc3_plus_sigma = 0;
    double count_pca = 0;
    double count_pca4 = 0;
    double count_mc_naive_1 = 0;
    double count_mc_naive_2 = 0;

    vector<double> lb_mv_vector(n2);
    vector<double> lb_mc1_vector(n2);
    vector<double> lb_mc2_vector(n2);
    vector<double> lb_mc3_vector(n2);
    vector<double> lb_mc3e_vector(n2);
    vector<double> lb_mc3e_divide_vector(n2);
    vector<double> lb_mc3_plus_vector(n2);
    vector<double> lb_mc3_plus_sort_vector(n2);
    vector<double> lb_mc3_plus_cov_vector(n2);
    vector<double> lb_mc3_plus_sigma_vector(n2);
    vector<double> lb_pca_vector(n2);
    vector<double> lb_pca4_vector(n2);
    vector<double> lb_mc_naive_1_vector(n2);
    vector<double> lb_mc_naive_2_vector(n2);

    vector<double> mu_Q(m);
    vector<double> sigma_Q(m);
    vector<double> mu_Q1(m);
    vector<double> sigma_Q1(m);
    vector<double> mu_Q2(m);
    vector<double> sigma_Q2(m);
    vector<double> mu_Q1_sort(m);
    vector<double> sigma_Q1_sort(m);
    vector<double> mu_Q2_sort(m);
    vector<double> sigma_Q2_sort(m);
    vector<double> mu_Q1_cov(m);
    vector<double> sigma_Q1_cov(m);
    vector<double> mu_Q2_cov(m);
    vector<double> sigma_Q2_cov(m);
    vector<double> mu_Q1_sigma(m);
    vector<double> sigma_Q1_sigma(m);
    vector<double> mu_Q2_sigma(m);
    vector<double> sigma_Q2_sigma(m);

    vector<double> mu_A(m);
    vector<double> sigma_A(m);
    vector<double> mu_A1(m);
    vector<double> sigma_A1(m);
    vector<double> mu_A2(m);
    vector<double> sigma_A2(m);
    vector<double> mu_A1_sort(m);
    vector<double> sigma_A1_sort(m);
    vector<double> mu_A2_sort(m);
    vector<double> sigma_A2_sort(m);
    vector<double> mu_A1_cov(m);
    vector<double> sigma_A1_cov(m);
    vector<double> mu_A2_cov(m);
    vector<double> sigma_A2_cov(m);
    vector<double> mu_A1_sigma(m);
    vector<double> sigma_A1_sigma(m);
    vector<double> mu_A2_sigma(m);
    vector<double> sigma_A2_sigma(m);

    vector<double> mu_Q_pca(m);
    vector<double> sigma_Q_pca(m);
    vector<double> mu_Q_pca4(m);
    vector<double> sigma_Q_pca4(m);
    vector<double> U_mu_Q, L_mu_Q, U_sigma_Q, L_sigma_Q;
    vector<double> U_mu_Q1, L_mu_Q1, U_sigma_Q1, L_sigma_Q1;
    vector<double> U_mu_Q2, L_mu_Q2, U_sigma_Q2, L_sigma_Q2;
    vector<double> U_mu_Q1_sort, L_mu_Q1_sort, U_sigma_Q1_sort, L_sigma_Q1_sort;
    vector<double> U_mu_Q2_sort, L_mu_Q2_sort, U_sigma_Q2_sort, L_sigma_Q2_sort;
    vector<double> U_mu_Q1_cov, L_mu_Q1_cov, U_sigma_Q1_cov, L_sigma_Q1_cov;
    vector<double> U_mu_Q2_cov, L_mu_Q2_cov, U_sigma_Q2_cov, L_sigma_Q2_cov;
    vector<double> U_mu_Q1_sigma, L_mu_Q1_sigma, U_sigma_Q1_sigma, L_sigma_Q1_sigma;
    vector<double> U_mu_Q2_sigma, L_mu_Q2_sigma, U_sigma_Q2_sigma, L_sigma_Q2_sigma;
    vector<double> U_mu_Q_pca, L_mu_Q_pca, U_sigma_Q_pca, L_sigma_Q_pca;
    vector<double> U_mu_Q_pca4_1, L_mu_Q_pca4_1, U_sigma_Q_pca4_1, L_sigma_Q_pca4_1;
    vector<double> U_mu_Q_pca4_2, L_mu_Q_pca4_2, U_sigma_Q_pca4_2, L_sigma_Q_pca4_2;
    vector<double> U_mu_Q_pca4_3, L_mu_Q_pca4_3, U_sigma_Q_pca4_3, L_sigma_Q_pca4_3;
    vector<double> U_mu_QE, L_mu_QE, U_sigma_QE, L_sigma_QE;
    vector<vector<double>> result,resulte,result1,result2,result1_sort,result2_sort,result1_cov,result2_cov,result1_sigma,result2_sigma,result_pca,result_pca4_1,result_pca4_2,result_pca4_3;
    vector<vector<vector<double>>>result_pca4;

    vector<double> intervalLength_mu(m);
    vector<double> intervalLength_sigma(m);
    vector<double> intervalLength_mu1(m);
    vector<double> intervalLength_sigma1(m);
    vector<double> intervalLength_mu2(m);
    vector<double> intervalLength_sigma2(m);
    vector<double> intervalLength_mu1_sort(m);
    vector<double> intervalLength_sigma1_sort(m);
    vector<double> intervalLength_mu2_sort(m);
    vector<double> intervalLength_sigma2_sort(m);
    vector<double> intervalLength_mu1_cov(m);
    vector<double> intervalLength_sigma1_cov(m);
    vector<double> intervalLength_mu2_cov(m);
    vector<double> intervalLength_sigma2_cov(m);
    vector<double> intervalLength_mu1_sigma(m);
    vector<double> intervalLength_sigma1_sigma(m);
    vector<double> intervalLength_mu2_sigma(m);
    vector<double> intervalLength_sigma2_sigma(m);
    vector<double> intervalLength_mu_pca(m);
    vector<double> intervalLength_sigma_pca(m);
    vector<double> intervalLength_mu_pca4_1(m);
    vector<double> intervalLength_sigma_pca4_1(m);
    vector<double> intervalLength_mu_pca4_2(m);
    vector<double> intervalLength_sigma_pca4_2(m);
    vector<double> intervalLength_mu_pca4_3(m);
    vector<double> intervalLength_sigma_pca4_3(m);
    vector<double> intervalLength_muE(m);
    vector<double> intervalLength_sigmaE(m);
    vector<int> num(m);

    vector<vector<double>> Q1(m), Q2(m);
    vector<vector<double>> A1(m), A2(m);
    vector<vector<double>> Q1_sort(m), Q2_sort(m);
    vector<vector<double>> A1_sort(m), A2_sort(m);
    vector<vector<double>> Q1_cov(m), Q2_cov(m);
    vector<vector<double>> A1_cov(m), A2_cov(m);
    vector<vector<double>> Q1_sigma(m), Q2_sigma(m);
    vector<vector<double>> A1_sigma(m), A2_sigma(m);
    vector<vector<double>> Q_pca, A_pca;
    vector<vector<double>> Q_pca4, A_pca4;
    vector<vector<double>> Q_pca4_1(m), Q_pca4_2(m),Q_pca4_3(m);
    vector<vector<double>> A_pca4_1(m), A_pca4_2(m),A_pca4_3(m);
    vector<int> split_indices;
    vector<vector<int>> split_indices_sort(m);
    vector<int> split_indices_cov;
    vector<int> split_indices_sigma;
    vector<int> split_indices_pca;
    Eigen::MatrixXd projection_matrix;
    Eigen::VectorXd mean_Q;
    Eigen::MatrixXd projection_matrix4;
    Eigen::VectorXd mean_Q4;
    vector<Eigen::MatrixXd> projection_matrix_vector(m);
    vector<Eigen::VectorXd> mean_Q_vector(m);
    vector<Eigen::MatrixXd> projection_matrix_vector4(m);
    vector<Eigen::VectorXd> mean_Q_vector4(m);
    vector<vector<vector<double>>> Q_pca_part(m);

    int K1,K2;
    int K1_sort,K2_sort;
    int K1_correlation,K2_correlation;
    int K1_sigma,K2_sigma;
    int K_pca;
    int K_pca4;

    vector<double> Q_2dim_0_mc1(m);
    vector<double> Q_2dim_1_mc1(m);
    vector<double> Q_2dim_0_mc2(m);
    vector<double> Q_2dim_1_mc2(m);
    vector<double> Q_2dim_0_mc3(m);
    vector<double> Q_2dim_1_mc3(m);
    vector<double> Q_2dim_0_naive1(m);
    vector<double> Q_2dim_1_naive1(m);
    vector<double> Q_2dim_0_naive2(m);
    vector<double> Q_2dim_1_naive2(m);
    vector<double> Q_2dim_0_mc3e(m);
    vector<double> Q_2dim_1_mc3e(m);
    vector<double> Q_2dim_0_half1(m);
    vector<double> Q_2dim_1_half1(m);
    vector<double> Q_2dim_0_half2(m);
    vector<double> Q_2dim_1_half2(m);
    vector<double> Q_2dim_0_sort1(m);
    vector<double> Q_2dim_1_sort1(m);
    vector<double> Q_2dim_0_sort2(m);
    vector<double> Q_2dim_1_sort2(m);
    vector<double> Q_2dim_0_cov1(m);
    vector<double> Q_2dim_1_cov1(m);
    vector<double> Q_2dim_0_cov2(m);
    vector<double> Q_2dim_1_cov2(m);
    vector<double> Q_2dim_0_sigma1(m);
    vector<double> Q_2dim_1_sigma1(m);
    vector<double> Q_2dim_0_sigma2(m);
    vector<double> Q_2dim_1_sigma2(m);
    vector<vector<double>> Q_2dim_0_pca(m);
    vector<vector<double>> Q_2dim_1_pca(m);
    vector<double> Q_2dim_0_pca4_1(m);
    vector<double> Q_2dim_1_pca4_1(m);
    vector<double> Q_2dim_0_pca4_2(m);
    vector<double> Q_2dim_1_pca4_2(m);
    vector<double> Q_2dim_0_pca4_3(m);
    vector<double> Q_2dim_1_pca4_3(m);



    for (int w = 0; w < r_num.size(); w++) {
        int r = ceil(r_num[w]*m);
        TLB_MV = 0;
        TLB_MC1 = 0;
        TLB_MC2 = 0;
        TLB_MC3 = 0;
        TLB_MC3E = 0;
        TLB_MC3_PLUS = 0;
        TLB_MC3_PLUS_SORT = 0;
        TLB_MC3_PLUS_COV = 0;
        TLB_MC3_PLUS_SIGMA = 0;
        TLB_PCA = 0;
        TLB_PCA4 = 0;
        TLB_MC_NAIVE1 = 0;
        TLB_MC_NAIVE2 = 0;

        Time_DTW = 0;
        Time_MV = 0;
        Time_MC1 = 0;
        Time_MC2 = 0;
        Time_MC3 = 0;
        Time_MC3E = 0;
        Time_MC3_PLUS = 0;
        Time_MC3_PLUS_SORT = 0;
        Time_MC3_PLUS_COV = 0;
        Time_MC3_PLUS_SIGMA = 0;
        Time_PCA = 0;
        Time_PCA4 = 0;
        Time_MC_NAIVE1 = 0;
        Time_MC_NAIVE2 = 0;

        pruning_rate_MV = 0;
        pruning_rate_MC1 = 0;
        pruning_rate_MC2 = 0;
        pruning_rate_MC3 = 0;
        pruning_rate_MC3E = 0;
        pruning_rate_MC3_PLUS = 0;
        pruning_rate_MC3_PLUS_SORT = 0;
        pruning_rate_MC3_PLUS_COV = 0;
        pruning_rate_MC3_PLUS_SIGMA = 0;
        pruning_rate_PCA = 0;
        pruning_rate_PCA4 = 0;
        pruning_rate_MC_NAIVE1 = 0;
        pruning_rate_MC_NAIVE2 = 0;

        //
        // //一起计算各个下界的TLB------------------------------------------------------------
        // for (int i = 0; i <n1; ++i) {
        //     U_mu_Q.clear(), L_mu_Q.clear(), U_sigma_Q.clear(), L_sigma_Q.clear();
        //     U_mu_QE.clear(), L_mu_QE.clear(), U_sigma_QE.clear(), L_sigma_QE.clear();
        //     U_mu_Q1.clear(), L_mu_Q1.clear(), U_sigma_Q1.clear(), L_sigma_Q1.clear();
        //     U_mu_Q2.clear(), L_mu_Q2.clear(), U_sigma_Q2.clear(), L_sigma_Q2.clear();
        //     U_mu_Q1_sort.clear(), L_mu_Q1_sort.clear(), U_sigma_Q1_sort.clear(), L_sigma_Q1_sort.clear();
        //     U_mu_Q2_sort.clear(), L_mu_Q2_sort.clear(), U_sigma_Q2_sort.clear(), L_sigma_Q2_sort.clear();
        //     U_mu_Q1_cov.clear(), L_mu_Q1_cov.clear(), U_sigma_Q1_cov.clear(), L_sigma_Q1_cov.clear();
        //     U_mu_Q2_cov.clear(), L_mu_Q2_cov.clear(), U_sigma_Q2_cov.clear(), L_sigma_Q2_cov.clear();
        //     U_mu_Q1_sigma.clear(), L_mu_Q1_sigma.clear(), U_sigma_Q1_sigma.clear(), L_sigma_Q1_sigma.clear();
        //     U_mu_Q2_sigma.clear(), L_mu_Q2_sigma.clear(), U_sigma_Q2_sigma.clear(), L_sigma_Q2_sigma.clear();
        //     U_mu_Q_pca.clear(), L_mu_Q_pca.clear(), U_sigma_Q_pca.clear(), L_sigma_Q_pca.clear();
        //     U_mu_Q_pca4_1.clear(), L_mu_Q_pca4_1.clear(), U_sigma_Q_pca4_1.clear(), L_sigma_Q_pca4_1.clear();
        //     U_mu_Q_pca4_2.clear(), L_mu_Q_pca4_2.clear(), U_sigma_Q_pca4_2.clear(), L_sigma_Q_pca4_2.clear();
        //     U_mu_Q_pca4_3.clear(), L_mu_Q_pca4_3.clear(), U_sigma_Q_pca4_3.clear(), L_sigma_Q_pca4_3.clear();
        //
        //
        //     mu_sigma(Test[i],mu_Q,sigma_Q);
        //     vector<vector<Bucket>> T_mu(m);
        //     vector<vector<Bucket>> T_sigma(m);
        //     vector<vector<vector<Cell>>> T(m);
        //     vector<vector<vector<Cell>>> T_pca(m);
        //     vector<vector<vector<Cell>>> T_pca4_1(m);
        //     vector<vector<vector<Cell>>> T_pca4_2(m);
        //     vector<vector<vector<Cell>>> T_pca4_3(m);
        //     //half
        //     vector<vector<vector<Cell>>> T1(m);
        //     vector<vector<vector<Cell>>> T2(m);
        //     //SORT
        //     vector<vector<vector<Cell>>> T1_sort(m);
        //     vector<vector<vector<Cell>>> T2_sort(m);
        //     //COV
        //     vector<vector<vector<Cell>>> T1_cov(m);
        //     vector<vector<vector<Cell>>> T2_cov(m);
        //     vector<vector<vector<Cell>>> T1_sigma(m);
        //     vector<vector<vector<Cell>>> T2_sigma(m);
        //
        //
        //
        //     if (MC1) {
        //         if (K==2) {
        //             split_two_dim(Test[i],Q_2dim_0_mc1,Q_2dim_1_mc1);
        //         }
        //     }
        //     if (MC2) {
        //         if (K==2) {
        //             split_two_dim(Test[i],Q_2dim_0_mc2,Q_2dim_1_mc2);
        //         }
        //     }
        //     if (NAIVE1) {
        //         if (K==2) {
        //             split_two_dim(Test[i],Q_2dim_0_naive1,Q_2dim_1_naive1);
        //         }
        //     }
        //     if (NAIVE2) {
        //         if (K==2) {
        //             split_two_dim(Test[i],Q_2dim_0_naive2,Q_2dim_1_naive2);
        //         }
        //     }
        //     if (MC3) {
        //         if (K==2) {
        //             split_two_dim(Test[i],Q_2dim_0_mc3,Q_2dim_1_mc3);
        //             result = Preprocessing_MC3(Test[i],Q_2dim_0_mc3,Q_2dim_1_mc3,U_mu_Q,L_mu_Q,U_sigma_Q,L_sigma_Q,T_mu,T_sigma,intervalLength_mu,intervalLength_sigma,num,r);
        //         }else {
        //             result = Preprocessing_MC3(Test[i],mu_Q,sigma_Q,U_mu_Q,L_mu_Q,U_sigma_Q,L_sigma_Q,T_mu,T_sigma,intervalLength_mu,intervalLength_sigma,num,r);
        //         }
        //     }
        //     if (MC3E) {
        //         if (DIVIDE) {
        //             if (K==2) {
        //                 split_two_dim(Test[i],Q_2dim_0_mc3e,Q_2dim_1_mc3e);
        //                 resulte = Preprocessing_MC3E_DIVIDE(Test[i],Q_2dim_0_mc3e,Q_2dim_1_mc3e,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,T,intervalLength_muE,intervalLength_sigmaE,num,r);
        //             }else {
        //                 resulte = Preprocessing_MC3E_DIVIDE(Test[i],mu_Q,sigma_Q,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,T,intervalLength_muE,intervalLength_sigmaE,num,r);
        //             }
        //         }else {
        //             if (K==2) {
        //                 split_two_dim(Test[i],Q_2dim_0_mc3e,Q_2dim_1_mc3e);
        //                 resulte = Preprocessing_MC3E(Test[i],Q_2dim_0_mc3e,Q_2dim_1_mc3e,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,T,intervalLength_muE,intervalLength_sigmaE,num,r);
        //             }else {
        //                 resulte = Preprocessing_MC3E(Test[i],mu_Q,sigma_Q,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,T,intervalLength_muE,intervalLength_sigmaE,num,r);
        //             }
        //         }
        //     }
        //     if (PLUS_HALF&&K>=4) {
        //         split_half(Test[i], Q1, Q2);
        //         K1 = Q1[0].size();
        //         K2 = Q2[0].size();
        //         mu_sigma(Q1, mu_Q1, sigma_Q1);
        //         mu_sigma(Q2, mu_Q2, sigma_Q2);
        //
        //         if (K1==2) {
        //             split_two_dim(Q1,Q_2dim_0_half1,Q_2dim_1_half1);
        //             result1 = Preprocessing_MC3E(Q1, Q_2dim_0_half1, Q_2dim_1_half1, U_mu_Q1, L_mu_Q1, U_sigma_Q1, L_sigma_Q1,T1,intervalLength_mu1,intervalLength_sigma1,num,r);
        //         }else {
        //             result1 = Preprocessing_MC3E(Q1, mu_Q1, sigma_Q1, U_mu_Q1, L_mu_Q1, U_sigma_Q1, L_sigma_Q1,T1,intervalLength_mu1,intervalLength_sigma1,num,r);
        //         }
        //
        //         if (K2==2) {
        //             split_two_dim(Q2,Q_2dim_0_half2,Q_2dim_1_half2);
        //             result2 = Preprocessing_MC3E(Q2, Q_2dim_0_half2, Q_2dim_1_half2, U_mu_Q2,L_mu_Q2,U_sigma_Q2,L_sigma_Q2,T2,intervalLength_mu2,intervalLength_sigma2,num,r);
        //         }else {
        //             result2 = Preprocessing_MC3E(Q2, mu_Q2, sigma_Q2, U_mu_Q2,L_mu_Q2,U_sigma_Q2,L_sigma_Q2,T2,intervalLength_mu2,intervalLength_sigma2,num,r);
        //         }
        //     }
        //     if (PLUS_SORT&&K>=4) {
        //         // splitQ_sort(Test[i],Q1_sort,Q2_sort,split_indices,r);
        //         // K1_sort = Q1_sort[0].size();
        //         // K2_sort = Q2_sort[0].size();
        //         // mu_sigma(Q1_sort, K1_sort, m, mu_Q1_sort, sigma_Q1_sort);
        //         // mu_sigma(Q2_sort, K2_sort, m, mu_Q2_sort, sigma_Q2_sort);
        //
        //         // if (K1_sort==2) {
        //         //     split_two_dim(Q1_sort,Q_2dim_0_sort1,Q_2dim_1_sort1);
        //         //     result1_sort = Preprocessing_MC3E(Q1_sort, Q_2dim_0_sort1, Q_2dim_1_sort1, U_mu_Q1_sort, L_mu_Q1_sort, U_sigma_Q1_sort, L_sigma_Q1_sort,T1_sort,intervalLength_mu1_sort,intervalLength_sigma1_sort,num,r);
        //         // }else {
        //         //     //sort_mu_mu
        //         //     result1_sort = Preprocessing_MC3E(Q1_sort, mu_Q1_sort, sigma_Q1_sort, U_mu_Q1_sort, L_mu_Q1_sort, U_sigma_Q1_sort, L_sigma_Q1_sort,T1_sort,intervalLength_mu1_sort,intervalLength_sigma1_sort,num,r);
        //         //     // result1_sort = Preprocessing_MC3E(Q1_sort, mu_Q1_sort, sigma_Q1_sort, U_mu_Q1_sort, L_mu_Q1_sort, U_sigma_Q1_sort, L_sigma_Q1_sort,T1_sort,intervalLength_mu1_sort,intervalLength_sigma1_sort,num,r);
        //         // }
        //         // if (K2_sort==2) {
        //         //     split_two_dim(Q2_sort,Q_2dim_0_sort2,Q_2dim_1_sort2);
        //         //     result2_sort = Preprocessing_MC3E(Q2_sort, Q_2dim_0_sort2, Q_2dim_1_sort2, U_mu_Q2_sort, L_mu_Q2_sort, U_sigma_Q2_sort, L_sigma_Q2_sort,T2_sort,intervalLength_mu2_sort,intervalLength_sigma2_sort,num,r);
        //         // }else {
        //         //     //sort_mu_mu
        //         //     result2_sort = Preprocessing_MC3E(Q2_sort, mu_Q2_sort, sigma_Q2_sort, U_mu_Q2_sort, L_mu_Q2_sort, U_sigma_Q2_sort, L_sigma_Q2_sort,T2_sort,intervalLength_mu2_sort,intervalLength_sigma2_sort,num,r);
        //         //     // result2_sort = Preprocessing_MC3E(Q2_sort, mu_Q2_sort, sigma_Q2_sort, U_mu_Q2_sort, L_mu_Q2_sort, U_sigma_Q2_sort, L_sigma_Q2_sort,T2_sort,intervalLength_mu2_sort,intervalLength_sigma2_sort,num,r);
        //         // }
        //
        //
        //         // //sort_mu_mu
        //         // result1_sort = Preprocessing_MC3E(Q1_sort, mu_Q1_sort, mu_Q2_sort, U_mu_Q1_sort, L_mu_Q1_sort, U_sigma_Q1_sort, L_sigma_Q1_sort,T1_sort,intervalLength_mu1_sort,intervalLength_sigma1_sort,num,r);
        //         //sort_mu_mu_strength
        //         result1_sort = Preprocessing_MC3E_SORT_strength(Test[i], U_mu_Q1_sort, L_mu_Q1_sort, U_sigma_Q1_sort, L_sigma_Q1_sort,T1_sort,intervalLength_mu1_sort,intervalLength_sigma1_sort,num,r,split_indices_sort);
        //
        //     }
        //     if (PLUS_COV&&K>=4) {
        //         splitQ_correlation(Test[i],Q1_cov,Q2_cov,split_indices_cov,r);
        //         K1_correlation = Q1_cov[0].size();
        //         K2_correlation = Q2_cov[0].size();
        //         mu_sigma(Q1_cov,  mu_Q1_cov, sigma_Q1_cov);
        //         mu_sigma(Q2_cov,  mu_Q2_cov, sigma_Q2_cov);
        //         if (K1_correlation==2) {
        //             split_two_dim(Q1_cov,Q_2dim_0_cov1,Q_2dim_1_cov1);
        //             result1_cov = Preprocessing_MC3E(Q1_cov, Q_2dim_0_cov1, Q_2dim_1_cov1, U_mu_Q1_cov, L_mu_Q1_cov, U_sigma_Q1_cov, L_sigma_Q1_cov,T1_cov,intervalLength_mu1_cov,intervalLength_sigma1_cov,num,r);
        //         }else {
        //             result1_cov = Preprocessing_MC3E(Q1_cov, mu_Q1_cov, sigma_Q1_cov, U_mu_Q1_cov, L_mu_Q1_cov, U_sigma_Q1_cov, L_sigma_Q1_cov,T1_cov,intervalLength_mu1_cov,intervalLength_sigma1_cov,num,r);
        //         }
        //         if (K2_correlation==2) {
        //             split_two_dim(Q2_cov,Q_2dim_0_cov2,Q_2dim_1_cov2);
        //             result2_cov = Preprocessing_MC3E(Q2_cov, Q_2dim_0_cov2, Q_2dim_1_cov2, U_mu_Q2_cov, L_mu_Q2_cov, U_sigma_Q2_cov, L_sigma_Q2_cov,T2_cov,intervalLength_mu2_cov,intervalLength_sigma2_cov,num,r);
        //         }else {
        //             result2_cov = Preprocessing_MC3E(Q2_cov, mu_Q2_cov, sigma_Q2_cov, U_mu_Q2_cov, L_mu_Q2_cov, U_sigma_Q2_cov, L_sigma_Q2_cov,T2_cov,intervalLength_mu2_cov,intervalLength_sigma2_cov,num,r);
        //         }
        //     }
        //     if (PLUS_SIGMA&&K>=4) {
        //         splitQ_min_metric_heuristic(Test[i],Q1_sigma,Q2_sigma,split_indices_sigma);
        //         // splitQ_min_metric(Test[i],Q1_sigma,Q2_sigma,split_indices_sigma);
        //         K1_sigma = Q1_sigma[0].size();
        //         K2_sigma = Q2_sigma[0].size();
        //         mu_sigma(Q1_sigma, mu_Q1_sigma, sigma_Q1_sigma);
        //         mu_sigma(Q2_sigma, mu_Q2_sigma, sigma_Q2_sigma);
        //         if (K1_sigma==2) {
        //             split_two_dim(Q1_sigma,Q_2dim_0_sigma1,Q_2dim_1_sigma1);
        //             result1_sigma = Preprocessing_MC3E(Q1_sigma, Q_2dim_0_sigma1, Q_2dim_1_sigma1, U_mu_Q1_sigma, L_mu_Q1_sigma, U_sigma_Q1_sigma, L_sigma_Q1_sigma,T1_sigma,intervalLength_mu1_sigma,intervalLength_sigma1_sigma,num,r);
        //         }else {
        //             result1_sigma = Preprocessing_MC3E(Q1_sigma, mu_Q1_sigma, sigma_Q1_sigma, U_mu_Q1_sigma, L_mu_Q1_sigma, U_sigma_Q1_sigma, L_sigma_Q1_sigma,T1_sigma,intervalLength_mu1_sigma,intervalLength_sigma1_sigma,num,r);
        //         }
        //         if (K2_sigma==2) {
        //             split_two_dim(Q2_sigma,Q_2dim_0_sigma2,Q_2dim_1_sigma2);
        //             result2_sigma = Preprocessing_MC3E(Q2_sigma, Q_2dim_0_sigma2, Q_2dim_1_sigma2, U_mu_Q2_sigma, L_mu_Q2_sigma, U_sigma_Q2_sigma, L_sigma_Q2_sigma,T2_sigma,intervalLength_mu2_sigma,intervalLength_sigma2_sigma,num,r);
        //         }else {
        //             result2_sigma = Preprocessing_MC3E(Q2_sigma, mu_Q2_sigma, sigma_Q2_sigma, U_mu_Q2_sigma, L_mu_Q2_sigma, U_sigma_Q2_sigma, L_sigma_Q2_sigma,T2_sigma,intervalLength_mu2_sigma,intervalLength_sigma2_sigma,num,r);
        //         }
        //     }
        //     if (PCA&&K>=2) {
        //         result_pca = Preprocessing_MC3E_PCA(Test[i],U_mu_Q_pca,L_mu_Q_pca,U_sigma_Q_pca,L_sigma_Q_pca,T_pca,intervalLength_mu_pca,intervalLength_sigma_pca,num,r,projection_matrix_vector,mean_Q_vector);
        //     }
        //     if (PCA4&&K>=4) {
        //         result_pca4 = Preprocessing_MC3E_PCA_4dim(Test[i],U_mu_Q_pca4_1,L_mu_Q_pca4_1,U_sigma_Q_pca4_1,L_sigma_Q_pca4_1,U_mu_Q_pca4_2,L_mu_Q_pca4_2,U_sigma_Q_pca4_2,L_sigma_Q_pca4_2,T_pca4_1,intervalLength_mu_pca4_1,intervalLength_sigma_pca4_1,T_pca4_2,intervalLength_mu_pca4_2,intervalLength_sigma_pca4_2,num,r,projection_matrix_vector4,mean_Q_vector4);
        //     }
        //
        //
        //     epsilon = MAX_VAL;
        //
        //     count_mv = 0;
        //     count_mc1 = 0;
        //     count_mc2 = 0;
        //     count_mc3 = 0;
        //     count_mc3e = 0;
        //     count_mc3_plus = 0;
        //     count_mc3_plus_sort = 0;
        //     count_mc3_plus_cov = 0;
        //     count_mc3_plus_sigma = 0;
        //     count_pca = 0;
        //     count_pca4 = 0;
        //     count_mc_naive_1 = 0;
        //     count_mc_naive_2 = 0;
        //     double lb_mv = 0;
        //     double lb_mc1 = 0;
        //     double lb_mc2 = 0;
        //     double lb_mc_naive1 = 0;
        //     double lb_mc_naive2 = 0;
        //     double lb_mc3 = 0;
        //     double lb_mc3e = 0;
        //     double lb_mc3_plus = 0;
        //     double lb_mc3_plus_sort = 0;
        //     double lb_mc3_plus_cov = 0;
        //     double lb_mc3_plus_sigma = 0;
        //     double lb_pca = 0;
        //     double lb_pca4 = 0;
        //
        //
        //
        //     for (int j = 0; j < n2; ++j) {
        //
        //         if (MV) {
        //             lb_mv = LB_MV(Test[i], Train[j], K, m, r);
        //             lb_mv_vector[j] = lb_mv;
        //         }
        //         if (MC1) {
        //             if (K==2) {
        //                 lb_mc1 = LB_MC1(Test[i], Train[j],Q_2dim_0_mc1, K, m, r);
        //             }else {
        //                 lb_mc1 = LB_MC1(Test[i], Train[j],mu_Q, K, m, r);
        //             }
        //             lb_mc1_vector[j] = lb_mc1;
        //         }
        //         if (MC2) {
        //             if (K==2) {
        //                 lb_mc2 = LB_MC2(Test[i], Train[j],Q_2dim_0_mc2,Q_2dim_1_mc2, K, m, r);
        //             }else {
        //                 lb_mc2 = LB_MC2(Test[i], Train[j],mu_Q,sigma_Q, K, m, r);
        //             }
        //             lb_mc2_vector[j] = lb_mc2;
        //         }
        //         if (NAIVE1) {
        //             if (K==2) {
        //                 lb_mc_naive1 = LB_MC_NAIVE1(Test[i], Train[j],Q_2dim_0_naive1, K, m, r);
        //             }else {
        //                 lb_mc_naive1 = LB_MC_NAIVE1(Test[i], Train[j],mu_Q, K, m, r);
        //             }
        //             lb_mc_naive_1_vector[j] = lb_mc_naive1;
        //         }
        //         if (NAIVE2) {
        //             if (K==2) {
        //                 lb_mc_naive2 = LB_MC_NAIVE2(Test[i], Train[j],Q_2dim_0_naive2,Q_2dim_1_naive2, K, m, r);
        //             }else {
        //                 lb_mc_naive2 = LB_MC_NAIVE2(Test[i], Train[j],mu_Q,sigma_Q, K, m, r);
        //             }
        //             lb_mc_naive_2_vector[j] = lb_mc_naive2;
        //         }
        //         if (MC3) {
        //             if (K==2) {
        //                 lb_mc3 = LB_MC3(Test[i], Train[j], Q_2dim_0_mc3,Q_2dim_1_mc3,U_mu_Q,L_mu_Q,U_sigma_Q,L_sigma_Q,result[0],result[1],result[2],result[3],T_mu,T_sigma,intervalLength_mu,intervalLength_sigma,num,K);
        //             }else {
        //                 lb_mc3 = LB_MC3(Test[i], Train[j], mu_Q,sigma_Q,U_mu_Q,L_mu_Q,U_sigma_Q,L_sigma_Q,result[0],result[1],result[2],result[3],T_mu,T_sigma,intervalLength_mu,intervalLength_sigma,num,K);
        //             }
        //             lb_mc3_vector[j] = lb_mc3;
        //         }
        //         if (MC3E) {
        //             if (K==2) {
        //                 lb_mc3e = LB_MC3E(Test[i], Train[j], Q_2dim_0_mc3e,Q_2dim_1_mc3e,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,resulte[0],resulte[1],resulte[2],resulte[3],T,intervalLength_muE,intervalLength_sigmaE,num,K);
        //             }else {
        //                 lb_mc3e = LB_MC3E(Test[i], Train[j], mu_Q,sigma_Q,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,resulte[0],resulte[1],resulte[2],resulte[3],T,intervalLength_muE,intervalLength_sigmaE,num,K);
        //             }
        //             lb_mc3e_vector[j] = lb_mc3e;
        //         }
        //         if (PLUS_HALF&&K>=4) {
        //             // if (i==0&&j==11) {
        //             //     cout<<1;
        //             // }
        //             split_half(Train[j], A1, A2);
        //             double lb_mc3_plus_1,lb_mc3_plus_2;
        //             if (K1==2) {
        //                 // lb_mc3_plus_1 = LB_MC_NAIVE2(Q1,A1,Q_2dim_0_half1,Q_2dim_1_half1,K1,m,r);
        //                 lb_mc3_plus_1 = LB_MC3E(Q1, A1, Q_2dim_0_half1,Q_2dim_1_half1,U_mu_Q1,L_mu_Q1,U_sigma_Q1,L_sigma_Q1,result1[0],result1[1],result1[2],result1[3],T1,intervalLength_mu1,intervalLength_sigma1,num,K1);
        //             }else {
        //                 // lb_mc3_plus_1 = LB_MC_NAIVE2(Q1,A1,mu_Q1,sigma_Q1,K1,m,r);
        //                 lb_mc3_plus_1 = LB_MC3E(Q1, A1, mu_Q1,sigma_Q1,U_mu_Q1,L_mu_Q1,U_sigma_Q1,L_sigma_Q1,result1[0],result1[1],result1[2],result1[3],T1,intervalLength_mu1,intervalLength_sigma1,num,K1);
        //             }
        //             if (K2==2) {
        //                 // lb_mc3_plus_2 = LB_MC_NAIVE2(Q2,A2,Q_2dim_0_half2,Q_2dim_1_half2,K2,m,r);
        //                 lb_mc3_plus_2 = LB_MC3E(Q2, A2, Q_2dim_0_half2,Q_2dim_1_half2,U_mu_Q2,L_mu_Q2,U_sigma_Q2,L_sigma_Q2,result2[0],result2[1],result2[2],result2[3],T2,intervalLength_mu2,intervalLength_sigma2,num,K2);
        //             }else {
        //                 // lb_mc3_plus_2 = LB_MC_NAIVE2(Q2,A2,mu_Q2,sigma_Q2,K2,m,r);
        //                 lb_mc3_plus_2 = LB_MC3E(Q2, A2, mu_Q2,sigma_Q2,U_mu_Q2,L_mu_Q2,U_sigma_Q2,L_sigma_Q2,result2[0],result2[1],result2[2],result2[3],T2,intervalLength_mu2,intervalLength_sigma2,num,K2);
        //             }
        //
        //             //debug
        //             // lb_mc3_plus_1 = LB_MC_NAIVE2(Q1,A1,mu_Q1,sigma_Q1,K1,m,r);
        //             // lb_mc3_plus_2 = LB_MC_NAIVE2(Q2,A2,mu_Q2,sigma_Q2,K2,m,r);
        //
        //             // lb_mc3_plus_1 = LB_MC3E(Q1, A1, mu_Q1,sigma_Q1,U_mu_Q1,L_mu_Q1,U_sigma_Q1,L_sigma_Q1,result1[0],result1[1],result1[2],result1[3],T1,intervalLength_mu1,intervalLength_sigma1,num,K1);
        //             // lb_mc3_plus_2 = LB_MC3E(Q2, A2, mu_Q2,sigma_Q2,U_mu_Q2,L_mu_Q2,U_sigma_Q2,L_sigma_Q2,result2[0],result2[1],result2[2],result2[3],T2,intervalLength_mu2,intervalLength_sigma2,num,K2);
        //
        //             lb_mc3_plus = sqrt(lb_mc3_plus_1*lb_mc3_plus_1+lb_mc3_plus_2*lb_mc3_plus_2);
        //             lb_mc3_plus_vector[j] = lb_mc3_plus;
        //         }
        //         if (PLUS_SORT&&K>=4) {
        //             // splitA_sort(Train[j],split_indices, A1_sort, A2_sort);
        //             // mu_sigma(A1_sort,K1_sort,m,mu_A1_sort,sigma_A1_sort);
        //             // mu_sigma(A2_sort,K2_sort,m,mu_A2_sort,sigma_A2_sort);
        //             //
        //             // //sort_mu_mu
        //             // lb_mc3_plus_sort = LB_MC3E_sort(mu_A1_sort,mu_A2_sort,U_mu_Q1_sort,L_mu_Q1_sort,U_sigma_Q1_sort,L_sigma_Q1_sort,result1_sort[0],result1_sort[1],result1_sort[2],result1_sort[3],T1_sort,intervalLength_mu1_sort,intervalLength_sigma1_sort,num,m,K);
        //             // lb_mc3_plus_sort_vector[j] = lb_mc3_plus_sort;
        //
        //             // double lb_mc3_plus_1,lb_mc3_plus_2;
        //             // if (K1_sort==2) {
        //             //     lb_mc3_plus_1 = LB_MC3E(Q1_sort, A1_sort, Q_2dim_0_sort1,Q_2dim_1_sort1,U_mu_Q1_sort,L_mu_Q1_sort,U_sigma_Q1_sort,L_sigma_Q1_sort,result1_sort[0],result1_sort[1],result1_sort[2],result1_sort[3],T1_sort,intervalLength_mu1_sort,intervalLength_sigma1_sort,num,K1_sort);
        //             // }else {
        //             //     lb_mc3_plus_1 = LB_MC3E(Q1_sort, A1_sort, mu_Q1_sort,sigma_Q1_sort,U_mu_Q1_sort,L_mu_Q1_sort,U_sigma_Q1_sort,L_sigma_Q1_sort,result1_sort[0],result1_sort[1],result1_sort[2],result1_sort[3],T1_sort,intervalLength_mu1_sort,intervalLength_sigma1_sort,num,K1_sort);
        //             // }
        //             // if (K2_sort==2) {
        //             //     lb_mc3_plus_2 = LB_MC3E(Q2_sort, A2_sort, Q_2dim_0_sort2,Q_2dim_1_sort2,U_mu_Q2_sort,L_mu_Q2_sort,U_sigma_Q2_sort,L_sigma_Q2_sort,result2_sort[0],result2_sort[1],result2_sort[2],result2_sort[3],T2_sort,intervalLength_mu2_sort,intervalLength_sigma2_sort,num,K2_sort);
        //             // }else {
        //             //     lb_mc3_plus_2 = LB_MC3E(Q2_sort, A2_sort, mu_Q2_sort,sigma_Q2_sort,U_mu_Q2_sort,L_mu_Q2_sort,U_sigma_Q2_sort,L_sigma_Q2_sort,result2_sort[0],result2_sort[1],result2_sort[2],result2_sort[3],T2_sort,intervalLength_mu2_sort,intervalLength_sigma2_sort,num,K2_sort);
        //             // }
        //             // lb_mc3_plus_sort = sqrt(lb_mc3_plus_1*lb_mc3_plus_1+lb_mc3_plus_2*lb_mc3_plus_2);
        //             // lb_mc3_plus_sort_vector[j] = lb_mc3_plus_sort;
        //
        //             //strength_mu_mu
        //             lb_mc3_plus_sort = LB_MC3E_sort_strength(Train[j],U_mu_Q1_sort,L_mu_Q1_sort,U_sigma_Q1_sort,L_sigma_Q1_sort,result1_sort[0],result1_sort[1],result1_sort[2],result1_sort[3],T1_sort,intervalLength_mu1_sort,intervalLength_sigma1_sort,num,split_indices_sort);
        //             lb_mc3_plus_sort_vector[j] = lb_mc3_plus_sort;
        //         }
        //
        //         if (PLUS_COV&&K>=4) {
        //             splitA_sort(Train[j],split_indices_cov, A1_cov, A2_cov);
        //             double lb_mc3_plus_1,lb_mc3_plus_2;
        //             if (K1_correlation==2) {
        //                 lb_mc3_plus_1 = LB_MC3E(Q1_cov, A1_cov, Q_2dim_0_cov1,Q_2dim_1_cov1,U_mu_Q1_cov,L_mu_Q1_cov,U_sigma_Q1_cov,L_sigma_Q1_cov,result1_cov[0],result1_cov[1],result1_cov[2],result1_cov[3],T1_cov,intervalLength_mu1_cov,intervalLength_sigma1_cov,num,K1_correlation);
        //             }else {
        //                 lb_mc3_plus_1 = LB_MC3E(Q1_cov, A1_cov, mu_Q1_cov,sigma_Q1_cov,U_mu_Q1_cov,L_mu_Q1_cov,U_sigma_Q1_cov,L_sigma_Q1_cov,result1_cov[0],result1_cov[1],result1_cov[2],result1_cov[3],T1_cov,intervalLength_mu1_cov,intervalLength_sigma1_cov,num,K1_correlation);
        //             }if (K2_correlation==2) {
        //                 lb_mc3_plus_2 = LB_MC3E(Q2_cov, A2_cov, Q_2dim_0_cov2,Q_2dim_1_cov2,U_mu_Q2_cov,L_mu_Q2_cov,U_sigma_Q2_cov,L_sigma_Q2_cov,result2_cov[0],result2_cov[1],result2_cov[2],result2_cov[3],T2_cov,intervalLength_mu2_cov,intervalLength_sigma2_cov,num,K2_correlation);
        //             }else {
        //                 lb_mc3_plus_2 = LB_MC3E(Q2_cov, A2_cov, mu_Q2_cov,sigma_Q2_cov,U_mu_Q2_cov,L_mu_Q2_cov,U_sigma_Q2_cov,L_sigma_Q2_cov,result2_cov[0],result2_cov[1],result2_cov[2],result2_cov[3],T2_cov,intervalLength_mu2_cov,intervalLength_sigma2_cov,num,K2_correlation);
        //             }
        //             lb_mc3_plus_cov = sqrt(lb_mc3_plus_1*lb_mc3_plus_1+lb_mc3_plus_2*lb_mc3_plus_2);
        //             lb_mc3_plus_cov_vector[j] = lb_mc3_plus_cov;
        //         }
        //         if (PLUS_SIGMA&&K>=4) {
        //             double lb_mc3_plus_1,lb_mc3_plus_2;
        //             splitA_sort(Train[j],split_indices_sigma, A1_sigma, A2_sigma);
        //             if (K1_sigma==2) {
        //                 lb_mc3_plus_1 = LB_MC3E(Q1_sigma, A1_sigma, Q_2dim_0_sigma1,Q_2dim_1_sigma1,U_mu_Q1_sigma,L_mu_Q1_sigma,U_sigma_Q1_sigma,L_sigma_Q1_sigma,result1_sigma[0],result1_sigma[1],result1_sigma[2],result1_sigma[3],T1_sigma,intervalLength_mu1_sigma,intervalLength_sigma1_sigma,num,K1_sigma);
        //             }else {
        //                 lb_mc3_plus_1 = LB_MC3E(Q1_sigma, A1_sigma, mu_Q1_sigma,sigma_Q1_sigma,U_mu_Q1_sigma,L_mu_Q1_sigma,U_sigma_Q1_sigma,L_sigma_Q1_sigma,result1_sigma[0],result1_sigma[1],result1_sigma[2],result1_sigma[3],T1_sigma,intervalLength_mu1_sigma,intervalLength_sigma1_sigma,num,K1_sigma);
        //             }
        //             if (K2_sigma==2) {
        //                 lb_mc3_plus_2 = LB_MC3E(Q2_sigma, A2_sigma, Q_2dim_0_sigma2,Q_2dim_1_sigma2,U_mu_Q2_sigma,L_mu_Q2_sigma,U_sigma_Q2_sigma,L_sigma_Q2_sigma,result2_sigma[0],result2_sigma[1],result2_sigma[2],result2_sigma[3],T2_sigma,intervalLength_mu2_sigma,intervalLength_sigma2_sigma,num,K2_sigma);
        //             }else {
        //                 lb_mc3_plus_2 = LB_MC3E(Q2_sigma, A2_sigma, mu_Q2_sigma,sigma_Q2_sigma,U_mu_Q2_sigma,L_mu_Q2_sigma,U_sigma_Q2_sigma,L_sigma_Q2_sigma,result2_sigma[0],result2_sigma[1],result2_sigma[2],result2_sigma[3],T2_sigma,intervalLength_mu2_sigma,intervalLength_sigma2_sigma,num,K2_sigma);
        //             }
        //             lb_mc3_plus_sigma = sqrt(lb_mc3_plus_1*lb_mc3_plus_1+lb_mc3_plus_2*lb_mc3_plus_2);
        //             lb_mc3_plus_sigma_vector[j] = lb_mc3_plus_sigma;
        //         }
        //         if (PCA&&K>=2) {
        //             lb_pca = LB_MC3E_2dim_pca(Train[j],U_mu_Q_pca,L_mu_Q_pca,U_sigma_Q_pca,L_sigma_Q_pca,result_pca[0],result_pca[1],result_pca[2],result_pca[3],T_pca,intervalLength_mu_pca,intervalLength_sigma_pca,num,projection_matrix_vector,mean_Q_vector);
        //             lb_pca_vector[j] = lb_pca;
        //         }
        //         if (PCA4&&K>=4) {
        //             lb_pca4 = LB_MC3E_4dim_pca(Train[j],U_mu_Q_pca4_1,L_mu_Q_pca4_1,U_sigma_Q_pca4_1,L_sigma_Q_pca4_1,U_mu_Q_pca4_2,L_mu_Q_pca4_2,U_sigma_Q_pca4_2,L_sigma_Q_pca4_2,result_pca4[0][0],result_pca4[0][1],result_pca4[0][2],result_pca4[0][3],result_pca4[1][0],result_pca4[1][1],result_pca4[1][2],result_pca4[1][3],T_pca4_1,intervalLength_mu_pca4_1,intervalLength_sigma_pca4_1,T_pca4_2,intervalLength_mu_pca4_2,intervalLength_sigma_pca4_2,num,projection_matrix_vector4,mean_Q_vector4);
        //             lb_pca4_vector[j] = lb_pca4;
        //         }
        //
        //         double dtwDistance = dtw(Test[i], Train[j], K, m, r);
        //
        //         //计算pruning_rate
        //         if (dtwDistance < epsilon) {
        //             epsilon  = dtwDistance;
        //         }
        //
        //
        //
        //
        //         // cout<<lb_mc3_plus_cov/dtwDistance<<endl;
        //         if (dtwDistance!=0) {
        //             TLB_MV += (lb_mv/dtwDistance);
        //             // cout<<TLB_MV<<" "<<i<<" "<<j<<endl;
        //             TLB_MC1 += (lb_mc1/dtwDistance);
        //             TLB_MC2 += (lb_mc2/dtwDistance);
        //             TLB_MC3 += (lb_mc3/dtwDistance);
        //             TLB_MC3E += (lb_mc3e/dtwDistance);
        //             TLB_MC3_PLUS += (lb_mc3_plus/dtwDistance);
        //             TLB_MC3_PLUS_SORT += (lb_mc3_plus_sort/dtwDistance);
        //             TLB_MC3_PLUS_COV += (lb_mc3_plus_cov/dtwDistance);
        //             TLB_MC3_PLUS_SIGMA += (lb_mc3_plus_sigma/dtwDistance);
        //             TLB_PCA += (lb_pca/dtwDistance);
        //             TLB_PCA4 += (lb_pca4/dtwDistance);
        //             TLB_MC_NAIVE1 += (lb_mc_naive1/dtwDistance);
        //             TLB_MC_NAIVE2 += (lb_mc_naive2/dtwDistance);
        //         }
        //     }
        //
        //     for (int j = 0; j < n2; ++j) {
        //         if (lb_mv_vector[j] > epsilon)count_mv+=1;
        //         if (lb_mc1_vector[j] > epsilon)count_mc1+=1;
        //         if (lb_mc2_vector[j] > epsilon)count_mc2+=1;
        //         if (lb_mc3_vector[j] > epsilon)count_mc3+=1;
        //         if (lb_mc3e_vector[j] > epsilon)count_mc3e+=1;
        //         if (lb_mc3_plus_vector[j] > epsilon)count_mc3_plus+=1;
        //         if (lb_mc3_plus_sort_vector[j] > epsilon)count_mc3_plus_sort+=1;
        //         if (lb_mc3_plus_cov_vector[j] > epsilon)count_mc3_plus_cov+=1;
        //         if (lb_mc3_plus_sigma_vector[j]>epsilon)count_mc3_plus_sigma+=1;
        //         if (lb_pca_vector[j] > epsilon)count_pca+=1;
        //         if (lb_pca4_vector[j] > epsilon)count_pca4+=1;
        //         if (lb_mc_naive_1_vector[j] > epsilon)count_mc_naive_1+=1;
        //         if (lb_mc_naive_2_vector[j] > epsilon)count_mc_naive_2+=1;
        //     }
        //
        //
        //     // cout<<count_mc3_plus_cov/n2<<endl;
        //     pruning_rate_MV += count_mv;
        //     pruning_rate_MC1 += count_mc1;
        //     pruning_rate_MC2 += count_mc2;
        //     pruning_rate_MC3 += count_mc3;
        //     pruning_rate_MC3E += count_mc3e;
        //     pruning_rate_MC3_PLUS += count_mc3_plus;
        //     pruning_rate_MC3_PLUS_SORT += count_mc3_plus_sort;
        //     pruning_rate_MC3_PLUS_COV += count_mc3_plus_cov;
        //     pruning_rate_MC3_PLUS_SIGMA += count_mc3_plus_sigma;
        //     pruning_rate_PCA += count_pca;
        //     pruning_rate_PCA4 += count_pca4;
        //     pruning_rate_MC_NAIVE1 += count_mc_naive_1;
        //     pruning_rate_MC_NAIVE2 += count_mc_naive_2;
        //
        // }
        //
        //
        // TLB_MV /= (n1*n2);
        // TLB_MC1 /= (n1*n2);
        // TLB_MC2 /= (n1*n2);
        // TLB_MC3 /= (n1*n2);
        // TLB_MC3E /= (n1*n2);
        // TLB_MC3_PLUS /= (n1*n2);
        // TLB_MC3_PLUS_SORT /= (n1*n2);
        // TLB_MC3_PLUS_COV /= (n1*n2);
        // TLB_MC3_PLUS_SIGMA /= (n1*n2);
        // TLB_PCA /= (n1*n2);
        // TLB_PCA4 /= (n1*n2);
        // TLB_MC_NAIVE1 /= (n1*n2);
        // TLB_MC_NAIVE2 /= (n1*n2);
        //
        //
        // pruning_rate_MV /= (n1*n2);
        // pruning_rate_MC1 /= (n1*n2);
        // pruning_rate_MC2 /= (n1*n2);
        // pruning_rate_MC3 /= (n1*n2);
        // pruning_rate_MC3E /= (n1*n2);
        // pruning_rate_MC3_PLUS /= (n1*n2);
        // pruning_rate_MC3_PLUS_SORT /= (n1*n2);
        // pruning_rate_MC3_PLUS_COV /= (n1*n2);
        // pruning_rate_MC3_PLUS_SIGMA /= (n1*n2);
        // pruning_rate_PCA /= (n1*n2);
        // pruning_rate_PCA4 /= (n1*n2);
        // pruning_rate_MC_NAIVE1 /= (n1*n2);
        // pruning_rate_MC_NAIVE2 /= (n1*n2);
        //


        // 计算DTW时间
        start = clock();
        //20倍
        // for (int k=0;k<20;k++) {
            for (int i = 0; i < n1; ++i) {
                only_mu(Test[i], mu_Q);
                for (int j = 0; j < n2; ++j) {
                    dtw(Test[i], Train[j], K, m, r);
                }
            }
        // }
        end = clock();
        Time_DTW = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        //计算各个下界的运行时间------------------------------------------------------------
        // 计算 LB_MV 时间 数据存储方式特殊 要特殊处理
        //
        lb_mv_time = 0;
        //20倍
        // for (int k=0;k<20;k++) {
            for (int i = 0; i < n1; ++i) {
                vector<vector<double>> new_test = transpose(Test[i]);
                int new_m = new_test.size();
                int new_k = new_test[0].size();
                int new_r = ceil(r_num[w]*new_k);
                for (int j = 0; j < n2; ++j) {
                    vector<vector<double>> new_train = transpose(Train[j]);
                    start = clock();
                    LB_MV_T(new_test, new_train, new_k, new_m, new_r);
                    end = clock();
                    lb_mv_time += (end - start);
                }
            }
        // }
        Time_MV = static_cast<double> (lb_mv_time) / CLOCKS_PER_SEC;
        //
        // // 计算 LB_MC1 时间
        // start = clock();
        // for (int i = 0; i < n1; ++i) {
        //     only_mu(Test[i], K, m, mu_Q);
        //     for (int j = 0; j < n2; ++j) {
        //         LB_MC1(Test[i], Train[j], mu_Q, K, m, r);
        //     }
        // }
        // end = clock();
        // Time_MC1 = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        //
        // // 计算 LB_MC2 时间
        // start = clock();
        // for (int i = 0; i < n1; ++i) {
        //     mu_sigma(Test[i], K, m, mu_Q, sigma_Q);
        //     for (int j = 0; j < n2; ++j) {
        //         LB_MC2(Test[i], Train[j], mu_Q, sigma_Q, K, m, r);
        //     }
        // }
        // end = clock();
        // Time_MC2 = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        //
        // // 计算 LB_MC3 时间
        // start = clock();
        // for (int i = 0; i < n1; ++i) {
        //  U_mu_Q.clear(), L_mu_Q.clear(), U_sigma_Q.clear(), L_sigma_Q.clear();
        //     vector<vector<Bucket>> T_mu(m);
        //     vector<vector<Bucket>> T_sigma(m);
        //     mu_sigma(Test[i], K, m, mu_Q, sigma_Q);
        //     result = Preprocessing_MC3(Test[i],mu_Q,sigma_Q,U_mu_Q,L_mu_Q,U_sigma_Q,L_sigma_Q,T_mu,T_sigma,intervalLength_mu,intervalLength_sigma,num,r);
        //     for (int j = 0; j < n2; ++j) {
        //         LB_MC3(Test[i], Train[j], mu_Q,sigma_Q,U_mu_Q,L_mu_Q,U_sigma_Q,L_sigma_Q,result[0],result[1],result[2],result[3],T_mu,T_sigma,intervalLength_mu,intervalLength_sigma,num,K);
        //     }
        //
        // }
        // end = clock();
        // Time_MC3 = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        //
        // // 计算 LB_MC_NAIVE1 时间
        // start = clock();
        // for (int i = 0; i < n1; ++i) {
        //     only_mu(Test[i], K, m, mu_Q);
        //     for (int j = 0; j < n2; ++j) {
        //         LB_MC_NAIVE1(Test[i], Train[j], mu_Q, K, m, r);
        //     }
        // }
        // end = clock();
        // Time_MC_NAIVE1 = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        //
        // // 计算 LB_MC_NAIVE2 时间
        // start = clock();
        // for (int i = 0; i < n1; ++i) {
        //     mu_sigma(Test[i], K, m, mu_Q, sigma_Q);
        //     for (int j = 0; j < n2; ++j) {
        //         LB_MC_NAIVE2(Test[i], Train[j], mu_Q, sigma_Q, K, m, r);
        //     }
        // }
        // end = clock();
        // Time_MC_NAIVE2 = static_cast<double>(end - start) / CLOCKS_PER_SEC;



        //
        // if (DIVIDE) {
        //     start = clock();
        //     for (int i = 0; i < n1; ++i) {
        //         U_mu_QE.clear(), L_mu_QE.clear(), U_sigma_QE.clear(), L_sigma_QE.clear();
        //         mu_sigma(Test[i],K,m,mu_Q,sigma_Q);
        //         vector<vector<Bucket>> T_mu(m);
        //         vector<vector<Bucket>> T_sigma(m);
        //         vector<vector<vector<Cell>>> T(m);
        //         resulte = Preprocessing_MC3E_DIVIDE(Test[i],mu_Q,sigma_Q,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,T,intervalLength_muE,intervalLength_sigmaE,num,r);
        //         for (int j = 0; j < n2; ++j) {
        //             LB_MC3E(Test[i], Train[j], mu_Q,sigma_Q,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,resulte[0],resulte[1],resulte[2],resulte[3],T,intervalLength_muE,intervalLength_sigmaE,num,K);
        //         }
        //
        //     }
        //     end = clock();
        //     Time_MC3E = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        // }else {
        //     start = clock();
        //     for (int i = 0; i < n1; ++i) {
        //         U_mu_QE.clear(), L_mu_QE.clear(), U_sigma_QE.clear(), L_sigma_QE.clear();
        //         mu_sigma(Test[i],K,m,mu_Q,sigma_Q);
        //         vector<vector<Bucket>> T_mu(m);
        //         vector<vector<Bucket>> T_sigma(m);
        //         vector<vector<vector<Cell>>> T(m);
        //         resulte = Preprocessing_MC3E(Test[i],mu_Q,sigma_Q,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,T,intervalLength_muE,intervalLength_sigmaE,num,r);
        //         for (int j = 0; j < n2; ++j) {
        //             LB_MC3E(Test[i], Train[j], mu_Q,sigma_Q,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,resulte[0],resulte[1],resulte[2],resulte[3],T,intervalLength_muE,intervalLength_sigmaE,num,K);
        //         }
        //
        //     }
        //     end = clock();
        //     Time_MC3E = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        // }



        // if (K>=4) {
        //     start = clock();
        //     for (int i = 0; i < n1; ++i) {
        //         U_mu_Q1.clear(), L_mu_Q1.clear(), U_sigma_Q1.clear(), L_sigma_Q1.clear();
        //         U_mu_Q2.clear(), L_mu_Q2.clear(), U_sigma_Q2.clear(), L_sigma_Q2.clear();
        //         vector<vector<vector<Cell>>> T1(m);
        //         vector<vector<vector<Cell>>> T2(m);
        //         split_half(Test[i], Q1, Q2);
        //         int K1 = Q1[0].size();
        //         int K2 = Q2[0].size();
        //         mu_sigma(Q1, K1, m, mu_Q1, sigma_Q1);
        //         mu_sigma(Q2, K2, m, mu_Q2, sigma_Q2);
        //         result1 = Preprocessing_MC3E(Q1, mu_Q1, sigma_Q1, U_mu_Q1, L_mu_Q1, U_sigma_Q1, L_sigma_Q1,T1,intervalLength_mu1,intervalLength_sigma1,num,r);
        //         result2 = Preprocessing_MC3E(Q2, mu_Q2, sigma_Q2, U_mu_Q2,L_mu_Q2,U_sigma_Q2,L_sigma_Q2,T2,intervalLength_mu2,intervalLength_sigma2,num,r);
        //         for (int j = 0; j < n2; ++j) {
        //             split_half(Train[j], A1, A2);
        //             LB_MC3E(Q1, A1, mu_Q1,sigma_Q1,U_mu_Q1,L_mu_Q1,U_sigma_Q1,L_sigma_Q1,result1[0],result1[1],result1[2],result1[3],T1,intervalLength_mu1,intervalLength_sigma1,num,K1);
        //             LB_MC3E(Q2, A2, mu_Q2,sigma_Q2,U_mu_Q2,L_mu_Q2,U_sigma_Q2,L_sigma_Q2,result2[0],result2[1],result2[2],result2[3],T2,intervalLength_mu2,intervalLength_sigma2,num,K2);
        //         }
        //     }
        //     end = clock();
        //     Time_MC3_PLUS = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        // }




        // if (K>=4) {
        //     start = clock();
        //     for (int i = 0; i < n1; ++i) {
        //         U_mu_Q1_sort.clear(), L_mu_Q1_sort.clear(), U_sigma_Q1_sort.clear(), L_sigma_Q1_sort.clear();
        //         U_mu_Q2_sort.clear(), L_mu_Q2_sort.clear(), U_sigma_Q2_sort.clear(), L_sigma_Q2_sort.clear();
        //         vector<vector<vector<Cell>>> T1_sort(m);
        //         vector<vector<vector<Cell>>> T2_sort(m);
        //         splitQ_sort(Test[i],Q1_sort,Q2_sort,split_indices,r);
        //         int K1_sort = Q1_sort[0].size();
        //         int K2_sort = Q2_sort[0].size();
        //         mu_sigma(Q1_sort, K1_sort, m, mu_Q1_sort, sigma_Q1_sort);
        //         mu_sigma(Q2_sort, K2_sort, m, mu_Q2_sort, sigma_Q2_sort);
        //         result1_sort = Preprocessing_MC3E(Q1_sort, mu_Q1_sort, sigma_Q1_sort, U_mu_Q1_sort, L_mu_Q1_sort, U_sigma_Q1_sort, L_sigma_Q1_sort,T1_sort,intervalLength_mu1_sort,intervalLength_sigma1_sort,num,r);
        //         result2_sort = Preprocessing_MC3E(Q2_sort, mu_Q2_sort, sigma_Q2_sort, U_mu_Q2_sort, L_mu_Q2_sort, U_sigma_Q2_sort, L_sigma_Q2_sort,T2_sort,intervalLength_mu2_sort,intervalLength_sigma2_sort,num,r);
        //         for (int j = 0; j < n2; ++j) {
        //             splitA_sort(Train[j],split_indices, A1_sort, A2_sort);
        //             LB_MC3E(Q1_sort, A1_sort, mu_Q1_sort,sigma_Q1_sort,U_mu_Q1_sort,L_mu_Q1_sort,U_sigma_Q1_sort,L_sigma_Q1_sort,result1_sort[0],result1_sort[1],result1_sort[2],result1_sort[3],T1_sort,intervalLength_mu1_sort,intervalLength_sigma1_sort,num,K1_sort);
        //             LB_MC3E(Q2_sort, A2_sort, mu_Q2_sort,sigma_Q2_sort,U_mu_Q2_sort,L_mu_Q2_sort,U_sigma_Q2_sort,L_sigma_Q2_sort,result2_sort[0],result2_sort[1],result2_sort[2],result2_sort[3],T2_sort,intervalLength_mu2_sort,intervalLength_sigma2_sort,num,K2_sort);
        //         }
        //     }
        //     end = clock();
        //     Time_MC3_PLUS_SORT = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        // }
        //


        //       if (K>=4) {
        //     start = clock();
        //     for (int i = 0; i < n1; ++i) {
        //         U_mu_Q1_cov.clear(), L_mu_Q1_cov.clear(), U_sigma_Q1_cov.clear(), L_sigma_Q1_cov.clear();
        //         U_mu_Q2_cov.clear(), L_mu_Q2_cov.clear(), U_sigma_Q2_cov.clear(), L_sigma_Q2_cov.clear();
        //         vector<vector<vector<Cell>>> T1_cov(m);
        //         vector<vector<vector<Cell>>> T2_cov(m);
        //         splitQ_correlation(Test[i],Q1_cov,Q2_cov,split_indices_cov,r);
        //         K1_correlation = Q1_cov[0].size();
        //         K2_correlation = Q2_cov[0].size();
        //         if (K1_correlation&&K2_correlation) {
        //             mu_sigma(Q1_cov, K1_correlation, m, mu_Q1_cov, sigma_Q1_cov);
        //             mu_sigma(Q2_cov, K2_correlation, m, mu_Q2_cov, sigma_Q2_cov);
        //             result1_cov = Preprocessing_MC3E(Q1_cov, mu_Q1_cov, sigma_Q1_cov, U_mu_Q1_cov, L_mu_Q1_cov, U_sigma_Q1_cov, L_sigma_Q1_cov,T1_cov,intervalLength_mu1_cov,intervalLength_sigma1_cov,num,r);
        //             result2_cov = Preprocessing_MC3E(Q2_cov, mu_Q2_cov, sigma_Q2_cov, U_mu_Q2_cov, L_mu_Q2_cov, U_sigma_Q2_cov, L_sigma_Q2_cov,T2_cov,intervalLength_mu2_cov,intervalLength_sigma2_cov,num,r);
        //             for (int j = 0; j < n2; ++j) {
        //                 splitA_sort(Train[j],split_indices_cov, A1_cov, A2_cov);
        //                 LB_MC3E(Q1_cov, A1_cov, mu_Q1_cov,sigma_Q1_cov,U_mu_Q1_cov,L_mu_Q1_cov,U_sigma_Q1_cov,L_sigma_Q1_cov,result1_cov[0],result1_cov[1],result1_cov[2],result1_cov[3],T1_cov,intervalLength_mu1_cov,intervalLength_sigma1_cov,num,K1_correlation);
        //                 LB_MC3E(Q2_cov, A2_cov, mu_Q2_cov,sigma_Q2_cov,U_mu_Q2_cov,L_mu_Q2_cov,U_sigma_Q2_cov,L_sigma_Q2_cov,result2_cov[0],result2_cov[1],result2_cov[2],result2_cov[3],T2_cov,intervalLength_mu2_cov,intervalLength_sigma2_cov,num,K2_correlation);
        //             }
        //         }
        //     }
        //     end = clock();
        //     Time_MC3_PLUS_COV = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        // }
        //
        // if (K>=4) {
        //     start = clock();
        //     for (int i = 0; i < n1; ++i) {
        //         U_mu_Q_pca.clear(), L_mu_Q_pca.clear(), U_sigma_Q_pca.clear(), L_sigma_Q_pca.clear();
        //         vector<vector<vector<Cell>>> T_pca(m);
        //         projection_matrix = pca_reduce_with_projection(Test[i],4,Q_pca,mean_Q);
        //         K_pca = Q_pca[0].size();
        //         mu_sigma(Q_pca, K_pca, m, mu_Q_pca, sigma_Q_pca);
        //         result_pca = Preprocessing_MC3E(Q_pca,mu_Q_pca,sigma_Q_pca,U_mu_Q_pca,L_mu_Q_pca,U_sigma_Q_pca,L_sigma_Q_pca,T_pca,intervalLength_mu_pca,intervalLength_sigma_pca,num,r);
        //         for (int j = 0; j < n2; ++j) {
        //             reduce_A_with_Q_projection(Train[j],projection_matrix,mean_Q,A_pca);
        //             LB_MC3E(Q_pca,A_pca,mu_Q_pca,sigma_Q_pca,U_mu_Q_pca,L_mu_Q_pca,U_sigma_Q_pca,L_sigma_Q_pca,result_pca[0],result_pca[1],result_pca[2],result_pca[3],T_pca,intervalLength_mu_pca,intervalLength_sigma_pca,num,K_pca);
        //         }
        //     }
        //     end = clock();
        //     Time_PCA = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        // }
        //
        // if (K>=6) {
        //     start = clock();
        //     for (int i = 0; i < n1; ++i) {
        //         U_mu_Q_pca4_1.clear(), L_mu_Q_pca4_1.clear(), U_sigma_Q_pca4_1.clear(), L_sigma_Q_pca4_1.clear();
        //         U_mu_Q_pca4_2.clear(), L_mu_Q_pca4_2.clear(), U_sigma_Q_pca4_2.clear(), L_sigma_Q_pca4_2.clear();
        //         U_mu_Q_pca4_3.clear(), L_mu_Q_pca4_3.clear(), U_sigma_Q_pca4_3.clear(), L_sigma_Q_pca4_3.clear();
        //         vector<vector<vector<Cell>>> T_pca4_1(m);
        //         vector<vector<vector<Cell>>> T_pca4_2(m);
        //         vector<vector<vector<Cell>>> T_pca4_3(m);
        //         projection_matrix4 = pca_reduce_with_projection(Test[i],6,Q_pca4,mean_Q4);
        //         // split_half(Q_pca4,Q_pca4_1,Q_pca4_2);
        //         split_third(Q_pca4,Q_pca4_1,Q_pca4_2,Q_pca4_3);
        //         // splitQ_min_metric_4_to_2(Q_pca4,Q_pca4_1,Q_pca4_2,split_indices_pca);
        //         split_two_dim(Q_pca4_1,Q_2dim_0_pca4_1,Q_2dim_1_pca4_1);
        //         split_two_dim(Q_pca4_2,Q_2dim_0_pca4_2,Q_2dim_1_pca4_2);
        //         split_two_dim(Q_pca4_3,Q_2dim_0_pca4_3,Q_2dim_1_pca4_3);
        //         result_pca4_1 = Preprocessing_MC3E(Q_pca4_1,Q_2dim_0_pca4_1,Q_2dim_1_pca4_1,U_mu_Q_pca4_1,L_mu_Q_pca4_1,U_sigma_Q_pca4_1,L_sigma_Q_pca4_1,T_pca4_1,intervalLength_mu_pca4_1,intervalLength_sigma_pca4_1,num,r);
        //         result_pca4_2 = Preprocessing_MC3E(Q_pca4_2,Q_2dim_0_pca4_2,Q_2dim_1_pca4_2,U_mu_Q_pca4_2,L_mu_Q_pca4_2,U_sigma_Q_pca4_2,L_sigma_Q_pca4_2,T_pca4_2,intervalLength_mu_pca4_2,intervalLength_sigma_pca4_2,num,r);
        //         result_pca4_3 = Preprocessing_MC3E(Q_pca4_3,Q_2dim_0_pca4_3,Q_2dim_1_pca4_3,U_mu_Q_pca4_3,L_mu_Q_pca4_3,U_sigma_Q_pca4_3,L_sigma_Q_pca4_3,T_pca4_3,intervalLength_mu_pca4_3,intervalLength_sigma_pca4_3,num,r);
        //         for (int j = 0; j < n2; ++j) {
        //             reduce_A_with_Q_projection(Train[j],projection_matrix4,mean_Q4,A_pca4);
        //             // split_half(A_pca4,A_pca4_1,A_pca4_2);
        //             split_third(A_pca4,A_pca4_1,A_pca4_2,A_pca4_3);
        //             // splitA_sort(A_pca4,split_indices_pca,A_pca4_1,A_pca4_2);
        //             LB_MC3E_2dim(Q_pca4_1,A_pca4_1,Q_2dim_0_pca4_1,Q_2dim_1_pca4_1,U_mu_Q_pca4_1,L_mu_Q_pca4_1,U_sigma_Q_pca4_1,L_sigma_Q_pca4_1,result_pca4_1[0],result_pca4_1[1],result_pca4_1[2],result_pca4_1[3],T_pca4_1,intervalLength_mu_pca4_1,intervalLength_sigma_pca4_1,num,2);
        //             LB_MC3E_2dim(Q_pca4_2,A_pca4_2,Q_2dim_0_pca4_2,Q_2dim_1_pca4_2,U_mu_Q_pca4_2,L_mu_Q_pca4_2,U_sigma_Q_pca4_2,L_sigma_Q_pca4_2,result_pca4_2[0],result_pca4_2[1],result_pca4_2[2],result_pca4_2[3],T_pca4_2,intervalLength_mu_pca4_2,intervalLength_sigma_pca4_2,num,2);
        //             LB_MC3E_2dim(Q_pca4_3,A_pca4_3,Q_2dim_0_pca4_3,Q_2dim_1_pca4_3,U_mu_Q_pca4_3,L_mu_Q_pca4_3,U_sigma_Q_pca4_3,L_sigma_Q_pca4_3,result_pca4_3[0],result_pca4_3[1],result_pca4_3[2],result_pca4_3[3],T_pca4_3,intervalLength_mu_pca4_3,intervalLength_sigma_pca4_3,num,2);
        //         }
        //     }
        //     end = clock();
        //     Time_PCA = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        // }
        //





        //统计结果------------------------------------------------------------
        outputFile.open("..//output.txt",std::ios::app);
        outputFile << name <<"\t"<<Time_DTW<< "\t" << TLB_MV << "\t" << Time_MV << "\t"
                   << pruning_rate_MV << "\t"
                   << TLB_MC1 << "\t" << Time_MC1 << "\t"
                   << pruning_rate_MC1 << "\t"
                   << TLB_MC_NAIVE1 << "\t" << Time_MC_NAIVE1 << "\t"
                   << pruning_rate_MC_NAIVE1 << "\t"
                   << TLB_MC2 << "\t" << Time_MC2 << "\t"
                   << pruning_rate_MC2 << "\t"
                   << TLB_MC_NAIVE2 << "\t" << Time_MC_NAIVE2 << "\t"
                   << pruning_rate_MC_NAIVE2 << "\t"
                   << TLB_MC3 << "\t" << Time_MC3 << "\t"
                   << pruning_rate_MC3 << "\t"
                   << TLB_MC3E << "\t" << Time_MC3E << "\t"
                   << pruning_rate_MC3E << "\t"
                   << TLB_MC3_PLUS << "\t" << Time_MC3_PLUS << "\t"
                   << pruning_rate_MC3_PLUS << "\t"
                   << TLB_MC3_PLUS_SORT << "\t" << Time_MC3_PLUS_SORT << "\t"
                   << pruning_rate_MC3_PLUS_SORT << "\t"
                   << TLB_MC3_PLUS_COV << "\t" << Time_MC3_PLUS_COV << "\t"
                   << pruning_rate_MC3_PLUS_COV << "\t"
                   << TLB_MC3_PLUS_SIGMA << "\t" << Time_MC3_PLUS_SIGMA << "\t"
                   << pruning_rate_MC3_PLUS_SIGMA << "\t"
                   << TLB_PCA << "\t" << Time_PCA << "\t"
                   << pruning_rate_PCA << "\t"
                   << TLB_PCA4 << "\t" << Time_PCA4 << "\t"
                   << pruning_rate_PCA4 << "\n";
        std::cout <<name<<" "<<"w="<<r_num[w]<<" "<< "Results have been written to output.txt" << std::endl;
        outputFile.close();
    }


    outputFile.open("..//output.txt",std::ios::app);
    outputFile << "\n";
    outputFile.close();
}


