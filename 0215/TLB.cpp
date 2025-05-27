#include <iostream>
#include <fstream>
#include "funcs.h"
#include <ctime>
#include <Eigen/Dense>
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
    double TLB_MC3E_DIVIDE = 0;
    double TLB_MC3_PLUS = 0;
    double TLB_MC3_PLUS_SORT = 0;
    double TLB_MC3_PLUS_COV = 0;
    double TLB_PCA = 0;
    double TLB_MC_NAIVE1 = 0;
    double TLB_MC_NAIVE2 = 0;

    double Time_MV = 0;
    double Time_MC1 = 0;
    double Time_MC2 = 0;
    double Time_MC3 = 0;
    double Time_MC3E = 0;
    double Time_MC3E_DIVIDE = 0;
    double Time_MC3_PLUS = 0;
    double Time_MC3_PLUS_SORT = 0;
    double Time_MC3_PLUS_COV = 0;
    double Time_PCA = 0;
    double Time_MC_NAIVE1 = 0;
    double Time_MC_NAIVE2 = 0;

    double pruning_rate_MV = 0;
    double pruning_rate_MC1 = 0;
    double pruning_rate_MC2 = 0;
    double pruning_rate_MC3 = 0;
    double pruning_rate_MC3E = 0;
    double pruning_rate_MC3E_DIVIDE = 0;
    double pruning_rate_MC3_PLUS = 0;
    double pruning_rate_MC3_PLUS_SORT = 0;
    double pruning_rate_MC3_PLUS_COV = 0;
    double pruning_rate_PCA = 0;
    double pruning_rate_MC_NAIVE1 = 0;
    double pruning_rate_MC_NAIVE2 = 0;

    double epsilon = MAX_VAL;

    double count_mv = 0;
    double count_mc1 = 0;
    double count_mc2 = 0;
    double count_mc3 = 0;
    double count_mc3e = 0;
    double count_mc3e_divide = 0;
    double count_mc3_plus = 0;
    double count_mc3_plus_sort = 0;
    double count_mc3_plus_cov = 0;
    double count_pca = 0;
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
    vector<double> lb_pca_vector(n2);
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
    vector<double> mu_Q_pca(m);
    vector<double> sigma_Q_pca(m);
    vector<double> U_mu_Q, L_mu_Q, U_sigma_Q, L_sigma_Q;
    vector<double> U_mu_Q1, L_mu_Q1, U_sigma_Q1, L_sigma_Q1;
    vector<double> U_mu_Q2, L_mu_Q2, U_sigma_Q2, L_sigma_Q2;
    vector<double> U_mu_Q1_sort, L_mu_Q1_sort, U_sigma_Q1_sort, L_sigma_Q1_sort;
    vector<double> U_mu_Q2_sort, L_mu_Q2_sort, U_sigma_Q2_sort, L_sigma_Q2_sort;
    vector<double> U_mu_Q1_cov, L_mu_Q1_cov, U_sigma_Q1_cov, L_sigma_Q1_cov;
    vector<double> U_mu_Q2_cov, L_mu_Q2_cov, U_sigma_Q2_cov, L_sigma_Q2_cov;
    vector<double> U_mu_Q_pca, L_mu_Q_pca, U_sigma_Q_pca, L_sigma_Q_pca;
    vector<double> U_mu_QE, L_mu_QE, U_sigma_QE, L_sigma_QE;
    vector<vector<double>> result,resulte,result1,result2,result1_sort,result2_sort,result1_cov,result2_cov,result_pca;

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
    vector<double> intervalLength_mu_pca(m);
    vector<double> intervalLength_sigma_pca(m);
    vector<double> intervalLength_muE(m);
    vector<double> intervalLength_sigmaE(m);
    vector<int> num(m);

    vector<vector<double>> Q1(m), Q2(m);
    vector<vector<double>> A1(m), A2(m);
    vector<vector<double>> Q1_sort(m), Q2_sort(m);
    vector<vector<double>> A1_sort(m), A2_sort(m);
    vector<vector<double>> Q1_cov(m), Q2_cov(m);
    vector<vector<double>> A1_cov(m), A2_cov(m);
    vector<vector<double>> Q_pca(m), A_pca(m);
    vector<int> split_indices;
    vector<int> split_indices_cov;
    Eigen::MatrixXd projection_matrix;

    int K1,K2;
    int K1_sort,K2_sort;
    int K1_correlation,K2_correlation;
    int K_pca;


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
        TLB_PCA = 0;
        TLB_MC_NAIVE1 = 0;
        TLB_MC_NAIVE2 = 0;

        Time_MV = 0;
        Time_MC1 = 0;
        Time_MC2 = 0;
        Time_MC3 = 0;
        Time_MC3E = 0;
        Time_MC3_PLUS = 0;
        Time_MC3_PLUS_SORT = 0;
        Time_MC3_PLUS_COV = 0;
        Time_PCA = 0;
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
        pruning_rate_PCA = 0;
        pruning_rate_MC_NAIVE1 = 0;
        pruning_rate_MC_NAIVE2 = 0;


        //一起计算各个下界的TLB------------------------------------------------------------
        for (int i = 0; i <n1; ++i) {
            U_mu_Q.clear(), L_mu_Q.clear(), U_sigma_Q.clear(), L_sigma_Q.clear();
            U_mu_QE.clear(), L_mu_QE.clear(), U_sigma_QE.clear(), L_sigma_QE.clear();
            U_mu_Q1.clear(), L_mu_Q1.clear(), U_sigma_Q1.clear(), L_sigma_Q1.clear();
            U_mu_Q2.clear(), L_mu_Q2.clear(), U_sigma_Q2.clear(), L_sigma_Q2.clear();
            U_mu_Q1_sort.clear(), L_mu_Q1_sort.clear(), U_sigma_Q1_sort.clear(), L_sigma_Q1_sort.clear();
            U_mu_Q2_sort.clear(), L_mu_Q2_sort.clear(), U_sigma_Q2_sort.clear(), L_sigma_Q2_sort.clear();
            U_mu_Q1_cov.clear(), L_mu_Q1_cov.clear(), U_sigma_Q1_cov.clear(), L_sigma_Q1_cov.clear();
            U_mu_Q2_cov.clear(), L_mu_Q2_cov.clear(), U_sigma_Q2_cov.clear(), L_sigma_Q2_cov.clear();
            U_mu_Q_pca.clear(), L_mu_Q_pca.clear(), U_sigma_Q_pca.clear(), L_sigma_Q_pca.clear();


            mu_sigma(Test[i],K,m,mu_Q,sigma_Q);
            vector<vector<Bucket>> T_mu(m);
            vector<vector<Bucket>> T_sigma(m);
            vector<vector<vector<Cell>>> T(m);
            vector<vector<vector<Cell>>> T_pca(m);
            //half
            vector<vector<vector<Cell>>> T1(m);
            vector<vector<vector<Cell>>> T2(m);
            //SORT
            vector<vector<vector<Cell>>> T1_sort(m);
            vector<vector<vector<Cell>>> T2_sort(m);
            //COV
            vector<vector<vector<Cell>>> T1_cov(m);
            vector<vector<vector<Cell>>> T2_cov(m);




            if (MC3) {
                result = Preprocessing_MC3(Test[i],mu_Q,sigma_Q,U_mu_Q,L_mu_Q,U_sigma_Q,L_sigma_Q,T_mu,T_sigma,intervalLength_mu,intervalLength_sigma,num,r);
            }
            if (MC3E) {
                if (DIVIDE) {
                    resulte = Preprocessing_MC3E_DIVIDE(Test[i],mu_Q,sigma_Q,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,T,intervalLength_muE,intervalLength_sigmaE,num,r);
                }else {
                    resulte = Preprocessing_MC3E(Test[i],mu_Q,sigma_Q,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,T,intervalLength_muE,intervalLength_sigmaE,num,r);
                }
            }
            if (PLUS_HALF&&K>=4) {
                split_half(Test[i], Q1, Q2);
                K1 = Q1[0].size();
                K2 = Q2[0].size();
                mu_sigma(Q1, K1, m, mu_Q1, sigma_Q1);
                mu_sigma(Q2, K2, m, mu_Q2, sigma_Q2);
                result1 = Preprocessing_MC3E(Q1, mu_Q1, sigma_Q1, U_mu_Q1, L_mu_Q1, U_sigma_Q1, L_sigma_Q1,T1,intervalLength_mu1,intervalLength_sigma1,num,r);
                result2 = Preprocessing_MC3E(Q2, mu_Q2, sigma_Q2, U_mu_Q2,L_mu_Q2,U_sigma_Q2,L_sigma_Q2,T2,intervalLength_mu2,intervalLength_sigma2,num,r);
            }
            if (PLUS_SORT&&K>=4) {
                splitQ_sort(Test[i],Q1_sort,Q2_sort,split_indices,r);
                K1_sort = Q1_sort[0].size();
                K2_sort = Q2_sort[0].size();
                mu_sigma(Q1_sort, K1_sort, m, mu_Q1_sort, sigma_Q1_sort);
                mu_sigma(Q2_sort, K2_sort, m, mu_Q2_sort, sigma_Q2_sort);
                result1_sort = Preprocessing_MC3E(Q1_sort, mu_Q1_sort, sigma_Q1_sort, U_mu_Q1_sort, L_mu_Q1_sort, U_sigma_Q1_sort, L_sigma_Q1_sort,T1_sort,intervalLength_mu1_sort,intervalLength_sigma1_sort,num,r);
                result2_sort = Preprocessing_MC3E(Q2_sort, mu_Q2_sort, sigma_Q2_sort, U_mu_Q2_sort, L_mu_Q2_sort, U_sigma_Q2_sort, L_sigma_Q2_sort,T2_sort,intervalLength_mu2_sort,intervalLength_sigma2_sort,num,r);
            }
            if (PLUS_COV&&K>=4) {
                splitQ_correlation(Test[i],Q1_cov,Q2_cov,split_indices_cov,r);
                K1_correlation = Q1_cov[0].size();
                K2_correlation = Q2_cov[0].size();
                mu_sigma(Q1_cov, K1_correlation, m, mu_Q1_cov, sigma_Q1_cov);
                mu_sigma(Q2_cov, K2_correlation, m, mu_Q2_cov, sigma_Q2_cov);
                result1_cov = Preprocessing_MC3E(Q1_cov, mu_Q1_cov, sigma_Q1_cov, U_mu_Q1_cov, L_mu_Q1_cov, U_sigma_Q1_cov, L_sigma_Q1_cov,T1_cov,intervalLength_mu1_cov,intervalLength_sigma1_cov,num,r);
                result2_cov = Preprocessing_MC3E(Q2_cov, mu_Q2_cov, sigma_Q2_cov, U_mu_Q2_cov, L_mu_Q2_cov, U_sigma_Q2_cov, L_sigma_Q2_cov,T2_cov,intervalLength_mu2_cov,intervalLength_sigma2_cov,num,r);
            }
            if (PCA&&K>=4) {
                projection_matrix = pca_reduce_with_projection(Test[i],2,Q_pca);
                K_pca = Q_pca[0].size();
                mu_sigma(Q_pca, K_pca, m, mu_Q_pca, sigma_Q_pca);
                result_pca = Preprocessing_MC3E(Q_pca,mu_Q_pca,sigma_Q_pca,U_mu_Q_pca,L_mu_Q_pca,U_sigma_Q_pca,L_sigma_Q_pca,T_pca,intervalLength_mu_pca,intervalLength_sigma_pca,num,r);
            }


            epsilon = MAX_VAL;

            count_mv = 0;
            count_mc1 = 0;
            count_mc2 = 0;
            count_mc3 = 0;
            count_mc3e = 0;
            count_mc3_plus = 0;
            count_mc3_plus_sort = 0;
            count_mc3_plus_cov = 0;
            count_pca = 0;
            count_mc_naive_1 = 0;
            count_mc_naive_2 = 0;
            double lb_mv = 0;
            double lb_mc1 = 0;
            double lb_mc2 = 0;
            double lb_mc_naive1 = 0;
            double lb_mc_naive2 = 0;
            double lb_mc3 = 0;
            double lb_mc3e = 0;
            double lb_mc3_plus = 0;
            double lb_mc3_plus_sort = 0;
            double lb_mc3_plus_cov = 0;
            double lb_pca = 0;



            for (int j = 0; j < n2; ++j) {
                if (MV) {
                    lb_mv = LB_MV(Test[i], Train[j], K, m, r);
                    lb_mv_vector[j] = lb_mv;
                }
                if (MC1) {
                    lb_mc1 = LB_MC1(Test[i], Train[j],mu_Q, K, m, r);
                    lb_mc1_vector[j] = lb_mc1;
                }
                if (MC2) {
                    lb_mc2 = LB_MC2(Test[i], Train[j],mu_Q,sigma_Q, K, m, r);
                    lb_mc2_vector[j] = lb_mc2;
                }
                if (NAIVE1) {
                    lb_mc_naive1 = LB_MC_NAIVE1(Test[i], Train[j],mu_Q, K, m, r);
                    lb_mc_naive_1_vector[j] = lb_mc_naive1;
                }
                if (NAIVE2) {
                    lb_mc_naive2 = LB_MC_NAIVE2(Test[i], Train[j],mu_Q,sigma_Q, K, m, r);
                    lb_mc_naive_2_vector[j] = lb_mc_naive2;
                }
                if (MC3) {
                    lb_mc3 = LB_MC3(Test[i], Train[j], mu_Q,sigma_Q,U_mu_Q,L_mu_Q,U_sigma_Q,L_sigma_Q,result[0],result[1],result[2],result[3],T_mu,T_sigma,intervalLength_mu,intervalLength_sigma,num,K);
                    lb_mc3_vector[j] = lb_mc3;
                }
                if (MC3E) {
                    lb_mc3e = LB_MC3E(Test[i], Train[j], mu_Q,sigma_Q,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,resulte[0],resulte[1],resulte[2],resulte[3],T,intervalLength_muE,intervalLength_sigmaE,num,K);
                    lb_mc3e_vector[j] = lb_mc3e;
                }
                if (PLUS_HALF&&K>=4) {
                    split_half(Train[j], A1, A2);
                    double lb_mc3_plus_1 = LB_MC3E(Q1, A1, mu_Q1,sigma_Q1,U_mu_Q1,L_mu_Q1,U_sigma_Q1,L_sigma_Q1,result1[0],result1[1],result1[2],result1[3],T1,intervalLength_mu1,intervalLength_sigma1,num,K1);
                    double lb_mc3_plus_2 = LB_MC3E(Q2, A2, mu_Q2,sigma_Q2,U_mu_Q2,L_mu_Q2,U_sigma_Q2,L_sigma_Q2,result2[0],result2[1],result2[2],result2[3],T2,intervalLength_mu2,intervalLength_sigma2,num,K2);
                    lb_mc3_plus = sqrt(lb_mc3_plus_1*lb_mc3_plus_1+lb_mc3_plus_2*lb_mc3_plus_2);
                    lb_mc3_plus_vector[j] = lb_mc3_plus;
                }
                if (PLUS_SORT&&K>=4) {
                    splitA_sort(Train[j],split_indices, A1_sort, A2_sort);
                    double lb_mc3_plus_1 = LB_MC3E(Q1_sort, A1_sort, mu_Q1_sort,sigma_Q1_sort,U_mu_Q1_sort,L_mu_Q1_sort,U_sigma_Q1_sort,L_sigma_Q1_sort,result1_sort[0],result1_sort[1],result1_sort[2],result1_sort[3],T1_sort,intervalLength_mu1_sort,intervalLength_sigma1_sort,num,K1_sort);
                    double lb_mc3_plus_2 = LB_MC3E(Q2_sort, A2_sort, mu_Q2_sort,sigma_Q2_sort,U_mu_Q2_sort,L_mu_Q2_sort,U_sigma_Q2_sort,L_sigma_Q2_sort,result2_sort[0],result2_sort[1],result2_sort[2],result2_sort[3],T2_sort,intervalLength_mu2_sort,intervalLength_sigma2_sort,num,K2_sort);
                    lb_mc3_plus_sort = sqrt(lb_mc3_plus_1*lb_mc3_plus_1+lb_mc3_plus_2*lb_mc3_plus_2);
                    lb_mc3_plus_sort_vector[j] = lb_mc3_plus_sort;
                }
                if (PLUS_COV&&K>=4) {
                    splitA_sort(Train[j],split_indices_cov, A1_cov, A2_cov);
                    double lb_mc3_plus_1 = LB_MC3E(Q1_cov, A1_cov, mu_Q1_cov,sigma_Q1_cov,U_mu_Q1_cov,L_mu_Q1_cov,U_sigma_Q1_cov,L_sigma_Q1_cov,result1_cov[0],result1_cov[1],result1_cov[2],result1_cov[3],T1_cov,intervalLength_mu1_cov,intervalLength_sigma1_cov,num,K1_correlation);
                    double lb_mc3_plus_2 = LB_MC3E(Q2_cov, A2_cov, mu_Q2_cov,sigma_Q2_cov,U_mu_Q2_cov,L_mu_Q2_cov,U_sigma_Q2_cov,L_sigma_Q2_cov,result2_cov[0],result2_cov[1],result2_cov[2],result2_cov[3],T2_cov,intervalLength_mu2_cov,intervalLength_sigma2_cov,num,K2_correlation);
                    lb_mc3_plus_cov = sqrt(lb_mc3_plus_1*lb_mc3_plus_1+lb_mc3_plus_2*lb_mc3_plus_2);
                    lb_mc3_plus_cov_vector[j] = lb_mc3_plus_cov;
                }
                if (PCA&&K>=4) {
                    reduce_A_with_Q_projection(Train[j],projection_matrix,A_pca);
                    lb_pca = LB_MC3E(Q_pca,A_pca,mu_Q_pca,sigma_Q_pca,U_mu_Q_pca,L_mu_Q_pca,U_sigma_Q_pca,L_sigma_Q_pca,result_pca[0],result_pca[1],result_pca[2],result_pca[3],T_pca,intervalLength_mu_pca,intervalLength_sigma_pca,num,K_pca);
                    lb_pca_vector[j] = lb_pca;
                }

                double dtwDistance = dtw(Test[i], Train[j], K, m, r);

                //计算pruning_rate
                if (dtwDistance < epsilon) {
                    epsilon  = dtwDistance;
                }

                TLB_MV += (lb_mv/dtwDistance);
                TLB_MC1 += (lb_mc1/dtwDistance);
                TLB_MC2 += (lb_mc2/dtwDistance);
                TLB_MC3 += (lb_mc3/dtwDistance);
                TLB_MC3E += (lb_mc3e/dtwDistance);
                TLB_MC3_PLUS += (lb_mc3_plus/dtwDistance);
                TLB_MC3_PLUS_SORT += (lb_mc3_plus_sort/dtwDistance);
                TLB_MC3_PLUS_COV += (lb_mc3_plus_cov/dtwDistance);
                TLB_PCA += (lb_pca/dtwDistance);
                TLB_MC_NAIVE1 += (lb_mc_naive1/dtwDistance);
                TLB_MC_NAIVE2 += (lb_mc_naive2/dtwDistance);
            }

            for (int j = 0; j < n2; ++j) {
                if (lb_mv_vector[j] > epsilon)count_mv+=1;
                if (lb_mc1_vector[j] > epsilon)count_mc1+=1;
                if (lb_mc2_vector[j] > epsilon)count_mc2+=1;
                if (lb_mc3_vector[j] > epsilon)count_mc3+=1;
                if (lb_mc3e_vector[j] > epsilon)count_mc3e+=1;
                if (lb_mc3_plus_vector[j] > epsilon)count_mc3_plus+=1;
                if (lb_mc3_plus_sort_vector[j] > epsilon)count_mc3_plus_sort+=1;
                if (lb_mc3_plus_cov_vector[j] > epsilon)count_mc3_plus_cov+=1;
                if (lb_pca_vector[j] > epsilon)count_pca+=1;
                if (lb_mc_naive_1_vector[j] > epsilon)count_mc_naive_1+=1;
                if (lb_mc_naive_2_vector[j] > epsilon)count_mc_naive_2+=1;
            }

            pruning_rate_MV += count_mv;
            pruning_rate_MC1 += count_mc1;
            pruning_rate_MC2 += count_mc2;
            pruning_rate_MC3 += count_mc3;
            pruning_rate_MC3E += count_mc3e;
            pruning_rate_MC3_PLUS += count_mc3_plus;
            pruning_rate_MC3_PLUS_SORT += count_mc3_plus_sort;
            pruning_rate_MC3_PLUS_COV += count_mc3_plus_cov;
            pruning_rate_PCA += count_pca;
            pruning_rate_MC_NAIVE1 += count_mc_naive_1;
            pruning_rate_MC_NAIVE2 += count_mc_naive_2;

        }


        TLB_MV /= (n1*n2);
        TLB_MC1 /= (n1*n2);
        TLB_MC2 /= (n1*n2);
        TLB_MC3 /= (n1*n2);
        TLB_MC3E /= (n1*n2);
        TLB_MC3_PLUS /= (n1*n2);
        TLB_MC3_PLUS_SORT /= (n1*n2);
        TLB_MC3_PLUS_COV /= (n1*n2);
        TLB_PCA /= (n1*n2);
        TLB_MC_NAIVE1 /= (n1*n2);
        TLB_MC_NAIVE2 /= (n1*n2);

        pruning_rate_MV /= (n1*n2);
        pruning_rate_MC1 /= (n1*n2);
        pruning_rate_MC2 /= (n1*n2);
        pruning_rate_MC3 /= (n1*n2);
        pruning_rate_MC3E /= (n1*n2);
        pruning_rate_MC3_PLUS /= (n1*n2);
        pruning_rate_MC3_PLUS_SORT /= (n1*n2);
        pruning_rate_MC3_PLUS_COV /= (n1*n2);
        pruning_rate_PCA /= (n1*n2);
        pruning_rate_MC_NAIVE1 /= (n1*n2);
        pruning_rate_MC_NAIVE2 /= (n1*n2);




        //计算各个下界的运行时间------------------------------------------------------------
        // 计算 LB_MV 时间 数据存储方式特殊 要特殊处理
        //
        // lb_mv_time = 0;
        // for (int i = 0; i < n1; ++i) {
        //     vector<vector<double>> new_test = transpose(Test[i]);
        //     int new_m = new_test.size();
        //     int new_k = new_test[0].size();
        //     int new_r = ceil(r_num[w]*new_m);
        //     for (int j = 0; j < n2; ++j) {
        //         vector<vector<double>> new_train = transpose(Train[j]);
        //         start = clock();
        //         LB_MV_T(new_test, new_train, new_k, new_m, new_r);
        //         end = clock();
        //         lb_mv_time += (end - start);
        //     }
        // }
        // Time_MV = static_cast<double> (lb_mv_time) / CLOCKS_PER_SEC;
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
        //         int K1_correlation = Q1_cov[0].size();
        //         int K2_correlation = Q2_cov[0].size();
        //         mu_sigma(Q1_cov, K1_correlation, m, mu_Q1_cov, sigma_Q1_cov);
        //         mu_sigma(Q2_cov, K2_correlation, m, mu_Q2_cov, sigma_Q2_cov);
        //         result1_cov = Preprocessing_MC3E(Q1_cov, mu_Q1_cov, sigma_Q1_cov, U_mu_Q1_cov, L_mu_Q1_cov, U_sigma_Q1_cov, L_sigma_Q1_cov,T1_cov,intervalLength_mu1_cov,intervalLength_sigma1_cov,num,r);
        //         result2_cov = Preprocessing_MC3E(Q2_cov, mu_Q2_cov, sigma_Q2_cov, U_mu_Q2_cov, L_mu_Q2_cov, U_sigma_Q2_cov, L_sigma_Q2_cov,T2_cov,intervalLength_mu2_cov,intervalLength_sigma2_cov,num,r);
        //         for (int j = 0; j < n2; ++j) {
        //             splitA_sort(Train[j],split_indices_cov, A1_cov, A2_cov);
        //             LB_MC3E(Q1_cov, A1_cov, mu_Q1_cov,sigma_Q1_cov,U_mu_Q1_cov,L_mu_Q1_cov,U_sigma_Q1_cov,L_sigma_Q1_cov,result1_cov[0],result1_cov[1],result1_cov[2],result1_cov[3],T1_cov,intervalLength_mu1_cov,intervalLength_sigma1_cov,num,K1_correlation);
        //             LB_MC3E(Q2_cov, A2_cov, mu_Q2_cov,sigma_Q2_cov,U_mu_Q2_cov,L_mu_Q2_cov,U_sigma_Q2_cov,L_sigma_Q2_cov,result2_cov[0],result2_cov[1],result2_cov[2],result2_cov[3],T2_cov,intervalLength_mu2_cov,intervalLength_sigma2_cov,num,K2_correlation);
        //         }
        //     }
        //     end = clock();
        //     Time_MC3_PLUS_COV = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        // }






        //统计结果------------------------------------------------------------
        outputFile.open("..//output.txt",std::ios::app);
        outputFile << name << "\t" << TLB_MV << "\t" << Time_MV << "\t"
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
                   << TLB_PCA << "\t" << Time_PCA << "\t"
                   << pruning_rate_PCA << "\n";
        std::cout <<name<<" "<<"w="<<r_num[w]<<" "<< "Results have been written to output.txt" << std::endl;
        outputFile.close();
    }


    outputFile.open("..//output.txt",std::ios::app);
    outputFile << "\n";
    outputFile.close();
}


