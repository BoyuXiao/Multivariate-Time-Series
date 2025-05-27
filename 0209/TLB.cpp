#include <iostream>
#include <fstream>
#include "funcs.h"
#include <ctime>
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
    std::ofstream outputFile_data;

    double TLB_MV = 0;
    double TLB_MC1 = 0;
    double TLB_MC2 = 0;
    double TLB_MC3 = 0;
    double TLB_MC3E = 0;
    double TLB_MC3_PLUS = 0;
    double TLB_MC_NAIVE1 = 0;
    double TLB_MC_NAIVE2 = 0;

    double Time_MV = 0;
    double Time_MC1 = 0;
    double Time_MC2 = 0;
    double Time_MC3 = 0;
    double Time_MC3E = 0;
    double Time_MC3_PLUS = 0;
    double Time_MC_NAIVE1 = 0;
    double Time_MC_NAIVE2 = 0;

    double pruning_rate_MV = 0;
    double pruning_rate_MC1 = 0;
    double pruning_rate_MC2 = 0;
    double pruning_rate_MC3 = 0;
    double pruning_rate_MC3E = 0;
    double pruning_rate_MC3_PLUS = 0;
    double pruning_rate_MC_NAIVE1 = 0;
    double pruning_rate_MC_NAIVE2 = 0;

    vector<double> lb_mv_vector(n2);
    vector<double> lb_mc1_vector(n2);
    vector<double> lb_mc2_vector(n2);
    vector<double> lb_mc3_vector(n2);
    vector<double> lb_mc3e_vector(n2);
    vector<double> lb_mc3_plus_vector(n2);
    vector<double> lb_mc_naive_1_vector(n2);
    vector<double> lb_mc_naive_2_vector(n2);

    vector<double> mu_Q(m);
    vector<double> sigma_Q(m);
    vector<double> mu_Q1(m);
    vector<double> sigma_Q1(m);
    vector<double> mu_Q2(m);
    vector<double> sigma_Q2(m);
    vector<double> U_mu_Q, L_mu_Q, U_sigma_Q, L_sigma_Q;
    vector<double> U_mu_Q1, L_mu_Q1, U_sigma_Q1, L_sigma_Q1;
    vector<double> U_mu_Q2, L_mu_Q2, U_sigma_Q2, L_sigma_Q2;
    vector<double> U_mu_QE, L_mu_QE, U_sigma_QE, L_sigma_QE;
    vector<vector<double>> result,resulte;

    vector<double> intervalLength_mu(m);
    vector<double> intervalLength_sigma(m);
    vector<double> intervalLength_mu1(m);
    vector<double> intervalLength_sigma1(m);
    vector<double> intervalLength_mu2(m);
    vector<double> intervalLength_sigma2(m);
    vector<double> intervalLength_muE(m);
    vector<double> intervalLength_sigmaE(m);
    vector<int> num(m);

    outputFile_data.open("..//data.txt",std::ios::app);

    for (int w = 0; w < r_num.size(); w++) {
        int r = ceil(r_num[w]*m);
        TLB_MV = 0;
        TLB_MC1 = 0;
        TLB_MC2 = 0;
        TLB_MC3 = 0;
        TLB_MC3E = 0;
        TLB_MC_NAIVE1 = 0;
        TLB_MC_NAIVE2 = 0;

        Time_MV = 0;
        Time_MC1 = 0;
        Time_MC2 = 0;
        Time_MC3 = 0;
        Time_MC3E = 0;
        Time_MC_NAIVE1 = 0;
        Time_MC_NAIVE2 = 0;

        pruning_rate_MV = 0;
        pruning_rate_MC1 = 0;
        pruning_rate_MC2 = 0;
        pruning_rate_MC3 = 0;
        pruning_rate_MC3E = 0;
        pruning_rate_MC_NAIVE1 = 0;
        pruning_rate_MC_NAIVE2 = 0;


        //一起计算各个下界的TLB------------------------------------------------------------
        for (int i = 0; i <n1; ++i) {
            mu_sigma(Test[i],K,m,mu_Q,sigma_Q);
            vector<vector<Bucket>> T_mu(m);
            vector<vector<Bucket>> T_sigma(m);
            vector<vector<vector<Cell>>> T(m);
            result = Preprocessing_MC3(Test[i],mu_Q,sigma_Q,U_mu_Q,L_mu_Q,U_sigma_Q,L_sigma_Q,T_mu,T_sigma,intervalLength_mu,intervalLength_sigma,num,r);
            // resulte = Preprocessing_MC3E(Test[i],mu_Q,sigma_Q,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,T,intervalLength_muE,intervalLength_sigmaE,num,r);
            // ---------------
            vector<vector<double>> Q1, Q2;
            vector<vector<int>> split_indices; // 记录每个时间点的拆分规则
            splitQ(Test[i], Q1, Q2, split_indices, r);
            int K1 = Q1[0].size();
            int K2 = Q2[0].size();
            // 步骤2：预处理 Q1 和 Q2
            vector<vector<Bucket>> T_mu1(m), T_sigma1(m);
            vector<vector<Bucket>> T_mu2(m), T_sigma2(m);
            mu_sigma(Q1, K1, m, mu_Q1, sigma_Q1);
            auto result1 = Preprocessing_MC3(Q1, mu_Q1, sigma_Q1, U_mu_Q1, L_mu_Q1, U_sigma_Q1, L_sigma_Q1,T_mu1,T_sigma1,intervalLength_mu1,intervalLength_sigma1,num,r);
            mu_sigma(Q2, K2, m, mu_Q2, sigma_Q2);
            auto result2 = Preprocessing_MC3(Q2, mu_Q2, sigma_Q2, U_mu_Q2,L_mu_Q2,U_sigma_Q2,L_sigma_Q2,T_mu2,T_sigma2,intervalLength_mu2,intervalLength_sigma2,num,r);
            // ----------
            double epsilon_MV = MAX_VAL;
            double epsilon_MC1 = MAX_VAL;
            double epsilon_MC2 = MAX_VAL;
            double epsilon_MC3 = MAX_VAL;
            double epsilon_MC3E = MAX_VAL;
            double epsilon_MC3_PLUS = MAX_VAL;
            double epsilon_MC_NAIVE1 = MAX_VAL;
            double epsilon_MC_NAIVE2 = MAX_VAL;

            double count_mv = 0;
            double count_mc1 = 0;
            double count_mc2 = 0;
            double count_mc3 = 0;
            double count_mc3e = 0;
            double count_mc3_plus = 0;
            double count_mc_naive_1 = 0;
            double count_mc_naive_2 = 0;



            for (int j = 0; j < n2; ++j) {

                // double lb_mv = LB_MV(Test[i], Train[j], K, m, r);
                // double lb_mc1 = LB_MC1(Test[i], Train[j],mu_Q, K, m, r);
                // double lb_mc2 = LB_MC2(Test[i], Train[j],mu_Q,sigma_Q, K, m, r);
                // double lb_mc_naive1 = LB_MC_NAIVE1(Test[i], Train[j],mu_Q, K, m, r);
                // double lb_mc_naive2 = LB_MC_NAIVE2(Test[i], Train[j],mu_Q,sigma_Q, K, m, r);
                // double lb_mc3 = LB_MC3(Test[i], Train[j], mu_Q,sigma_Q,U_mu_Q,L_mu_Q,U_sigma_Q,L_sigma_Q,result[0],result[1],result[2],result[3],T_mu,T_sigma,intervalLength_mu,intervalLength_sigma,num,K);
                // double lb_mc3e = LB_MC3E(Test[i], Train[j], mu_Q,sigma_Q,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,resulte[0],resulte[1],resulte[2],resulte[3],T,intervalLength_muE,intervalLength_sigmaE,num,K);

                vector<vector<double>> A1, A2;
                splitA(Train[j], split_indices, A1, A2);
                double lb_mc3_plus = LB_MC3(Q1, A1, mu_Q1,sigma_Q1,U_mu_Q1,L_mu_Q1,U_sigma_Q1,L_sigma_Q1,result1[0],result1[1],result1[2],result1[3],T_mu1,T_sigma1,intervalLength_mu1,intervalLength_sigma1,num,K1)+
                    LB_MC3(Q2, A2, mu_Q2,sigma_Q2,U_mu_Q2,L_mu_Q2,U_sigma_Q2,L_sigma_Q2,result2[0],result2[1],result2[2],result2[3],T_mu2,T_sigma2,intervalLength_mu2,intervalLength_sigma2,num,K2);

                double dtwDistance = dtw(Test[i], Train[j], K, m, r);

                // lb_mv_vector[j] = lb_mv;
                // lb_mc1_vector[j] = lb_mc1;
                // lb_mc2_vector[j] = lb_mc2;
                // lb_mc3_vector[j] = lb_mc3;
                // lb_mc3e_vector[j] = lb_mc3e;
                lb_mc3_plus_vector[j] = lb_mc3_plus;
                // lb_mc_naive_1_vector[j] = lb_mc_naive1;
                // lb_mc_naive_2_vector[j] = lb_mc_naive2;

                //计算pruning_rate
                if (dtwDistance < epsilon_MV) {
                    epsilon_MV  = dtwDistance;
                }
                if (dtwDistance < epsilon_MC1) {
                    epsilon_MC1  = dtwDistance;
                }
                if (dtwDistance < epsilon_MC2) {
                    epsilon_MC2  = dtwDistance;
                }
                if (dtwDistance < epsilon_MC3) {
                    epsilon_MC3  = dtwDistance;
                }
                if (dtwDistance < epsilon_MC3E) {
                    epsilon_MC3E  = dtwDistance;
                }
                if (dtwDistance < epsilon_MC3_PLUS) {
                    epsilon_MC3_PLUS  = dtwDistance;
                }
                if (dtwDistance < epsilon_MC_NAIVE1) {
                    epsilon_MC_NAIVE1 = dtwDistance;
                }
                if (dtwDistance < epsilon_MC_NAIVE2) {
                    epsilon_MC_NAIVE2 = dtwDistance;
                }


                // TLB_MV += (lb_mv/dtwDistance);
                // TLB_MC1 += (lb_mc1/dtwDistance);
                // TLB_MC2 += (lb_mc2/dtwDistance);
                // TLB_MC3 += (lb_mc3/dtwDistance);
                // TLB_MC3E += (lb_mc3e/dtwDistance);
                TLB_MC3_PLUS += (lb_mc3_plus/dtwDistance);
                // TLB_MC_NAIVE1 += (lb_mc_naive1/dtwDistance);
                // TLB_MC_NAIVE2 += (lb_mc_naive2/dtwDistance);
            }

            for (int j = 0; j < n2; ++j) {
                // cout<<lb_mv_vector[j]<<" "<<epsilon_MV<<endl;
                if (lb_mv_vector[j] > epsilon_MV)count_mv+=1;
                if (lb_mc1_vector[j] > epsilon_MC1)count_mc1+=1;
                if (lb_mc2_vector[j] > epsilon_MC2)count_mc2+=1;
                if (lb_mc3_vector[j] > epsilon_MC3)count_mc3+=1;
                if (lb_mc3e_vector[j] > epsilon_MC3E)count_mc3e+=1;
                if (lb_mc3_plus_vector[j] > epsilon_MC3)count_mc3_plus+=1;
                if (lb_mc_naive_1_vector[j] > epsilon_MC_NAIVE1)count_mc_naive_1+=1;
                if (lb_mc_naive_2_vector[j] > epsilon_MC_NAIVE2)count_mc_naive_2+=1;
            }
            count_mv /= n2;
            // cout<<"count_mv:"<<count_mv<<endl;
            count_mc1 /= n2;
            count_mc2 /= n2;
            count_mc3 /= n2;
            count_mc3e /= n2;
            count_mc3_plus /= n2;
            count_mc_naive_1 /= n2;
            count_mc_naive_2 /= n2;

            pruning_rate_MV += count_mv;
            // cout<<"pruning_rate_MV:"<<pruning_rate_MV<<endl;
            pruning_rate_MC1 += count_mc1;
            pruning_rate_MC2 += count_mc2;
            pruning_rate_MC3 += count_mc3;
            pruning_rate_MC3E += count_mc3e;
            pruning_rate_MC3_PLUS += count_mc3_plus;
            pruning_rate_MC_NAIVE1 += count_mc_naive_1;
            pruning_rate_MC_NAIVE2 += count_mc_naive_2;

            U_mu_Q.clear(), L_mu_Q.clear(), U_sigma_Q.clear(), L_sigma_Q.clear();
            U_mu_QE.clear(), L_mu_QE.clear(), U_sigma_QE.clear(), L_sigma_QE.clear();

        }

        TLB_MV /= (n1*n2);
        TLB_MC1 /= (n1*n2);
        TLB_MC2 /= (n1*n2);
        TLB_MC3 /= (n1*n2);
        TLB_MC3E /= (n1*n2);
        TLB_MC3_PLUS /= (n1*n2);
        TLB_MC_NAIVE1 /= (n1*n2);
        TLB_MC_NAIVE2 /= (n1*n2);

        pruning_rate_MV /= (n1);
        pruning_rate_MC1 /= (n1);
        pruning_rate_MC2 /= (n1);
        pruning_rate_MC3 /= (n1);
        pruning_rate_MC3E /= (n1);
        pruning_rate_MC3_PLUS /= (n1);
        pruning_rate_MC_NAIVE1 /= (n1);
        pruning_rate_MC_NAIVE2 /= (n1);




        //计算各个下界的运行时间------------------------------------------------------------
        // 计算 LB_MV 时间 数据存储方式特殊 要特殊处理
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
        //     vector<vector<Bucket>> T_mu(m);
        //     vector<vector<Bucket>> T_sigma(m);
        //     mu_sigma(Test[i], K, m, mu_Q, sigma_Q);
        //     result = Preprocessing_MC3(Test[i],mu_Q,sigma_Q,U_mu_Q,L_mu_Q,U_sigma_Q,L_sigma_Q,T_mu,T_sigma,intervalLength_mu,intervalLength_sigma,num,r);
        //     for (int j = 0; j < n2; ++j) {
        //         LB_MC3(Test[i], Train[j], mu_Q,sigma_Q,U_mu_Q,L_mu_Q,U_sigma_Q,L_sigma_Q,result[0],result[1],result[2],result[3],T_mu,T_sigma,intervalLength_mu,intervalLength_sigma,num,K);
        //     }
        //     U_mu_Q.clear(), L_mu_Q.clear(), U_sigma_Q.clear(), L_sigma_Q.clear();
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
        // start = clock();
        // for (int i = 0; i < n1; ++i) {
        //     mu_sigma(Test[i],K,m,mu_Q,sigma_Q);
        //     vector<vector<Bucket>> T_mu(m);
        //     vector<vector<Bucket>> T_sigma(m);
        //     vector<vector<vector<Cell>>> T(m);
        //     resulte = Preprocessing_MC3E(Test[i],mu_Q,sigma_Q,U_mu_Q,L_mu_Q,U_sigma_Q,L_sigma_Q,T,intervalLength_muE,intervalLength_sigmaE,num,r);
        //     for (int j = 0; j < n2; ++j) {
        //         LB_MC3E(Test[i], Train[j], mu_Q,sigma_Q,U_mu_Q,L_mu_Q,U_sigma_Q,L_sigma_Q,resulte[0],resulte[1],resulte[2],resulte[3],T,intervalLength_muE,intervalLength_sigmaE,num,K);
        //     }
        //     U_mu_QE.clear(), L_mu_QE.clear(), U_sigma_QE.clear(), L_sigma_QE.clear();
        // }
        // end = clock();
        // Time_MC3E = static_cast<double>(end - start) / CLOCKS_PER_SEC;

        // start = clock();
        // for (int i = 0; i < n1; ++i) {
        //     vector<vector<Bucket>> T_mu(m);
        //     vector<vector<Bucket>> T_sigma(m);
        //     mu_sigma(Test[i], K, m, mu_Q, sigma_Q);
        //     result = Preprocessing_MC3(Test[i],mu_Q,sigma_Q,U_mu_Q,L_mu_Q,U_sigma_Q,L_sigma_Q,T_mu,T_sigma,intervalLength_mu,intervalLength_sigma,num,r);
        //     for (int j = 0; j < n2; ++j) {
        //         LB_MC3_plus(Test[i], Train[j], mu_Q,sigma_Q,U_mu_Q,L_mu_Q,U_sigma_Q,L_sigma_Q,result[0],result[1],result[2],result[3],T_mu,T_sigma,intervalLength_mu,intervalLength_sigma,num,K);
        //     }
        //     U_mu_Q.clear(), L_mu_Q.clear(), U_sigma_Q.clear(), L_sigma_Q.clear();
        // }
        // end = clock();
        // Time_MC3_PLUS = static_cast<double>(end - start) / CLOCKS_PER_SEC;


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
                   << pruning_rate_MC3_PLUS << "\n";
        std::cout <<name<<" "<<"w="<<r_num[w]<<" "<< "Results have been written to output.txt" << std::endl;
        outputFile.close();
    }

    outputFile_data.close();

    outputFile.open("..//output.txt",std::ios::app);
    outputFile << "\n";
    outputFile.close();
}


