
#include <iostream>
#include <fstream>
#include "funcs.h"
#include <ctime>

//test GRAPH
void TLB(vector<vector<vector<double>>>& Test, vector<vector<vector<double>>>& Train,string name) {
    int n1 = Test.size();
    int n2 = Train.size();
    int m = Test[0].size();
    int K = Test[0][0].size();


    vector<vector<double>> new_test;
    vector<vector<double>> new_train;
    int new_m;
    int new_k;
    int new_r;
    //
    // vector<vector<double>> new_test_mv;
    // vector<vector<double>> new_train_mv;
    // int new_m_mv;
    // int new_k_mv;
    // int new_r_mv;


    vector<double> r_num = {0.01,0.02,0.05,0.1,0.2};

    clock_t start;
    clock_t end;
    clock_t lb_mv_time;
    clock_t lb_mv_time_pre;

    std::ofstream outputFile;
    std::ofstream outputFile111;

    double TLB_MV = 0;
    double TLB_MV_DOUBLE = 0;
    double TLB_GRAPH = 0;
    double TLB_MC3E = 0;
    double TLB_MV_P = 0;
    double TLB_MV_WEB = 0;


    double Time_DTW = 0;
    double Time_MV = 0;
    double Time_MC3E = 0;
    double Time_MC3E_PRE = 0;
    double Time_D_MV = 0;
    double Time_GRAPH = 0;
    double Time_PC = 0;
    double Time_GRAPH_PRE = 0;
    double Time_PC_PRE = 0;
    double Time_MV_P = 0;
    double Time_MV_WEB = 0;


    vector<double> lb_mv_vector(n2);
    vector<double> lb_mv_double_vector(n2);

    vector<double> mu_Q(m);
    vector<double> sigma_Q(m);

    vector<vector<double>> mu_A(n2,vector<double>(m));
    vector<vector<double>> sigma_A(n2,vector<double>(m));


    vector<double> U_mu_QE(m), L_mu_QE(m), U_sigma_QE(m), L_sigma_QE(m);
    vector<vector<double>> result,resulte;
    // vector<vector<pair<double,int>>> result,resulte;

    vector<double> intervalLength_muE(m);
    vector<double> intervalLength_sigmaE(m);
    vector<int> num(m);


    vector<double> Q_2dim_0_mc3e(m);
    vector<double> Q_2dim_1_mc3e(m);

    vector<vector<double>>Q_L(K,vector<double>(m)),Q_U(K,vector<double>(m));

    vector<vector<double>>L_UQ(K,vector<double>(m)),U_LQ(K,vector<double>(m));



    for (int w = 0; w < r_num.size(); w++) {
        int r = ceil(r_num[w]*m);
        TLB_MV = 0;
        TLB_MV_DOUBLE = 0;
        TLB_GRAPH = 0;
        TLB_MC3E = 0;
        TLB_MV_P = 0;
        TLB_MV_WEB = 0;


        //
        // //一起计算各个下界的TLB------------------------------------------------------------
        for (int i = 0; i <n1; ++i) {

            mu_sigma(Test[i],mu_Q,sigma_Q);
            vector<vector<vector<Cell>>> T(m);

            new_test = transpose(Test[i]);
            new_m = new_test.size();
            new_k = new_test[0].size();
            new_r = ceil(r_num[w]*new_k);

            if (MV||MV_DOUBLE||MV_P||MV_WEB) {
                for (int k =0;k<K;k++) {
                    lower_upper_lemire(new_test[k],m,r,Q_L[k],Q_U[k]);
                }
            }

            if (GRAPH||MC3E) {
                    if (K==2) {
                        split_two_dim(Test[i],Q_2dim_0_mc3e,Q_2dim_1_mc3e);
                        resulte = Preprocessing_MC3E(Test[i],Q_2dim_0_mc3e,Q_2dim_1_mc3e,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,T,intervalLength_muE,intervalLength_sigmaE,num,r);
                    }else {
                        resulte = Preprocessing_MC3E(Test[i],mu_Q,sigma_Q,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,T,intervalLength_muE,intervalLength_sigmaE,num,r);
                    }
                for (int j=0;j<n2;++j) {
                    if (K==2) {
                        split_two_dim(Train[j],mu_A[j],sigma_A[j]);
                    }else {
                        mu_sigma(Train[j],mu_A[j],sigma_A[j]);
                    }
                }
            }

            if (MV_WEB) {
                for (int k =0;k<K;k++) {
                    lower_lemire(Q_U[k],m,r,L_UQ[k]);
                    upper_lemire(Q_L[k],m,r,U_LQ[k]);
                }
            }

            double lb_mv = 0;
            double lb_mv_double = 0;
            double lb_graph = 0;
            double lb_mc3e = 0;
            double lb_mv_p = 0;
            double lb_mv_web = 0;

            for (int j = 0; j < n2; ++j) {
                new_train = transpose(Train[j]);
                if (MV) {
                    // new_train_mv = transpose(Train[j]);
                    lb_mv = LB_MV_T(Q_L,Q_U, new_train, new_k, new_m);
                }
                if (MC3E) {
                    lb_mc3e = LB_MC3E(Train[j],U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,mu_A[j],sigma_A[j],resulte[0],resulte[1],resulte[2],resulte[3],T,intervalLength_muE,intervalLength_sigmaE,num);
                }
                if (GRAPH) {
                    vector<double>lb(m,0);
                    // new_train = transpose(Train[j]);
                    lb_graph = LB_GRAPH(Train[j],Test[i],new_train,new_test,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,mu_A[j],sigma_A[j],resulte[0],resulte[1],resulte[2],resulte[3],T,intervalLength_muE,intervalLength_sigmaE,num,lb,r);
                    // lb_graph = LB_GRAPH(Train[j],Test[i],U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,resulte[0],resulte[1],resulte[2],resulte[3],T,intervalLength_muE,intervalLength_sigmaE,num,lb,r);
                }
                if (MV_DOUBLE) {
                    // lb_mv_double = LB_MV_double_T(new_train,new_test,Q_L,Q_U,new_r);
                    // lb_mv_double = LB_MV_double_KIM_T(new_train,new_test,Q_L,Q_U,new_r);
                    // lb_mv_double = LB_MV_double_T_V2(new_train,new_test,Q_L,Q_U,new_r);
                    // lb_mv_double = LB_MV_double_IDEA_T(new_train,new_test,Q_L,Q_U,new_r);
                    // lb_mv_double = LB_MV_double_IDEA_T_V2(new_train,new_test,Train[j],Test[i],Q_L,Q_U,new_r);

                    //test
                    double t1 = LB_MV_double_IDEA_T(new_train,new_test,Q_L,Q_U,new_r);
                    double t2 = LB_MV_double_IDEA_T_V2(new_train,new_test,Train[j],Test[i],Q_L,Q_U,new_r);
                    if (t1>t2){cout<<t1<<" > "<<t2<<endl;};
                    // lb_mv_double = LB_TI(Train[j],Test[i],r,5);
                }
                if (MV_P) {
                    lb_mv_p = LB_P_MV(new_train,new_test,Q_L,Q_U,new_r);
                }
                if (MV_WEB) {
                    lb_mv_web = LB_WEB_MV(new_train,new_test,Q_L,Q_U,L_UQ,U_LQ,new_r);
                }
                double dtwDistance = dtw(Train[j], Test[i], K, m, r);
                // if (lb_mv_double>dtwDistance) {
                    // cout<<i<<" "<<j<<" "<<lb_mv_double<<" "<<dtwDistance<<endl;
                // }
                // if (lb_graph>dtwDistance) {
                //     cout<<i<<" "<<j<<" "<<lb_graph<<" "<<dtwDistance<<endl;
                // }
                // if (lb_mv_p>dtwDistance) {
                //     cout<<i<<" "<<j<<" "<<lb_mv_p<<" "<<dtwDistance<<endl;
                // }
                // if (lb_mv_web>dtwDistance) {
                //     cout<<i<<" "<<j<<" "<<lb_mv_web<<" "<<dtwDistance<<endl;
                // }
                if (dtwDistance!=0) {
                    TLB_MV += (lb_mv/dtwDistance);
                    TLB_MV_DOUBLE += (lb_mv_double/dtwDistance);
                    TLB_GRAPH += (lb_graph/dtwDistance);
                    TLB_MC3E += (lb_mc3e/dtwDistance);
                    TLB_MV_P += (lb_mv_p/dtwDistance);
                    TLB_MV_WEB += (lb_mv_web/dtwDistance);
                }
            }
        }


        TLB_MV /= (n1*n2);
        TLB_MV_DOUBLE /= (n1*n2);
        TLB_GRAPH /= (n1*n2);
        TLB_MC3E /= (n1*n2);
        TLB_MV_P /= (n1*n2);
        TLB_MV_WEB /= (n1*n2);

      //   // -------------------------------------------------------------------------
      if (Time_DTW_TEST) {
          // 计算DTW时间
          start = clock();
          for (int i = 0; i < n1; ++i) {
              for (int j = 0; j < n2; ++j) {
                  dtw(Test[i], Train[j], K, m, r);
              }
          }
          end = clock();
          Time_DTW = static_cast<double>(end - start) / CLOCKS_PER_SEC;
      }
      //
        if (Time_MV_TEST) {
            // // 计算 LB_MV 时间 数据存储方式特殊 要特殊处理
            lb_mv_time = 0;
            for (int i = 0; i < n1; ++i) {
                new_test = transpose(Test[i]);
                new_m = new_test.size();
                new_k = new_test[0].size();
                new_r = ceil(r_num[w]*new_k);


                start = clock();
                for (int k =0;k<K;k++) {
                    lower_upper_lemire(new_test[k],m,r,Q_L[k],Q_U[k]);
                }
                end = clock();
                lb_mv_time += (end - start);


                for (int j = 0; j < n2; ++j) {
                    new_train = transpose(Train[j]);

                    start = clock();
                    LB_MV_T(Q_L,Q_U, new_train, new_k, new_m);
                    end = clock();
                    lb_mv_time += (end - start);
                }
            }
            Time_MV = static_cast<double> (lb_mv_time) / CLOCKS_PER_SEC;
        }

        if (Time_GRAPH_TEST) {
            lb_mv_time = 0;
            lb_mv_time_pre = 0;

            //先对A压缩 即Q和A的压缩都保证只计算一次
            start = clock();
            for (int j=0;j<n2;++j) {
                if (K==2) {
                    split_two_dim(Train[j],mu_A[j],sigma_A[j]);
                }else {
                    mu_sigma(Train[j],mu_A[j],sigma_A[j]);
                }
            }
            end = clock();
            lb_mv_time_pre += (end - start);


            for (int i = 0; i < n1; ++i) {
                vector<vector<vector<Cell>>> T(m);
                start = clock();
                mu_sigma(Test[i],mu_Q,sigma_Q);
                if (K==2) {
                    split_two_dim(Test[i],Q_2dim_0_mc3e,Q_2dim_1_mc3e);
                    resulte = Preprocessing_MC3E(Test[i],Q_2dim_0_mc3e,Q_2dim_1_mc3e,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,T,intervalLength_muE,intervalLength_sigmaE,num,r);
                }else {
                    resulte = Preprocessing_MC3E(Test[i],mu_Q,sigma_Q,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,T,intervalLength_muE,intervalLength_sigmaE,num,r);
                }
                end = clock();
                lb_mv_time_pre += (end - start);
                for (int j = 0; j < n2; ++j) {
                    vector<double>lb(m,0);
                    start = clock();
                    LB_GRAPH(Train[j],Test[i],new_train,new_test,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,mu_A[j],sigma_A[j],resulte[0],resulte[1],resulte[2],resulte[3],T,intervalLength_muE,intervalLength_sigmaE,num,lb,r);
                    end = clock();
                    lb_mv_time += (end - start);
                }
            }
            // end = clock();
            // Time_GRAPH = static_cast<double>(end - start) / CLOCKS_PER_SEC;
            Time_GRAPH = static_cast<double> (lb_mv_time) / CLOCKS_PER_SEC;
            Time_GRAPH_PRE = static_cast<double> (lb_mv_time_pre) / CLOCKS_PER_SEC;
        }
    if (Time_MC3E_TEST) {
        lb_mv_time = 0;
        lb_mv_time_pre = 0;


        //先对A压缩 即Q和A的压缩都保证只计算一次
        start = clock();
        for (int j=0;j<n2;++j) {
            if (K==2) {
                split_two_dim(Train[j],mu_A[j],sigma_A[j]);
            }else {
                mu_sigma(Train[j],mu_A[j],sigma_A[j]);
            }
        }
        end = clock();
        lb_mv_time_pre += (end - start);


        for (int i = 0; i < n1; ++i) {
            vector<vector<vector<Cell>>> T(m);
            start = clock();
            mu_sigma(Test[i],mu_Q,sigma_Q);
            if (K==2) {
                split_two_dim(Test[i],Q_2dim_0_mc3e,Q_2dim_1_mc3e);
                resulte = Preprocessing_MC3E(Test[i],Q_2dim_0_mc3e,Q_2dim_1_mc3e,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,T,intervalLength_muE,intervalLength_sigmaE,num,r);
            }else {
                resulte = Preprocessing_MC3E(Test[i],mu_Q,sigma_Q,U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,T,intervalLength_muE,intervalLength_sigmaE,num,r);
            }
            end = clock();
            lb_mv_time_pre += (end - start);
            for (int j = 0; j < n2; ++j) {
                start = clock();
                LB_MC3E(Train[j],U_mu_QE,L_mu_QE,U_sigma_QE,L_sigma_QE,mu_A[j],sigma_A[j],resulte[0],resulte[1],resulte[2],resulte[3],T,intervalLength_muE,intervalLength_sigmaE,num);
                end = clock();
                lb_mv_time += (end - start);
            }
        }
        // end = clock();
        // Time_GRAPH = static_cast<double>(end - start) / CLOCKS_PER_SEC;
        Time_MC3E = static_cast<double> (lb_mv_time) / CLOCKS_PER_SEC;
        Time_MC3E_PRE = static_cast<double> (lb_mv_time_pre) / CLOCKS_PER_SEC;
    }

      if (Time_DOUBLE_TEST) {
          lb_mv_time = 0;
          for (int i = 0; i < n1; ++i) {
              new_test = transpose(Test[i]);
              new_m = new_test.size();
              new_k = new_test[0].size();
              new_r = ceil(r_num[w]*new_k);

              start = clock();
              for (int k =0;k<K;k++) {
                  lower_upper_lemire(new_test[k],m,r,Q_L[k],Q_U[k]);
              }
              end = clock();
              lb_mv_time += (end - start);

              for (int j = 0; j < n2; ++j) {
                  new_train = transpose(Train[j]);

                  start = clock();
                  // LB_MV_double_T(new_train,new_test,Q_L,Q_U,new_r);
                  // LB_MV_double_KIM_T(new_train,new_test,Q_L,Q_U,new_r);
                  // LB_MV_double_IDEA_T(new_train,new_test,Q_L,Q_U,new_r);
                  LB_MV_double_IDEA_T_V2(new_train,new_test,Train[j],Test[i],Q_L,Q_U,new_r);
                  end = clock();
                  lb_mv_time += (end - start);

              }

          }
          Time_D_MV = static_cast<double> (lb_mv_time) / CLOCKS_PER_SEC;
      }

        if (Time_MV_P_TEST) {
            lb_mv_time = 0;
            for (int i = 0; i < n1; ++i) {
                new_test = transpose(Test[i]);
                new_m = new_test.size();
                new_k = new_test[0].size();
                new_r = ceil(r_num[w]*new_k);

                start = clock();
                for (int k =0;k<K;k++) {
                    lower_upper_lemire(new_test[k],m,r,Q_L[k],Q_U[k]);
                }
                for (int k =0;k<K;k++) {
                    lower_lemire(Q_U[k],m,r,L_UQ[k]);
                    upper_lemire(Q_L[k],m,r,U_LQ[k]);
                }
                end = clock();
                lb_mv_time += (end - start);

                for (int j = 0; j < n2; ++j) {
                    new_train = transpose(Train[j]);

                    start = clock();
                    LB_P_MV(new_train,new_test,Q_L,Q_U,new_r);
                    end = clock();
                    lb_mv_time += (end - start);

                }

            }
            Time_MV_P = static_cast<double> (lb_mv_time) / CLOCKS_PER_SEC;
        }

        if (Time_MV_WEB_TEST) {
            lb_mv_time = 0;
            for (int i = 0; i < n1; ++i) {
                new_test = transpose(Test[i]);
                new_m = new_test.size();
                new_k = new_test[0].size();
                new_r = ceil(r_num[w]*new_k);

                start = clock();
                for (int k =0;k<K;k++) {
                    lower_upper_lemire(new_test[k],m,r,Q_L[k],Q_U[k]);
                }
                end = clock();
                lb_mv_time += (end - start);

                for (int j = 0; j < n2; ++j) {
                    new_train = transpose(Train[j]);

                    start = clock();
                    LB_WEB_MV(new_train,new_test,Q_L,Q_U,L_UQ,U_LQ,new_r);
                    end = clock();
                    lb_mv_time += (end - start);

                }

            }
            Time_MV_WEB = static_cast<double> (lb_mv_time) / CLOCKS_PER_SEC;
        }

        if (Time_PC_TEST) {
            // 计算 LB_PC 时间
            lb_mv_time = 0;
            lb_mv_time_pre = 0;
            vector<vector<vector<vector<double>>>>bboxes;
            for (int i = 0; i < n1; ++i) {
                start = clock();
                bboxes = LB_PC_Preprocess(Test[i],r,6,2);
                end = clock();
                lb_mv_time_pre += (end - start);
                for (int j = 0; j < n2; ++j) {
                    start = clock();
                    LB_PC(Test[i],Train[j],bboxes);
                    end = clock();
                    lb_mv_time += (end - start);
                }
            }
            Time_PC = static_cast<double> (lb_mv_time) / CLOCKS_PER_SEC;
            Time_PC_PRE = static_cast<double> (lb_mv_time_pre) / CLOCKS_PER_SEC;
        }


        //统计结果------------------------------------------------------------
        outputFile.open("..//output.txt",std::ios::app);
        outputFile << name <<"\t"<<Time_DTW<< "\t" << TLB_MV << "\t" << Time_MV << "\t"
                   << TLB_MV_DOUBLE << "\t" << Time_D_MV << "\t"<< TLB_GRAPH << "\t" << Time_GRAPH<<"\t"
                   <<  Time_GRAPH_PRE<<"\t"<<Time_PC<<"\t" <<Time_PC_PRE<<"\t" <<TLB_MC3E<<"\t" <<Time_MC3E<<"\t" <<Time_MC3E_PRE<<"\t"
                   <<  TLB_MV_P<<"\t"<<Time_MV_P<<"\t"<<TLB_MV_WEB<<"\t"<<Time_MV_WEB<<"\n";
        std::cout <<name<<" "<<"w="<<r_num[w]<<" "<< "Results have been written to output.txt" << std::endl;
        outputFile.close();
    }


    outputFile.open("..//output.txt",std::ios::app);
    outputFile << "\n";
    outputFile.close();
}

//
// //to get the hot map
// void TLB(vector<vector<vector<double>>>& Test, vector<vector<vector<double>>>& Train,string name) {
//     int n1 = Test.size();
//     int n2 = Train.size();
//     int m = Test[0].size();
//     int K = Test[0][0].size();
//
//
//     vector<vector<double>> new_test;
//     vector<vector<double>> new_train;
//     int new_m;
//     int new_k;
//     int new_r;
//     //
//     // vector<vector<double>> new_test_mv;
//     // vector<vector<double>> new_train_mv;
//     // int new_m_mv;
//     // int new_k_mv;
//     // int new_r_mv;
//
//
//     vector<double> r_num = {0.01,0.02,0.05,0.1,0.2};
//
//     clock_t start;
//     clock_t end;
//     clock_t lb_mv_time;
//     clock_t lb_mv_time_pre;
//
//     std::ofstream outputFile;
//     std::ofstream outputFile111;
//
//     double TLB_MV = 0;
//     double TLB_MV_DOUBLE = 0;
//     double TLB_GRAPH = 0;
//     double TLB_MC3E = 0;
//     double TLB_MV_P = 0;
//     double TLB_MV_WEB = 0;
//
//
//     double Time_DTW = 0;
//     double Time_MV = 0;
//     double Time_MC3E = 0;
//     double Time_MC3E_PRE = 0;
//     double Time_D_MV = 0;
//     double Time_GRAPH = 0;
//     double Time_PC = 0;
//     double Time_GRAPH_PRE = 0;
//     double Time_PC_PRE = 0;
//     double Time_MV_P = 0;
//     double Time_MV_WEB = 0;
//
//
//     vector<double> lb_mv_vector(n2);
//     vector<double> lb_mv_double_vector(n2);
//
//     vector<double> mu_Q(m);
//     vector<double> sigma_Q(m);
//
//     vector<double> mu_A(m);
//     vector<double> sigma_A(m);
//     vector<double> U_mu_QE(m), L_mu_QE(m), U_sigma_QE(m), L_sigma_QE(m);
//     vector<vector<double>> result,resulte;
//     // vector<vector<pair<double,int>>> result,resulte;
//
//     vector<double> intervalLength_muE(m);
//     vector<double> intervalLength_sigmaE(m);
//     vector<int> num(m);
//
//
//     vector<double> Q_2dim_0_mc3e(m);
//     vector<double> Q_2dim_1_mc3e(m);
//
//     vector<vector<double>>Q_L(K,vector<double>(m)),Q_U(K,vector<double>(m));
//
//     vector<vector<double>>L_UQ(K,vector<double>(m)),U_LQ(K,vector<double>(m));
//
//
//
//     for (int w = 0; w < r_num.size(); w++) {
//         int r = ceil(r_num[w]*m);
//         TLB_MV = 0;
//         TLB_MV_DOUBLE = 0;
//         TLB_GRAPH = 0;
//         TLB_MC3E = 0;
//         TLB_MV_P = 0;
//         TLB_MV_WEB = 0;
//
//
//         //
//         // //一起计算各个下界的TLB------------------------------------------------------------
//         for (int i = 0; i <n1; ++i) {
//
//             mu_sigma(Test[i],mu_Q,sigma_Q);
//             vector<vector<vector<Cell>>> T(m);
//
//             new_test = transpose(Test[i]);
//             new_m = new_test.size();
//             new_k = new_test[0].size();
//             new_r = ceil(r_num[w]*new_k);
//
//
//             for (int j = 0; j < n2; ++j) {
//                 new_train = transpose(Train[j]);
//                 if (MV) {
//                     // new_train_mv = transpose(Train[j]);
//                     LB_MV_T_hot_map(Q_L,Q_U, new_train, new_k, new_m);
//                 }
//             }
//         }
//     }
// }
