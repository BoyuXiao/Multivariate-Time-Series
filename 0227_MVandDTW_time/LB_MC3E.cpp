
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <utility>
#include "funcs.h"
using namespace std;

// 优化版本
  pair<double,double> initializeT(std::vector<vector<Cell>>& T, std::pair<vector<double>,vector<double>>& Q, double L_mu_Q, double U_mu_Q, double L_sigma_Q, double U_sigma_Q, int lambda, int n) {
      int numIntervals = lambda * static_cast<int>(std::sqrt(n));
      double intervalLength_mu = (U_mu_Q - L_mu_Q) / numIntervals;
      double intervalLength_sigma = (U_sigma_Q - L_sigma_Q) / numIntervals;
      T.resize(numIntervals+1);
      for (int i = 0; i < numIntervals+1; ++i) {
          double sigmaBottom = L_sigma_Q + i * intervalLength_sigma;
          double sigmaTop = sigmaBottom + intervalLength_sigma;
          for (int j = 0; j < numIntervals+1; ++j) {
              double muLeft = L_mu_Q + j * intervalLength_mu;
              double muRight = muLeft + intervalLength_mu;
              T[i].emplace_back(muLeft, muRight, sigmaBottom, sigmaTop);
          }
      }
      T[0][0].numIntervals_sigma = numIntervals;
      T[0][0].numIntervals_mu = numIntervals;
      if (intervalLength_mu!=0)intervalLength_mu = 1/intervalLength_mu;
      if (intervalLength_sigma!=0)intervalLength_sigma = 1/intervalLength_sigma;
      // if (intervalLength_mu==0)intervalLength_mu = 1;
      // if (intervalLength_sigma==0)intervalLength_sigma = 1;
      for (int i = 0; i < n; ++i) {
          double mu = Q.first[i];
          double sigma = Q.second[i];
          int muIndex ;
          int sigmaIndex ;
          if (mu!=U_mu_Q) muIndex = static_cast<int>((mu - L_mu_Q) * intervalLength_mu);
          else muIndex = numIntervals -1;
          // if (mu!=U_mu_Q) muIndex = static_cast<int>((mu - L_mu_Q) / intervalLength_mu);

          if (sigma!=U_sigma_Q) sigmaIndex = static_cast<int>((sigma - L_sigma_Q) * intervalLength_sigma);
          else sigmaIndex = numIntervals -1;
          // if (sigma!=U_sigma_Q) sigmaIndex = static_cast<int>((sigma - L_sigma_Q) / intervalLength_sigma);


          T[sigmaIndex][muIndex].dataList.emplace_back(mu,sigma);
     }

      //leftbottom point and left,bottom edge
      for (int i = 0; i < numIntervals; ++i) {
          for (int j = 0; j < numIntervals; ++j) {
              double mu,sigma,dist_X;
              //继承两个点
              int row = i;
              int col = MAX(j-1,0);
              // if (T[row][col].left_bottom_nearest_dist!=MAX_VAL) {
                  mu = T[row][col].left_bottom_nearest_mu;
                  sigma = T[row][col].left_bottom_nearest_sigma;
                  dist_X = sqrt(dist(T[i][j].Lleft_bottom.first,mu)+dist(T[i][j].Lleft_bottom.second,sigma));
                  T[i][j].left_bottom_nearest_dist = dist_X;
                  T[i][j].left_bottom_nearest_mu = mu;
                  T[i][j].left_bottom_nearest_sigma = sigma;
              // }
              row = MAX(i-1,0);
              col = j;
              // if (T[row][col].left_bottom_nearest_dist!=MAX_VAL) {
                  mu = T[row][col].left_bottom_nearest_mu;
                  sigma = T[row][col].left_bottom_nearest_sigma;
                  dist_X = sqrt(dist(T[i][j].Lleft_bottom.first,mu)+dist(T[i][j].Lleft_bottom.second,sigma));
                  if (dist_X < T[i][j].left_bottom_nearest_dist) {
                      T[i][j].left_bottom_nearest_dist = dist_X;
                      T[i][j].left_bottom_nearest_mu = mu;
                      T[i][j].left_bottom_nearest_sigma = sigma;
                  }
              // }
              row = MAX(i-1,0);
              col = MAX(j-1,0);
              //遍历一个格子来给点赋值
              for (int k=0;k<T[row][col].dataList.size();k++) {
                  if (i<1||j<1) break;
                  mu = T[row][col].dataList[k].first;
                  sigma = T[row][col].dataList[k].second;
                  dist_X = sqrt(dist(T[i][j].Lleft_bottom.first,mu)+dist(T[i][j].Lleft_bottom.second,sigma));
                  if (dist_X < T[i][j].left_bottom_nearest_dist) {
                      T[i][j].left_bottom_nearest_dist = dist_X;
                      T[i][j].left_bottom_nearest_mu = mu;
                      T[i][j].left_bottom_nearest_sigma = sigma;
                  }
              }
              //遍历一个格子来给left edge赋值
              row = i;
              col = MAX(0,j-1);
              if (T[row][col].dataList.size()!=0) {
                  for (int k=0;k<T[row][col].dataList.size();k++) {
                      if (j<1)break;
                      mu = T[row][col].dataList[k].first;
                      sigma = T[row][col].dataList[k].second;
                      dist_X = T[i][j].Lleft_bottom.first - mu;
                      if (dist_X < T[i][j].minDistToLleft) {
                          T[i][j].minDistToLleft = dist_X;
                          T[i][j].left_nearest_mu = mu;
                          T[i][j].left_nearest_sigma = sigma;
                      }
                  }
              }else {
                  // if (T[row][col].minDistToLleft!=MAX_VAL) {
                      mu = T[row][col].left_nearest_mu;
                      sigma = T[row][col].left_nearest_sigma;
                      dist_X = T[i][j].Lleft_bottom.first - mu;
                      T[i][j].minDistToLleft = dist_X;
                      T[i][j].left_nearest_mu = mu;
                      T[i][j].left_nearest_sigma = sigma;
                  // }
              }
              row = MAX(i-1,0);
              col = j;
              //遍历一个格子来给bottom edge赋值
              if (T[row][col].dataList.size()!=0) {
                  for (int k=0;k<T[row][col].dataList.size();k++) {
                      if (i<1)break;
                      mu = T[row][col].dataList[k].first;
                      sigma = T[row][col].dataList[k].second;
                      dist_X = T[i][j].Lright_bottom.second - sigma;
                      if (dist_X < T[i][j].minDistToLbottom) {
                          T[i][j].minDistToLbottom = dist_X;
                          T[i][j].bottom_nearest_mu = mu;
                          T[i][j].bottom_nearest_sigma = sigma;
                      }
                  }
              }else {
                  // if (T[row][col].minDistToLbottom!=MAX_VAL) {
                      mu = T[row][col].bottom_nearest_mu;
                      sigma = T[row][col].bottom_nearest_sigma;
                      dist_X = T[i][j].Lright_bottom.second - sigma;
                      T[i][j].minDistToLbottom = dist_X;
                      T[i][j].bottom_nearest_mu = mu;
                      T[i][j].bottom_nearest_sigma = sigma;
                  // }
              }
          }
      }

      //righttop point and right,top edge
      for (int i = numIntervals - 1; i >= 0; --i) {
          for (int j = numIntervals - 1; j >= 0; --j) {
              double mu,sigma,dist_X;
              //继承两个点
              int row = MIN(i+1,numIntervals-1);
              int col  = j;
              // if (T[row][col].right_top_nearest_dist!=MAX_VAL) {
                  mu = T[row][col].right_top_nearest_mu;
                  sigma = T[row][col].right_top_nearest_sigma;
                  dist_X = sqrt(dist(T[i][j].Lright_top.first,mu)+dist(T[i][j].Lright_top.second,sigma));
                  T[i][j].right_top_nearest_dist = dist_X;
                  T[i][j].right_top_nearest_mu = mu;
                  T[i][j].right_top_nearest_sigma = sigma;
              // }
              row = i;
              col = MIN(j+1,numIntervals-1);
              // if (T[row][col].right_top_nearest_dist!=MAX_VAL) {
                  mu = T[row][col].right_top_nearest_mu;
                  sigma = T[row][col].right_top_nearest_sigma;
                  dist_X = sqrt(dist(T[i][j].Lright_top.first,mu)+dist(T[i][j].Lright_top.second,sigma));
                  if (dist_X < T[i][j].right_top_nearest_dist) {
                      T[i][j].right_top_nearest_dist = dist_X;
                      T[i][j].right_top_nearest_mu = mu;
                      T[i][j].right_top_nearest_sigma = sigma;
                  }
              // }
              row = MIN(i+1,numIntervals-1);
              col = MIN(j+1,numIntervals-1);
              //遍历一个格子来给点赋值
              for (int k=0;k<T[row][col].dataList.size();k++) {
                  if (i+2>numIntervals||j+2>numIntervals)break;
                  mu = T[row][col].dataList[k].first;
                  sigma = T[row][col].dataList[k].second;
                  dist_X = sqrt(dist(T[i][j].Lright_top.first,mu)+dist(T[i][j].Lright_top.second,sigma));
                  if (dist_X < T[i][j].right_top_nearest_dist) {
                      T[i][j].right_top_nearest_dist = dist_X;
                      T[i][j].right_top_nearest_mu = mu;
                      T[i][j].right_top_nearest_sigma = sigma;
                  }
              }
              //遍历一个格子来给right edge赋值
              row = i;
              col = MIN(j+1,numIntervals - 1);
              if (T[row][col].dataList.size()!=0) {
                  for (int k=0;k<T[row][col].dataList.size();k++) {
                      if (j+2>numIntervals)break;
                      mu = T[row][col].dataList[k].first;
                      sigma = T[row][col].dataList[k].second;
                      dist_X = mu - T[i][j].Lright_bottom.first;
                      if (dist_X < T[i][j].minDistToLright) {
                          T[i][j].minDistToLright = dist_X;
                          T[i][j].right_nearest_mu = mu;
                          T[i][j].right_nearest_sigma = sigma;
                      }
                  }
              }else {
                  // if (T[row][col].minDistToLright!=MAX_VAL) {
                      mu = T[row][col].right_nearest_mu;
                      sigma = T[row][col].right_nearest_sigma;
                      dist_X = mu - T[i][j].Lright_bottom.first;
                      T[i][j].minDistToLright = dist_X;
                      T[i][j].right_nearest_mu = mu;
                      T[i][j].right_nearest_sigma = sigma;
                  // }
              }
              row = MIN(i+1,numIntervals-1);
              col = j;
              //遍历一个格子来给top edge赋值
              if (T[row][col].dataList.size()!=0) {
                  for (int k=0;k<T[row][col].dataList.size();k++) {
                      if (i+2>numIntervals)break;
                      mu = T[row][col].dataList[k].first;
                      sigma = T[row][col].dataList[k].second;
                      dist_X = sigma - T[i][j].Lright_top.second;
                      if (dist_X < T[i][j].minDistToLtop) {
                          T[i][j].minDistToLtop = dist_X;
                          T[i][j].top_nearest_mu = mu;
                          T[i][j].top_nearest_sigma = sigma;
                      }
                  }
              }else {
                  // if (T[row][col].minDistToLtop!=MAX_VAL) {
                      mu = T[row][col].top_nearest_mu;
                      sigma = T[row][col].top_nearest_sigma;
                      dist_X = sigma - T[i][j].Lright_top.second;
                      T[i][j].minDistToLtop = dist_X;
                      T[i][j].top_nearest_mu = mu;
                      T[i][j].top_nearest_sigma = sigma;
                  // }
              }
          }
      }


       //rightbottom point
       for (int i = 0; i < numIntervals; ++i) {
           for (int j = numIntervals - 1; j >= 0; --j) {
               double mu,sigma,dist_X;
               //继承两个点
               int row = i;
               int col = MIN(j+1,numIntervals - 1);
               // if (T[row][col].right_bottom_nearest_dist!=MAX_VAL) {
                   mu = T[row][col].right_bottom_nearest_mu;
                   sigma = T[row][col].right_bottom_nearest_sigma;
                   dist_X = sqrt(dist(T[i][j].Lright_bottom.first,mu)+dist(T[i][j].Lright_bottom.second,sigma));
                   T[i][j].right_bottom_nearest_dist = dist_X;
                   T[i][j].right_bottom_nearest_mu = mu;
                   T[i][j].right_bottom_nearest_sigma = sigma;
               // }
               row = MAX(i-1,0);
               col = j;
               // if (T[row][col].right_bottom_nearest_dist!=MAX_VAL) {
                   mu = T[row][col].right_bottom_nearest_mu;
                   sigma = T[row][col].right_bottom_nearest_sigma;
                   dist_X = sqrt(dist(T[i][j].Lright_bottom.first,mu)+dist(T[i][j].Lright_bottom.second,sigma));
                   if (dist_X < T[i][j].right_bottom_nearest_dist) {
                       T[i][j].right_bottom_nearest_dist = dist_X;
                       T[i][j].right_bottom_nearest_mu = mu;
                       T[i][j].right_bottom_nearest_sigma = sigma;
                   }
               // }
               row = MAX(i-1,0);
               col = MIN(j+1,numIntervals - 1);
               //遍历一个格子来给点赋值
               for (int k=0;k<T[row][col].dataList.size();k++) {
                   if (i<1||j+2>numIntervals)break;
                   mu = T[row][col].dataList[k].first;
                   sigma = T[row][col].dataList[k].second;
                   dist_X = sqrt(dist(T[i][j].Lright_bottom.first,mu)+dist(T[i][j].Lright_bottom.second,sigma));
                   if (dist_X < T[i][j].right_bottom_nearest_dist) {
                       T[i][j].right_bottom_nearest_dist = dist_X;
                       T[i][j].right_bottom_nearest_mu = mu;
                       T[i][j].right_bottom_nearest_sigma = sigma;
                   }
               }
           }
       }

       //lefttop point
       for (int i = numIntervals - 1; i >= 0; --i) {
           for (int j = 0; j < numIntervals; ++j) {
               double mu,sigma,dist_X;
               //继承两个点
               int row = i;
               int col = MAX(j-1,0);
               // if (T[row][col].left_top_nearest_dist!=MAX_VAL) {
                   mu = T[row][col].left_top_nearest_mu;
                   sigma = T[row][col].left_top_nearest_sigma;
                   dist_X = sqrt(dist(T[i][j].Lleft_top.first,mu)+dist(T[i][j].Lleft_top.second,sigma));
                   T[i][j].left_top_nearest_dist = dist_X;
                   T[i][j].left_top_nearest_mu = mu;
                   T[i][j].left_top_nearest_sigma = sigma;
               // }
               row = MIN(i+1,numIntervals-1);
               col = j;
               // if (T[row][col].left_top_nearest_dist!=MAX_VAL) {
                   mu = T[row][col].left_top_nearest_mu;
                   sigma = T[row][col].left_top_nearest_sigma;
                   dist_X = sqrt(dist(T[i][j].Lleft_top.first,mu)+dist(T[i][j].Lleft_top.second,sigma));
                   if (dist_X < T[i][j].left_top_nearest_dist) {
                       T[i][j].left_top_nearest_dist = dist_X;
                       T[i][j].left_top_nearest_mu = mu;
                       T[i][j].left_top_nearest_sigma = sigma;
                   }
               // }
               row = MIN(i+1,numIntervals-1);
               col = MAX(j-1,0);
               //遍历一个格子来给点赋值
               for (int k=0;k<T[row][col].dataList.size();k++) {
                   if (j<1||i+2>numIntervals)break;
                   mu = T[row][col].dataList[k].first;
                   sigma = T[row][col].dataList[k].second;
                   dist_X = sqrt(dist(T[i][j].Lleft_top.first,mu)+dist(T[i][j].Lleft_top.second,sigma));
                   if (dist_X < T[i][j].left_top_nearest_dist) {
                       T[i][j].left_top_nearest_dist = dist_X;
                       T[i][j].left_top_nearest_mu = mu;
                       T[i][j].left_top_nearest_sigma = sigma;
                   }
               }
           }
       }

       //update
       for (int i = 0; i < numIntervals; ++i) {
           for (int j = 0; j < numIntervals; ++j) {
               T[i][j].minDistToLleft = MIN(MIN(T[i][j].minDistToLleft,T[i][j].left_bottom_nearest_dist),T[i][j].left_top_nearest_dist);
               T[i][j].minDistToLright = MIN(MIN(T[i][j].minDistToLright,T[i][j].right_bottom_nearest_dist),T[i][j].right_top_nearest_dist);
               T[i][j].minDistToLbottom = MIN(MIN(T[i][j].minDistToLbottom,T[i][j].left_bottom_nearest_dist),T[i][j].right_bottom_nearest_dist);
               T[i][j].minDistToLtop = MIN(MIN(T[i][j].minDistToLtop,T[i][j].left_top_nearest_dist),T[i][j].right_top_nearest_dist);
           }
       }

      //copy
      for (int i = 0; i < numIntervals; ++i) {
          T[i][numIntervals]=T[i][numIntervals-1];
      }
      for (int j = 0; j < numIntervals; ++j) {
          T[numIntervals][j]=T[numIntervals-1][j];
      }
      T[numIntervals][numIntervals]=T[numIntervals-1][numIntervals-1];

       return {intervalLength_mu, intervalLength_sigma};
   }
//
// //

pair<double,double> initializeT_divide(std::vector<vector<Cell>>& T, std::pair<vector<double>,vector<double>>& Q, double L_mu_Q, double U_mu_Q, double L_sigma_Q, double U_sigma_Q, int lambda, int n) {
    int numIntervals = 2 * lambda * static_cast<int>(std::sqrt(n));
    int numIntervals_mu ;
    int numIntervals_sigma ;
    double range_mu = U_mu_Q - L_mu_Q;
    double range_sigma = U_sigma_Q - L_sigma_Q;
    if (range_mu == 0 && range_sigma == 0) {
        numIntervals_mu = 1;
        numIntervals_sigma = 1;
    }
    else if (range_mu > range_sigma) {
        numIntervals_mu = static_cast<int>(numIntervals * (range_mu / (range_mu + range_sigma)));
        numIntervals_sigma = numIntervals - numIntervals_mu;
    } else {
        numIntervals_sigma = static_cast<int>(numIntervals * (range_sigma / (range_mu + range_sigma)));
        numIntervals_mu = numIntervals - numIntervals_sigma;
    }
    double intervalLength_mu = (U_mu_Q - L_mu_Q) / numIntervals_mu;
    double intervalLength_sigma = (U_sigma_Q - L_sigma_Q) / numIntervals_sigma;
    // cout<<range_mu<<" "<<range_sigma<<""<<numIntervals_sigma<<" "<<numIntervals_mu<<" "<<intervalLength_sigma<<" "<<intervalLength_mu<<endl;
    T.resize(numIntervals_sigma+1);
    for (int i = 0; i < numIntervals_sigma+1; ++i) {
        double sigmaBottom = L_sigma_Q + i * intervalLength_sigma;
        double sigmaTop = sigmaBottom + intervalLength_sigma;
        for (int j = 0; j < numIntervals_mu+1; ++j) {
            double muLeft = L_mu_Q + j * intervalLength_mu;
            double muRight = muLeft + intervalLength_mu;
            T[i].emplace_back(muLeft, muRight, sigmaBottom, sigmaTop);
        }
    }
      T[0][0].numIntervals_sigma = numIntervals_sigma;
      T[0][0].numIntervals_mu = numIntervals_mu;
      if (intervalLength_mu!=0)intervalLength_mu = 1/intervalLength_mu;
      if (intervalLength_sigma!=0)intervalLength_sigma = 1/intervalLength_sigma;

    for (int i = 0; i < n; ++i) {
        double mu = Q.first[i];
        double sigma = Q.second[i];
        int muIndex;
        int sigmaIndex;
        if (mu!=U_mu_Q) muIndex = static_cast<int>((mu - L_mu_Q) * intervalLength_mu);
        else muIndex = numIntervals_mu -1;
        if (sigma!=U_sigma_Q) sigmaIndex = static_cast<int>((sigma - L_sigma_Q) * intervalLength_sigma);
        else sigmaIndex = numIntervals_sigma -1;
        // if (muIndex >= 0 && muIndex < numIntervals_mu && sigmaIndex >= 0 && sigmaIndex < numIntervals_sigma) {
            T[sigmaIndex][muIndex].dataList.emplace_back(mu,sigma);
        // }
    }
    //leftbottom point and left,bottom edge
      for (int i = 0; i < numIntervals_sigma; ++i) {
          for (int j = 0; j < numIntervals_mu; ++j) {
              double mu,sigma,dist_X;
              //继承两个点
              int row = i;
              int col = MAX(j-1,0);
              // if (T[row][col].left_bottom_nearest_dist!=MAX_VAL) {
                  mu = T[row][col].left_bottom_nearest_mu;
                  sigma = T[row][col].left_bottom_nearest_sigma;
                  dist_X = sqrt(dist(T[i][j].Lleft_bottom.first,mu)+dist(T[i][j].Lleft_bottom.second,sigma));
                  T[i][j].left_bottom_nearest_dist = dist_X;
                  T[i][j].left_bottom_nearest_mu = mu;
                  T[i][j].left_bottom_nearest_sigma = sigma;
              // }
              row = MAX(i-1,0);
              col = j;
              // if (T[row][col].left_bottom_nearest_dist!=MAX_VAL) {
                  mu = T[row][col].left_bottom_nearest_mu;
                  sigma = T[row][col].left_bottom_nearest_sigma;
                  dist_X = sqrt(dist(T[i][j].Lleft_bottom.first,mu)+dist(T[i][j].Lleft_bottom.second,sigma));
                  if (dist_X < T[i][j].left_bottom_nearest_dist) {
                      T[i][j].left_bottom_nearest_dist = dist_X;
                      T[i][j].left_bottom_nearest_mu = mu;
                      T[i][j].left_bottom_nearest_sigma = sigma;
                  }
              // }
              row = MAX(i-1,0);
              col = MAX(j-1,0);
              //遍历一个格子来给点赋值
              for (int k=0;k<T[row][col].dataList.size();k++) {
                  if (i<1||j<1) break;
                  mu = T[row][col].dataList[k].first;
                  sigma = T[row][col].dataList[k].second;
                  dist_X = sqrt(dist(T[i][j].Lleft_bottom.first,mu)+dist(T[i][j].Lleft_bottom.second,sigma));
                  if (dist_X < T[i][j].left_bottom_nearest_dist) {
                      T[i][j].left_bottom_nearest_dist = dist_X;
                      T[i][j].left_bottom_nearest_mu = mu;
                      T[i][j].left_bottom_nearest_sigma = sigma;
                  }
              }
              //遍历一个格子来给left edge赋值
              row = i;
              col = MAX(0,j-1);
              if (T[row][col].dataList.size()!=0) {
                  for (int k=0;k<T[row][col].dataList.size();k++) {
                      if (j<1)break;
                      mu = T[row][col].dataList[k].first;
                      sigma = T[row][col].dataList[k].second;
                      dist_X = T[i][j].Lleft_bottom.first - mu;
                      if (dist_X < T[i][j].minDistToLleft) {
                          T[i][j].minDistToLleft = dist_X;
                          T[i][j].left_nearest_mu = mu;
                          T[i][j].left_nearest_sigma = sigma;
                      }
                  }
              }else {
                  // if (T[row][col].minDistToLleft!=MAX_VAL) {
                      mu = T[row][col].left_nearest_mu;
                      sigma = T[row][col].left_nearest_sigma;
                      dist_X = T[i][j].Lleft_bottom.first - mu;
                      T[i][j].minDistToLleft = dist_X;
                      T[i][j].left_nearest_mu = mu;
                      T[i][j].left_nearest_sigma = sigma;
                  // }
              }
              row = MAX(i-1,0);
              col = j;
              //遍历一个格子来给bottom edge赋值
              if (T[row][col].dataList.size()!=0) {
                  for (int k=0;k<T[row][col].dataList.size();k++) {
                      if (i<1)break;
                      mu = T[row][col].dataList[k].first;
                      sigma = T[row][col].dataList[k].second;
                      dist_X = T[i][j].Lright_bottom.second - sigma;
                      if (dist_X < T[i][j].minDistToLbottom) {
                          T[i][j].minDistToLbottom = dist_X;
                          T[i][j].bottom_nearest_mu = mu;
                          T[i][j].bottom_nearest_sigma = sigma;
                      }
                  }
              }else {
                  // if (T[row][col].minDistToLbottom!=MAX_VAL) {
                      mu = T[row][col].bottom_nearest_mu;
                      sigma = T[row][col].bottom_nearest_sigma;
                      dist_X = T[i][j].Lright_bottom.second - sigma;
                      T[i][j].minDistToLbottom = dist_X;
                      T[i][j].bottom_nearest_mu = mu;
                      T[i][j].bottom_nearest_sigma = sigma;
                  // }
              }
          }
      }

      //righttop point and right,top edge
      for (int i = numIntervals_sigma - 1; i >= 0; --i) {
          for (int j = numIntervals_mu - 1; j >= 0; --j) {
              double mu,sigma,dist_X;
              //继承两个点
              int row = MIN(i+1,numIntervals_sigma-1);
              int col  = j;
              // if (T[row][col].right_top_nearest_dist!=MAX_VAL) {
                  mu = T[row][col].right_top_nearest_mu;
                  sigma = T[row][col].right_top_nearest_sigma;
                  dist_X = sqrt(dist(T[i][j].Lright_top.first,mu)+dist(T[i][j].Lright_top.second,sigma));
                  T[i][j].right_top_nearest_dist = dist_X;
                  T[i][j].right_top_nearest_mu = mu;
                  T[i][j].right_top_nearest_sigma = sigma;
              // }
              row = i;
              col = MIN(j+1,numIntervals_mu-1);
              // if (T[row][col].right_top_nearest_dist!=MAX_VAL) {
                  mu = T[row][col].right_top_nearest_mu;
                  sigma = T[row][col].right_top_nearest_sigma;
                  dist_X = sqrt(dist(T[i][j].Lright_top.first,mu)+dist(T[i][j].Lright_top.second,sigma));
                  if (dist_X < T[i][j].right_top_nearest_dist) {
                      T[i][j].right_top_nearest_dist = dist_X;
                      T[i][j].right_top_nearest_mu = mu;
                      T[i][j].right_top_nearest_sigma = sigma;
                  }
              // }
              row = MIN(i+1,numIntervals_sigma-1);
              col = MIN(j+1,numIntervals_mu-1);
              //遍历一个格子来给点赋值
              for (int k=0;k<T[row][col].dataList.size();k++) {
                  if (i+2>numIntervals_sigma||j+2>numIntervals_mu)break;
                  mu = T[row][col].dataList[k].first;
                  sigma = T[row][col].dataList[k].second;
                  dist_X = sqrt(dist(T[i][j].Lright_top.first,mu)+dist(T[i][j].Lright_top.second,sigma));
                  if (dist_X < T[i][j].right_top_nearest_dist) {
                      T[i][j].right_top_nearest_dist = dist_X;
                      T[i][j].right_top_nearest_mu = mu;
                      T[i][j].right_top_nearest_sigma = sigma;
                  }
              }
              //遍历一个格子来给right edge赋值
              row = i;
              col = MIN(j+1,numIntervals_mu - 1);
              if (T[row][col].dataList.size()!=0) {
                  for (int k=0;k<T[row][col].dataList.size();k++) {
                      if (j+2>numIntervals_mu)break;
                      mu = T[row][col].dataList[k].first;
                      sigma = T[row][col].dataList[k].second;
                      dist_X = mu - T[i][j].Lright_bottom.first;
                      if (dist_X < T[i][j].minDistToLright) {
                          T[i][j].minDistToLright = dist_X;
                          T[i][j].right_nearest_mu = mu;
                          T[i][j].right_nearest_sigma = sigma;
                      }
                  }
              }else {
                  // if (T[row][col].minDistToLright!=MAX_VAL) {
                      mu = T[row][col].right_nearest_mu;
                      sigma = T[row][col].right_nearest_sigma;
                      dist_X = mu - T[i][j].Lright_bottom.first;
                      T[i][j].minDistToLright = dist_X;
                      T[i][j].right_nearest_mu = mu;
                      T[i][j].right_nearest_sigma = sigma;
                  // }
              }
              row = MIN(i+1,numIntervals_sigma-1);
              col = j;
              //遍历一个格子来给top edge赋值
              if (T[row][col].dataList.size()!=0) {
                  for (int k=0;k<T[row][col].dataList.size();k++) {
                      if (i+2>numIntervals_sigma)break;
                      mu = T[row][col].dataList[k].first;
                      sigma = T[row][col].dataList[k].second;
                      dist_X = sigma - T[i][j].Lright_top.second;
                      if (dist_X < T[i][j].minDistToLtop) {
                          T[i][j].minDistToLtop = dist_X;
                          T[i][j].top_nearest_mu = mu;
                          T[i][j].top_nearest_sigma = sigma;
                      }
                  }
              }else {
                  // if (T[row][col].minDistToLtop!=MAX_VAL) {
                      mu = T[row][col].top_nearest_mu;
                      sigma = T[row][col].top_nearest_sigma;
                      dist_X = sigma - T[i][j].Lright_top.second;
                      T[i][j].minDistToLtop = dist_X;
                      T[i][j].top_nearest_mu = mu;
                      T[i][j].top_nearest_sigma = sigma;
                  // }
              }
          }
      }


       //rightbottom point
       for (int i = 0; i < numIntervals_sigma; ++i) {
           for (int j = numIntervals_mu - 1; j >= 0; --j) {
               double mu,sigma,dist_X;
               //继承两个点
               int row = i;
               int col = MIN(j+1,numIntervals_mu - 1);
               // if (T[row][col].right_bottom_nearest_dist!=MAX_VAL) {
                   mu = T[row][col].right_bottom_nearest_mu;
                   sigma = T[row][col].right_bottom_nearest_sigma;
                   dist_X = sqrt(dist(T[i][j].Lright_bottom.first,mu)+dist(T[i][j].Lright_bottom.second,sigma));
                   T[i][j].right_bottom_nearest_dist = dist_X;
                   T[i][j].right_bottom_nearest_mu = mu;
                   T[i][j].right_bottom_nearest_sigma = sigma;
               // }
               row = MAX(i-1,0);
               col = j;
               // if (T[row][col].right_bottom_nearest_dist!=MAX_VAL) {
                   mu = T[row][col].right_bottom_nearest_mu;
                   sigma = T[row][col].right_bottom_nearest_sigma;
                   dist_X = sqrt(dist(T[i][j].Lright_bottom.first,mu)+dist(T[i][j].Lright_bottom.second,sigma));
                   if (dist_X < T[i][j].right_bottom_nearest_dist) {
                       T[i][j].right_bottom_nearest_dist = dist_X;
                       T[i][j].right_bottom_nearest_mu = mu;
                       T[i][j].right_bottom_nearest_sigma = sigma;
                   }
               // }
               row = MAX(i-1,0);
               col = MIN(j+1,numIntervals_mu - 1);
               //遍历一个格子来给点赋值
               for (int k=0;k<T[row][col].dataList.size();k++) {
                   if (i<1||j+2>numIntervals_mu)break;
                   mu = T[row][col].dataList[k].first;
                   sigma = T[row][col].dataList[k].second;
                   dist_X = sqrt(dist(T[i][j].Lright_bottom.first,mu)+dist(T[i][j].Lright_bottom.second,sigma));
                   if (dist_X < T[i][j].right_bottom_nearest_dist) {
                       T[i][j].right_bottom_nearest_dist = dist_X;
                       T[i][j].right_bottom_nearest_mu = mu;
                       T[i][j].right_bottom_nearest_sigma = sigma;
                   }
               }
           }
       }

       //lefttop point
       for (int i = numIntervals_sigma - 1; i >= 0; --i) {
           for (int j = 0; j < numIntervals_mu; ++j) {
               double mu,sigma,dist_X;
               //继承两个点
               int row = i;
               int col = MAX(j-1,0);
               // if (T[row][col].left_top_nearest_dist!=MAX_VAL) {
                   mu = T[row][col].left_top_nearest_mu;
                   sigma = T[row][col].left_top_nearest_sigma;
                   dist_X = sqrt(dist(T[i][j].Lleft_top.first,mu)+dist(T[i][j].Lleft_top.second,sigma));
                   T[i][j].left_top_nearest_dist = dist_X;
                   T[i][j].left_top_nearest_mu = mu;
                   T[i][j].left_top_nearest_sigma = sigma;
               // }
               row = MIN(i+1,numIntervals_sigma-1);
               col = j;
               // if (T[row][col].left_top_nearest_dist!=MAX_VAL) {
                   mu = T[row][col].left_top_nearest_mu;
                   sigma = T[row][col].left_top_nearest_sigma;
                   dist_X = sqrt(dist(T[i][j].Lleft_top.first,mu)+dist(T[i][j].Lleft_top.second,sigma));
                   if (dist_X < T[i][j].left_top_nearest_dist) {
                       T[i][j].left_top_nearest_dist = dist_X;
                       T[i][j].left_top_nearest_mu = mu;
                       T[i][j].left_top_nearest_sigma = sigma;
                   }
               // }
               row = MIN(i+1,numIntervals_sigma-1);
               col = MAX(j-1,0);
               //遍历一个格子来给点赋值
               for (int k=0;k<T[row][col].dataList.size();k++) {
                   if (j<1||i+2>numIntervals_sigma)break;
                   mu = T[row][col].dataList[k].first;
                   sigma = T[row][col].dataList[k].second;
                   dist_X = sqrt(dist(T[i][j].Lleft_top.first,mu)+dist(T[i][j].Lleft_top.second,sigma));
                   if (dist_X < T[i][j].left_top_nearest_dist) {
                       T[i][j].left_top_nearest_dist = dist_X;
                       T[i][j].left_top_nearest_mu = mu;
                       T[i][j].left_top_nearest_sigma = sigma;
                   }
               }
           }
       }

       //update
       for (int i = 0; i < numIntervals_sigma; ++i) {
           for (int j = 0; j < numIntervals_mu; ++j) {
               T[i][j].minDistToLleft = MIN(MIN(T[i][j].minDistToLleft,T[i][j].left_bottom_nearest_dist),T[i][j].left_top_nearest_dist);
               T[i][j].minDistToLright = MIN(MIN(T[i][j].minDistToLright,T[i][j].right_bottom_nearest_dist),T[i][j].right_top_nearest_dist);
               T[i][j].minDistToLbottom = MIN(MIN(T[i][j].minDistToLbottom,T[i][j].left_bottom_nearest_dist),T[i][j].right_bottom_nearest_dist);
               T[i][j].minDistToLtop = MIN(MIN(T[i][j].minDistToLtop,T[i][j].left_top_nearest_dist),T[i][j].right_top_nearest_dist);
           }
       }

      //copy
      for (int i = 0; i < numIntervals_sigma; ++i) {
          T[i][numIntervals_mu]=T[i][numIntervals_mu-1];
      }
      for (int j = 0; j < numIntervals_mu; ++j) {
          T[numIntervals_sigma][j]=T[numIntervals_sigma-1][j];
      }
      T[numIntervals_sigma][numIntervals_mu]=T[numIntervals_sigma-1][numIntervals_mu-1];
    return {intervalLength_mu, intervalLength_sigma};
}


double CornerDistE(const vector<double>& X, const vector<double>& V, double dNN) {
    double dMu = (X[0] - V[0])*(X[0] - V[0]);
    double dSigma = (X[1] - V[1])*(X[1] - V[1]);
    return dMu + dSigma + dNN*dNN + 2 * dNN * sqrt(MIN(dMu, dSigma));
}


double CenterDist2D(double mu_X, double sigma_X, double L_mu_Q, double U_mu_Q, double L_sigma_Q, double U_sigma_Q, std::vector<vector<Cell>>& T,double intervalLength_mu,double intervalLength_sigma,int lambda,int num) {
    double dist_x = 0.0;
     int col = round((mu_X - L_mu_Q) * intervalLength_mu);
     int row = round((sigma_X - L_sigma_Q) * intervalLength_sigma);
      // int col = round((mu_X - L_mu_Q) / intervalLength_mu);
      // int row = round((sigma_X - L_sigma_Q) / intervalLength_sigma);
    if (T[row][col].dataList.size() >= 2) {
        return dist_x;
    } else {
        double leftDist = dist(T[row][col].minDistToLleft+mu_X,T[row][col].Lleft_bottom.first);
        double rightDist = dist(T[row][col].minDistToLright+T[row][col].Lright_bottom.first,mu_X);
        double bottomDist = dist(T[row][col].minDistToLbottom+sigma_X,T[row][col].Lleft_bottom.second);
        double topDist = dist(T[row][col].minDistToLtop+T[row][col].Lleft_top.second,sigma_X);
        double epsilon = MIN(MIN(leftDist, rightDist),MIN(bottomDist,topDist));
        if (T[row][col].dataList.size()== 1)
            epsilon = MIN(epsilon, dist(mu_X,T[row][col].dataList[0].first)+dist(sigma_X,T[row][col].dataList[0].second));
        dist_x += epsilon;
    }
    return dist_x;
}


double EdgeDistE(const vector<double>& X, const vector<double>& V, double dNN,bool MiddleLeftOrRight,bool TopOrBottomMiddle,double L_mu_Q, double U_mu_Q,
    double L_sigma_Q, double U_sigma_Q,vector<vector<Cell>>& T,double intervalLength_mu,double intervalLength_sigma,int lambda,int num) {
    double dMu = (X[0] - V[0])*(X[0] - V[0]);
    double dSigma = (X[1] - V[1])*(X[1] - V[1]);
    if (MiddleLeftOrRight) {
        if (dNN*dNN <= dSigma) {
            return dMu + CenterDist2D(V[0],X[1],L_mu_Q,U_mu_Q,L_sigma_Q,U_sigma_Q,T,intervalLength_mu,intervalLength_sigma,lambda,num);
        } else {
            // return dNN*dNN + dMu + dSigma - 2 * dNN * sqrt(dSigma);
            return MAX(dMu + CenterDist2D(V[0],X[1],L_mu_Q,U_mu_Q,L_sigma_Q,U_sigma_Q,T,intervalLength_mu,intervalLength_sigma,lambda,num),dNN*dNN + dMu + dSigma - 2 * dNN * sqrt(dSigma));
        }
    } else if (TopOrBottomMiddle) {
        if (dNN*dNN <= dMu) {
            return dSigma + CenterDist2D(X[0],V[1],L_mu_Q,U_mu_Q,L_sigma_Q,U_sigma_Q,T,intervalLength_mu,intervalLength_sigma,lambda,num);
            // return dSigma+CenterDist(X[0],L_mu_Q,U_mu_Q,T_mu,intervalLength_mu,lambda,num);
        } else {
            // return dNN*dNN + dMu + dSigma - 2 * dNN * sqrt(dMu);
            return MAX(dSigma + CenterDist2D(X[0],V[1],L_mu_Q,U_mu_Q,L_sigma_Q,U_sigma_Q,T,intervalLength_mu,intervalLength_sigma,lambda,num),dNN*dNN + dMu + dSigma - 2 * dNN * sqrt(dMu));
        }
    }

    return -1;
}



vector<vector<double>> Preprocessing_MC3E(vector<vector<double>>& Q,
    vector<double>& mu_Q,vector<double>& sigma_Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,int r) {
    int lambda =3;
    int m = Q.size();
    lower_upper_lemire(mu_Q, m, r, L_mu_Q, U_mu_Q);
    lower_upper_lemire(sigma_Q, m, r, L_sigma_Q, U_sigma_Q);

    vector<vector<double>> result(4, vector<double>(m));


    double minDbl,minDbr,minDtl,minDtr;
    vector<int> start(m);
    vector<int> end(m);
    vector<double> Qt(2);
    vector<double> corner(2);
    double dist;

    for (int i = 0; i < m; ++i) {
        start[i] = MAX(0,i-r);
        end[i] = MIN(m-1,i+r);
    }
    //计算四个dNN预选值
    for (int i = 0; i < m; ++i) {
        minDbl = MAX_VAL;
        for (int t = start[i]; t <= end[i]; ++t) {
            Qt = {mu_Q[t], sigma_Q[t]};
            corner = {L_mu_Q[i], L_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDbl = MIN(minDbl, dist);
        }
        result[0][i] = minDbl;
    }
    for (int i = 0; i < m; ++i) {
        minDbr = MAX_VAL;
        for (int t = start[i]; t <= end[i]; ++t) {
            Qt = {mu_Q[t], sigma_Q[t]};
            corner = {U_mu_Q[i], L_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDbr = MIN(minDbr, dist);
        }
        result[1][i] = minDbr;
    }
    for (int i = 0; i < m; ++i) {
        minDtl = MAX_VAL;
        for (int t = start[i]; t <= end[i]; ++t) {
            Qt = {mu_Q[t], sigma_Q[t]};
            corner = {L_mu_Q[i], U_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDtl = MIN(minDtl, dist);
        }
        result[2][i] = minDtl;
    }
    for (int i = 0; i < m; ++i) {
        minDtr = MAX_VAL;
        for (int t = start[i]; t <= end[i]; ++t) {
            Qt = {mu_Q[t], sigma_Q[t]};
            corner = {U_mu_Q[i], U_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDtr = MIN(minDtr, dist);
        }
        result[3][i] = minDtr;
    }



    std::pair<double, double> intervalLength;
    //计算Q的m个T表
    for (int i = 0; i < m; ++i) {
        num[i] = end[i]-start[i]+1;
        vector<double> mu_Q_part(mu_Q.begin() + start[i], mu_Q.begin() + end[i] + 1);
        vector<double> sigma_Q_part(sigma_Q.begin() + start[i], sigma_Q.begin() + end[i] + 1);
        pair<vector<double>,vector<double>> Q_part = {mu_Q_part,sigma_Q_part};
        intervalLength = initializeT(T[i], Q_part, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i],U_sigma_Q[i],lambda, num[i]);
        intervalLength_mu[i]=intervalLength.first;
        intervalLength_sigma[i]=intervalLength.second;
    }

    for (int i = 0; i < result.size(); ++i) {
        for (int j = 0; j < result[i].size(); ++j) {
            result[i][j] = sqrt(result[i][j]);
        }
    }

    return result;
}

vector<vector<double>> Preprocessing_MC3E_PCA(vector<vector<double>>& Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,int r,
    vector<Eigen::MatrixXd>& projection_matrix_vector,vector<Eigen::VectorXd>& mean_Q_vector) {

    int lambda =3;
    int m = Q.size();
    vector<vector<double>> result(4, vector<double>(m));
    double minDbl,minDbr,minDtl,minDtr;
    vector<int> start(m);
    vector<int> end(m);
    vector<double> Qt(2);
    vector<double> corner(2);
    double dist;

    vector<vector<double>> mu_Q_pca(m);
    vector<vector<double>> sigma_Q_pca(m);
    vector<vector<double>> Q_pca_part;

    for (int i = 0; i < m; ++i) {
        start[i] = MAX(0,i-r);
        end[i] = MIN(m-1,i+r);
    }

    //把每一组的Q降维并记录下投影矩阵和均值
      std::pair<double, double> intervalLength;
      int dim = 2;
      //计算Q的m个T表
      for (int i = 0; i < m; ++i) {
          num[i] = end[i]-start[i]+1;
          //split
          vector<vector<double>> Q_original_part(Q.begin() + start[i], Q.begin() + end[i] + 1);
          //PCA
          projection_matrix_vector[i] = pca_reduce_with_projection(Q_original_part,dim,Q_pca_part,mean_Q_vector[i]);
          // debug------half
          // split_first_and_second(Q_original_part,Q_pca_part);
          // mu_sigma(Q_original_part,Q_original_part[0].size(),num[i],mu_Q_pca[i],sigma_Q_pca[i]);

          //mu and sigma
          // mu_sigma(Q_pca_part,dim,num[i],mu_Q_pca[i],sigma_Q_pca[i]);
          split_two_dim(Q_pca_part,mu_Q_pca[i],sigma_Q_pca[i]);
          //U and L
          lower_upper_lemire_pca(num[i],mu_Q_pca[i],sigma_Q_pca[i],U_mu_Q,L_mu_Q,U_sigma_Q,L_sigma_Q);
          //T
          pair<vector<double>,vector<double>> Q_part = {mu_Q_pca[i],sigma_Q_pca[i]};
          intervalLength = initializeT(T[i], Q_part, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i],U_sigma_Q[i],lambda, num[i]);
          intervalLength_mu[i]=intervalLength.first;
          intervalLength_sigma[i]=intervalLength.second;
      }


    //计算四个dNN预选值
    for (int i = 0; i < m; ++i) {
        minDbl = MAX_VAL;
        for (int t = 0; t < num[i]; ++t) {
            Qt = {mu_Q_pca[i][t], sigma_Q_pca[i][t]};
            corner = {L_mu_Q[i], L_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDbl = MIN(minDbl, dist);
        }
        result[0][i] = minDbl;
    }
    for (int i = 0; i < m; ++i) {
        minDbr = MAX_VAL;
        for (int t = 0; t < num[i]; ++t) {
            Qt = {mu_Q_pca[i][t], sigma_Q_pca[i][t]};
            corner = {U_mu_Q[i], L_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDbr = MIN(minDbr, dist);
        }
        result[1][i] = minDbr;
    }
    for (int i = 0; i < m; ++i) {
        minDtl = MAX_VAL;
        for (int t = 0; t < num[i]; ++t) {
            Qt = {mu_Q_pca[i][t], sigma_Q_pca[i][t]};
            corner = {L_mu_Q[i], U_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDtl = MIN(minDtl, dist);
        }
        result[2][i] = minDtl;
    }
    for (int i = 0; i < m; ++i) {
        minDtr = MAX_VAL;
        for (int t = 0; t < num[i]; ++t) {
            Qt = {mu_Q_pca[i][t], sigma_Q_pca[i][t]};
            corner = {U_mu_Q[i], U_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDtr = MIN(minDtr, dist);
        }
        result[3][i] = minDtr;
    }



    for (int i = 0; i < result.size(); ++i) {
        for (int j = 0; j < result[i].size(); ++j) {
            result[i][j] = sqrt(result[i][j]);
        }
    }

    return result;
}


vector<vector<double>> Preprocessing_MC3E_SORT_strength(vector<vector<double>>& Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,int r,
    vector<vector<int>> &split_indices_sort) {

    int lambda =3;
    int m = Q.size();
    vector<vector<double>> result(4, vector<double>(m));
    double minDbl,minDbr,minDtl,minDtr;
    vector<int> start(m);
    vector<int> end(m);
    vector<double> Qt(2);
    vector<double> corner(2);
    double dist;

    vector<vector<double>> mu_Q_sort1(m);
    vector<vector<double>> sigma_Q_sort1(m);
      vector<vector<double>> mu_Q_sort2(m);
      vector<vector<double>> sigma_Q_sort2(m);

    for (int i = 0; i < m; ++i) {
        start[i] = MAX(0,i-r);
        end[i] = MIN(m-1,i+r);
    }

    //把每一组的Q降维并记录下投影矩阵和均值
      std::pair<double, double> intervalLength;
      int dim = 2;
      //计算Q的m个T表
      for (int i = 0; i < m; ++i) {

          num[i] = end[i]-start[i]+1;
          vector<vector<double>> Q_sort1(num[i]);
          vector<vector<double>> Q_sort2(num[i]);
          //split
          vector<vector<double>> Q_original_part(Q.begin() + start[i], Q.begin() + end[i] + 1);
          //sort split
          splitQ_sort(Q_original_part,Q_sort1,Q_sort2,split_indices_sort[i],r);
          // debug------half
          // split_first_and_second(Q_original_part,Q_pca_part);
          // mu_sigma(Q_original_part,Q_original_part[0].size(),num[i],mu_Q_pca[i],sigma_Q_pca[i]);

          //mu and sigma
          // mu_sigma(Q_pca_part,dim,num[i],mu_Q_pca[i],sigma_Q_pca[i]);
          // split_two_dim(Q_pca_part,mu_Q_pca[i],sigma_Q_pca[i]);
          mu_sigma(Q_sort1,mu_Q_sort1[i],sigma_Q_sort1[i]);
          mu_sigma(Q_sort2,mu_Q_sort2[i],sigma_Q_sort2[i]);
          //U and L
          lower_upper_lemire_pca(num[i],mu_Q_sort1[i],mu_Q_sort2[i],U_mu_Q,L_mu_Q,U_sigma_Q,L_sigma_Q);
          //T
          pair<vector<double>,vector<double>> Q_part = {mu_Q_sort1[i],mu_Q_sort2[i]};
          intervalLength = initializeT(T[i], Q_part, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i],U_sigma_Q[i],lambda, num[i]);
          intervalLength_mu[i]=intervalLength.first;
          intervalLength_sigma[i]=intervalLength.second;
      }


    //计算四个dNN预选值
    for (int i = 0; i < m; ++i) {
        minDbl = MAX_VAL;
        for (int t = 0; t < num[i]; ++t) {
            Qt = {mu_Q_sort1[i][t], mu_Q_sort2[i][t]};
            corner = {L_mu_Q[i], L_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDbl = MIN(minDbl, dist);
        }
        result[0][i] = minDbl;
    }
    for (int i = 0; i < m; ++i) {
        minDbr = MAX_VAL;
        for (int t = 0; t < num[i]; ++t) {
            Qt = {mu_Q_sort1[i][t], mu_Q_sort2[i][t]};
            corner = {U_mu_Q[i], L_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDbr = MIN(minDbr, dist);
        }
        result[1][i] = minDbr;
    }
    for (int i = 0; i < m; ++i) {
        minDtl = MAX_VAL;
        for (int t = 0; t < num[i]; ++t) {
            Qt = {mu_Q_sort1[i][t], mu_Q_sort2[i][t]};
            corner = {L_mu_Q[i], U_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDtl = MIN(minDtl, dist);
        }
        result[2][i] = minDtl;
    }
    for (int i = 0; i < m; ++i) {
        minDtr = MAX_VAL;
        for (int t = 0; t < num[i]; ++t) {
            Qt = {mu_Q_sort1[i][t], mu_Q_sort2[i][t]};
            corner = {U_mu_Q[i], U_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDtr = MIN(minDtr, dist);
        }
        result[3][i] = minDtr;
    }



    for (int i = 0; i < result.size(); ++i) {
        for (int j = 0; j < result[i].size(); ++j) {
            result[i][j] = sqrt(result[i][j]);
        }
    }

    return result;
}

vector<vector<vector<double>>> Preprocessing_MC3E_PCA_4dim(vector<vector<double>>& Q,
    vector<double>& U_mu_Q1,vector<double>& L_mu_Q1,
    vector<double>& U_sigma_Q1,vector<double>& L_sigma_Q1,
    vector<double>& U_mu_Q2,vector<double>& L_mu_Q2,
    vector<double>& U_sigma_Q2,vector<double>& L_sigma_Q2,
    vector<vector<vector<Cell>>>& T1,
    vector <double>& intervalLength_mu1,vector <double>& intervalLength_sigma1 ,
    vector<vector<vector<Cell>>>& T2,
    vector <double>& intervalLength_mu2,vector <double>& intervalLength_sigma2 ,vector <int>& num,int r,
    vector<Eigen::MatrixXd>& projection_matrix_vector,vector<Eigen::VectorXd>& mean_Q_vector) {

    int lambda =3;
    int m = Q.size();
    vector<vector<double>> result1(4, vector<double>(m));
      vector<vector<double>> result2(4, vector<double>(m));
      vector<vector<vector<double>>> result(2);
    double minDbl,minDbr,minDtl,minDtr;
    vector<int> start(m);
    vector<int> end(m);
    vector<double> Qt(2);
    vector<double> corner(2);

    vector<vector<double>> mu_Q_pca1(m);
    vector<vector<double>> mu_Q_pca2(m);
    vector<vector<double>> sigma_Q_pca1(m);
    vector<vector<double>> sigma_Q_pca2(m);
    vector<vector<double>> Q_pca_part;
    vector<vector<double>> Q_pca_part1;
    vector<vector<double>> Q_pca_part2;

    for (int i = 0; i < m; ++i) {
        start[i] = MAX(0,i-r);
        end[i] = MIN(m-1,i+r);
    }

    //把每一组的Q降维并记录下投影矩阵和均值
      std::pair<double, double> intervalLength;
      int dim = 4;
      //计算Q的m个T表
      for (int i = 0; i < m; ++i) {
          num[i] = end[i]-start[i]+1;
          //split
          vector<vector<double>> Q_original_part(Q.begin() + start[i], Q.begin() + end[i] + 1);
          //PCA
          projection_matrix_vector[i] = pca_reduce_with_projection(Q_original_part,dim,Q_pca_part,mean_Q_vector[i]);
          //取一半维度分别处理
          split_half(Q_pca_part,Q_pca_part1,Q_pca_part2);

          //前两维
          //mu and sigma
          // mu_sigma(Q_pca_part,dim,num[i],mu_Q_pca1[i],sigma_Q_pca1[i]);
          split_two_dim(Q_pca_part,mu_Q_pca1[i],sigma_Q_pca1[i]);
          //U and L
          lower_upper_lemire_pca(num[i],mu_Q_pca1[i],sigma_Q_pca1[i],U_mu_Q1,L_mu_Q1,U_sigma_Q1,L_sigma_Q1);
          //T
          pair<vector<double>,vector<double>> Q_part = {mu_Q_pca1[i],sigma_Q_pca1[i]};
          intervalLength = initializeT(T1[i], Q_part, L_mu_Q1[i], U_mu_Q1[i], L_sigma_Q1[i],U_sigma_Q1[i],lambda, num[i]);
          intervalLength_mu1[i]=intervalLength.first;
          intervalLength_sigma1[i]=intervalLength.second;

          //后两维
          //mu and sigma
          // mu_sigma(Q_pca_part,dim,num[i],mu_Q_pca2[i],sigma_Q_pca2[i]);
          split_two_dim(Q_pca_part,mu_Q_pca2[i],sigma_Q_pca2[i]);
          //U and L
          lower_upper_lemire_pca(num[i],mu_Q_pca2[i],sigma_Q_pca2[i],U_mu_Q2,L_mu_Q2,U_sigma_Q2,L_sigma_Q2);
          //T
          Q_part = {mu_Q_pca2[i],sigma_Q_pca2[i]};
          intervalLength = initializeT(T2[i], Q_part, L_mu_Q2[i], U_mu_Q2[i], L_sigma_Q2[i],U_sigma_Q2[i],lambda, num[i]);
          intervalLength_mu2[i]=intervalLength.first;
          intervalLength_sigma2[i]=intervalLength.second;
      }


    //计算四个dNN预选值(1)
      for (int i = 0; i < m; ++i) {
          minDbl = MAX_VAL;
          minDbr = MAX_VAL;
          minDtl = MAX_VAL;
          minDtr = MAX_VAL;

          for (int t = 0; t < num[i]; ++t) {
              Qt = {mu_Q_pca1[i][t], sigma_Q_pca1[i][t]};

              // 左下角
              std::vector<double> corner_bl = {L_mu_Q1[i], L_sigma_Q1[i]};
              double dist_bl = euclideanDistance_2(Qt, corner_bl, Qt.size());
              minDbl = MIN(minDbl, dist_bl);

              // 右下角
              std::vector<double> corner_br = {U_mu_Q1[i], L_sigma_Q1[i]};
              double dist_br = euclideanDistance_2(Qt, corner_br, Qt.size());
              minDbr = MIN(minDbr, dist_br);

              // 左上角
              std::vector<double> corner_tl = {L_mu_Q1[i], U_sigma_Q1[i]};
              double dist_tl = euclideanDistance_2(Qt, corner_tl, Qt.size());
              minDtl = MIN(minDtl, dist_tl);

              // 右上角
              std::vector<double> corner_tr = {U_mu_Q1[i], U_sigma_Q1[i]};
              double dist_tr = euclideanDistance_2(Qt, corner_tr, Qt.size());
              minDtr = MIN(minDtr, dist_tr);
          }

          result1[0][i] = minDbl;
          result1[1][i] = minDbr;
          result1[2][i] = minDtl;
          result1[3][i] = minDtr;
      }

      for (int i = 0; i < result1.size(); ++i) {
          for (int j = 0; j < result1[i].size(); ++j) {
              result1[i][j] = std::sqrt(result1[i][j]);
          }
      }

      //计算四个dNN预选值(2)
      for (int i = 0; i < m; ++i) {
          minDbl = MAX_VAL;
          minDbr = MAX_VAL;
          minDtl = MAX_VAL;
          minDtr = MAX_VAL;

          for (int t = 0; t < num[i]; ++t) {
              Qt = {mu_Q_pca2[i][t], sigma_Q_pca2[i][t]};

              // 左下角
              std::vector<double> corner_bl = {L_mu_Q2[i], L_sigma_Q2[i]};
              double dist_bl = euclideanDistance_2(Qt, corner_bl, Qt.size());
              minDbl = MIN(minDbl, dist_bl);

              // 右下角
              std::vector<double> corner_br = {U_mu_Q2[i], L_sigma_Q2[i]};
              double dist_br = euclideanDistance_2(Qt, corner_br, Qt.size());
              minDbr = MIN(minDbr, dist_br);

              // 左上角
              std::vector<double> corner_tl = {L_mu_Q2[i], U_sigma_Q2[i]};
              double dist_tl = euclideanDistance_2(Qt, corner_tl, Qt.size());
              minDtl = MIN(minDtl, dist_tl);

              // 右上角
              std::vector<double> corner_tr = {U_mu_Q2[i], U_sigma_Q2[i]};
              double dist_tr = euclideanDistance_2(Qt, corner_tr, Qt.size());
              minDtr = MIN(minDtr, dist_tr);
          }

          result2[0][i] = minDbl;
          result2[1][i] = minDbr;
          result2[2][i] = minDtl;
          result2[3][i] = minDtr;
      }

      for (int i = 0; i < result2.size(); ++i) {
          for (int j = 0; j < result2[i].size(); ++j) {
              result2[i][j] = std::sqrt(result2[i][j]);
          }
      }
    result[0]=result1;
      result[1]=result2;
    return result;
}

vector<vector<double>> Preprocessing_MC3E_DIVIDE(vector<vector<double>>& Q,
    vector<double>& mu_Q,vector<double>& sigma_Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,int r) {
    int lambda =3;
    int m = Q.size();
    lower_upper_lemire(mu_Q, m, r, L_mu_Q, U_mu_Q);
    lower_upper_lemire(sigma_Q, m, r, L_sigma_Q, U_sigma_Q);

    vector<vector<double>> result(4, vector<double>(m));


    double minDbl,minDbr,minDtl,minDtr;
    vector<int> start(m);
    vector<int> end(m);
    vector<double> Qt(2);
    vector<double> corner(2);
    double dist;

    for (int i = 0; i < m; ++i) {
        start[i] = MAX(0,i-r);
        end[i] = MIN(m-1,i+r);
    }
    //计算四个dNN预选值
    for (int i = 0; i < m; ++i) {
        minDbl = MAX_VAL;
        for (int t = start[i]; t <= end[i]; ++t) {
            Qt = {mu_Q[t], sigma_Q[t]};
            corner = {L_mu_Q[i], L_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDbl = MIN(minDbl, dist);
        }
        result[0][i] = minDbl;
    }
    for (int i = 0; i < m; ++i) {
        minDbr = MAX_VAL;
        for (int t = start[i]; t <= end[i]; ++t) {
            Qt = {mu_Q[t], sigma_Q[t]};
            corner = {U_mu_Q[i], L_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDbr = MIN(minDbr, dist);
        }
        result[1][i] = minDbr;
    }
    for (int i = 0; i < m; ++i) {
        minDtl = MAX_VAL;
        for (int t = start[i]; t <= end[i]; ++t) {
            Qt = {mu_Q[t], sigma_Q[t]};
            corner = {L_mu_Q[i], U_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDtl = MIN(minDtl, dist);
        }
        result[2][i] = minDtl;
    }
    for (int i = 0; i < m; ++i) {
        minDtr = MAX_VAL;
        for (int t = start[i]; t <= end[i]; ++t) {
            Qt = {mu_Q[t], sigma_Q[t]};
            corner = {U_mu_Q[i], U_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            minDtr = MIN(minDtr, dist);
        }
        result[3][i] = minDtr;
    }



    std::pair<double, double> intervalLength;
    //计算Q的m个T表
    for (int i = 0; i < m; ++i) {
        num[i] = end[i]-start[i]+1;
        vector<double> mu_Q_part(mu_Q.begin() + start[i], mu_Q.begin() + end[i] + 1);
        vector<double> sigma_Q_part(sigma_Q.begin() + start[i], sigma_Q.begin() + end[i] + 1);
        pair<vector<double>,vector<double>> Q_part = {mu_Q_part,sigma_Q_part};
        intervalLength = initializeT_divide(T[i], Q_part, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i],U_sigma_Q[i],lambda, num[i]);
        intervalLength_mu[i]=intervalLength.first;
        intervalLength_sigma[i]=intervalLength.second;
    }

    for (int i = 0; i < result.size(); ++i) {
        for (int j = 0; j < result[i].size(); ++j) {
            result[i][j] = sqrt(result[i][j]);
        }
    }

    return result;
}

double LB_MC3E(vector<vector<double>>& Q, vector<vector<double>>& A,
    vector<double>& mu_Q,vector<double>& sigma_Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<double>& Dbl, vector<double>& Dbr,
    vector<double>& Dtl, vector<double>& Dtr,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,
    int K)
{

    // std::ofstream outputFile_data;
    // outputFile_data.open("..//data.txt",std::ios::app);

    int m = A.size();
    int lambda  = 3;
    vector<double> mu_A(m);
    vector<double> sigma_A(m);
      if (K==2) {
          split_two_dim(A,mu_A,sigma_A);
      }else {
      mu_sigma(A,mu_A,sigma_A);
      }
    double lb = 0.0;

    vector<double> Ai(2);
    vector<double> V(2);
    double dNN;
    vector<double> V_plus(2);
    vector<double> V_minus(2);
    double d_plus;
    double d_minus;
    bool MiddleLeftOrRight=false;
    bool TopOrBottomMiddle=false;


    for (int i = 0; i < m; ++i) {
        Ai = {mu_A[i],sigma_A[i]};

        //concer
        if (Ai[0] < L_mu_Q[i] && Ai[1] < L_sigma_Q[i]) {
            V = {L_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbl[i];
            lb += CornerDistE(Ai, V, dNN);
        } else if (Ai[0] > U_mu_Q[i] && Ai[1] < L_sigma_Q[i]) {
            V = {U_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbr[i];
            lb += CornerDistE(Ai, V, dNN);
        } else if (Ai[0] < L_mu_Q[i] && Ai[1] > U_sigma_Q[i]) {
            V = {L_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtl[i];
            lb += CornerDistE(Ai, V, dNN);
        } else if (Ai[0] > U_mu_Q[i] && Ai[1] > U_sigma_Q[i]){
            V = {U_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtr[i];
            lb += CornerDistE(Ai, V, dNN);
        }
        //edge
        else if (Ai[0] <= L_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i]) {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {L_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbl[i];
            d_minus = Dtl[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle=false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        } else if (Ai[0] >= U_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i]) {
            V_plus = {U_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbr[i];
            d_minus = Dtr[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle=false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        } else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] <= L_sigma_Q[i]) {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], L_sigma_Q[i]};
            d_plus = Dbl[i];
            d_minus = Dbr[i];
            MiddleLeftOrRight=false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        } else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] >= U_sigma_Q[i]){
            V_plus = {L_mu_Q[i], U_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dtl[i];
            d_minus = Dtr[i];
            MiddleLeftOrRight=false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        }
        //center
         else {
            lb += CenterDist2D(mu_A[i],sigma_A[i],L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
        }
    }
      if (K==2) {
          return sqrt(lb);
      }
    return sqrt(K*lb);
}



double LB_MC3E_sort(
    vector<double>& mu_A,vector<double>& sigma_A,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<double>& Dbl, vector<double>& Dbr,
    vector<double>& Dtl, vector<double>& Dtr,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,
    int m,int K)
{

    int lambda  = 3;
    double lb = 0.0;

    vector<double> Ai(2);
    vector<double> V(2);
    double dNN;
    vector<double> V_plus(2);
    vector<double> V_minus(2);
    double d_plus;
    double d_minus;
    bool MiddleLeftOrRight=false;
    bool TopOrBottomMiddle=false;


    for (int i = 0; i < m; ++i) {
        Ai = {mu_A[i],sigma_A[i]};

        //concer
        if (Ai[0] < L_mu_Q[i] && Ai[1] < L_sigma_Q[i]) {
            V = {L_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbl[i];
            lb += CornerDistE(Ai, V, dNN);
        } else if (Ai[0] > U_mu_Q[i] && Ai[1] < L_sigma_Q[i]) {
            V = {U_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbr[i];
            lb += CornerDistE(Ai, V, dNN);
        } else if (Ai[0] < L_mu_Q[i] && Ai[1] > U_sigma_Q[i]) {
            V = {L_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtl[i];
            lb += CornerDistE(Ai, V, dNN);
        } else if (Ai[0] > U_mu_Q[i] && Ai[1] > U_sigma_Q[i]){
            V = {U_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtr[i];
            lb += CornerDistE(Ai, V, dNN);
        }
        //edge
        else if (Ai[0] <= L_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i]) {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {L_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbl[i];
            d_minus = Dtl[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle=false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        } else if (Ai[0] >= U_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i]) {
            V_plus = {U_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbr[i];
            d_minus = Dtr[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle=false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        } else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] <= L_sigma_Q[i]) {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], L_sigma_Q[i]};
            d_plus = Dbl[i];
            d_minus = Dbr[i];
            MiddleLeftOrRight=false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        } else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] >= U_sigma_Q[i]){
            V_plus = {L_mu_Q[i], U_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dtl[i];
            d_minus = Dtr[i];
            MiddleLeftOrRight=false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        }
        //center
         else {
            lb += CenterDist2D(mu_A[i],sigma_A[i],L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
        }
    }

    return sqrt(K*lb);
}



double LB_MC3E_2dim(vector<vector<double>>& Q, vector<vector<double>>& A,
    vector<double>& mu_Q,vector<double>& sigma_Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<double>& Dbl, vector<double>& Dbr,
    vector<double>& Dtl, vector<double>& Dtr,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,
    int K)
{

    // std::ofstream outputFile_data;
    // outputFile_data.open("..//data.txt",std::ios::app);

    int m = A.size();
    int lambda  = 3;
    vector<double> mu_A(m);
    vector<double> sigma_A(m);
    split_two_dim(A,mu_A,sigma_A);
    double lb = 0.0;

    vector<double> Ai(2);
    vector<double> V(2);
    double dNN;
    vector<double> V_plus(2);
    vector<double> V_minus(2);
    double d_plus;
    double d_minus;
    bool MiddleLeftOrRight=false;
    bool TopOrBottomMiddle=false;


    for (int i = 0; i < m; ++i) {
        Ai = {mu_A[i],sigma_A[i]};

        //concer
        if (Ai[0] < L_mu_Q[i] && Ai[1] < L_sigma_Q[i]) {
            V = {L_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbl[i];
            lb += CornerDistE(Ai, V, dNN);
        } else if (Ai[0] > U_mu_Q[i] && Ai[1] < L_sigma_Q[i]) {
            V = {U_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbr[i];
            lb += CornerDistE(Ai, V, dNN);
        } else if (Ai[0] < L_mu_Q[i] && Ai[1] > U_sigma_Q[i]) {
            V = {L_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtl[i];
            lb += CornerDistE(Ai, V, dNN);
        } else if (Ai[0] > U_mu_Q[i] && Ai[1] > U_sigma_Q[i]){
            V = {U_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtr[i];
            lb += CornerDistE(Ai, V, dNN);
        }
        //edge
        else if (Ai[0] <= L_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i]) {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {L_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbl[i];
            d_minus = Dtl[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle=false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        } else if (Ai[0] >= U_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i]) {
            V_plus = {U_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbr[i];
            d_minus = Dtr[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle=false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        } else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] <= L_sigma_Q[i]) {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], L_sigma_Q[i]};
            d_plus = Dbl[i];
            d_minus = Dbr[i];
            MiddleLeftOrRight=false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        } else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] >= U_sigma_Q[i]){
            V_plus = {L_mu_Q[i], U_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dtl[i];
            d_minus = Dtr[i];
            MiddleLeftOrRight=false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        }
        //center
         else {
            lb += CenterDist2D(mu_A[i],sigma_A[i],L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
        }
    }
    return sqrt(lb);
}


double LB_MC3E_2dim_pca(vector<vector<double>>& A,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<double>& Dbl, vector<double>& Dbr,
    vector<double>& Dtl, vector<double>& Dtr,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,
    vector<Eigen::MatrixXd>& projection_matrix_vector,vector<Eigen::VectorXd>& mean_Q_vector)
{

    int m = A.size();
      int K = A[0].size();
    int lambda  = 3;
    double lb = 0.0;
    vector<double> Ai(2);
    vector<double> V(2);
    double dNN;
    vector<double> V_plus(2);
    vector<double> V_minus(2);
    double d_plus;
    double d_minus;
    bool MiddleLeftOrRight=false;
    bool TopOrBottomMiddle=false;
    vector<double> A_pca;
      //debug
      // vector<double> A_mu(m);
      // vector<double> A_sigma(m);
      // mu_sigma(A,K,m,A_mu,A_sigma);

    for (int i = 0; i < m; ++i) {
        //PCA
        reduce_single_A_with_Q_projection(A[i],projection_matrix_vector[i],mean_Q_vector[i],A_pca);
        //debug----half
        // A_pca={A[i][0],A[i][1]};
        //mu and sigma
        Ai = {A_pca[0],A_pca[1]};
        // double mu = (A_pca[0]+A_pca[1])/2;
        // double sigma = std::sqrt((std::pow(A_pca[0] - mu, 2) + std::pow(A_pca[1] - mu, 2)) / 2);
        // Ai = {mu,sigma};
        // Ai={A_mu[i],A_sigma[i]};
        //corner
        if (Ai[0] < L_mu_Q[i] && Ai[1] < L_sigma_Q[i]) {
            V = {L_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbl[i];
            lb += CornerDistE(Ai, V, dNN);
        } else if (Ai[0] > U_mu_Q[i] && Ai[1] < L_sigma_Q[i]) {
            V = {U_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbr[i];
            lb += CornerDistE(Ai, V, dNN);
        } else if (Ai[0] < L_mu_Q[i] && Ai[1] > U_sigma_Q[i]) {
            V = {L_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtl[i];
            lb += CornerDistE(Ai, V, dNN);
        } else if (Ai[0] > U_mu_Q[i] && Ai[1] > U_sigma_Q[i]){
            V = {U_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtr[i];
            lb += CornerDistE(Ai, V, dNN);
        }
        //edge
        else if (Ai[0] <= L_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i]) {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {L_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbl[i];
            d_minus = Dtl[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle=false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        } else if (Ai[0] >= U_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i]) {
            V_plus = {U_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbr[i];
            d_minus = Dtr[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle=false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        } else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] <= L_sigma_Q[i]) {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], L_sigma_Q[i]};
            d_plus = Dbl[i];
            d_minus = Dbr[i];
            MiddleLeftOrRight=false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        } else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] >= U_sigma_Q[i]){
            V_plus = {L_mu_Q[i], U_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dtl[i];
            d_minus = Dtr[i];
            MiddleLeftOrRight=false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        }
         else {//center
            lb += CenterDist2D(Ai[0],Ai[1],L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
        }
    }
    return sqrt(lb);
}


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
    vector<Eigen::MatrixXd>& projection_matrix_vector,vector<Eigen::VectorXd>& mean_Q_vector)
{

    int m = A.size();
    int lambda  = 3;
    double lb1 = 0.0;
    double lb2 = 0.0;
    vector<double> Ai1(2);
      vector<double> Ai2(2);
    vector<double> V(2);
    double dNN;
    vector<double> V_plus(2);
    vector<double> V_minus(2);
    double d_plus;
    double d_minus;
    bool MiddleLeftOrRight=false;
    bool TopOrBottomMiddle=false;
    vector<double> A_pca;

    for (int i = 0; i < m; ++i) {
        //PCA
        reduce_single_A_with_Q_projection(A[i],projection_matrix_vector[i],mean_Q_vector[i],A_pca);
        //mu and sigma
        Ai1 = {A_pca[0],A_pca[1]};
        Ai2 = {A_pca[2],A_pca[3]};
        // double mu = (A_pca[0]+A_pca[1])/2;
        // double sigma = std::sqrt((std::pow(A_pca[0] - mu, 2) + std::pow(A_pca[1] - mu, 2)) / 2);
        // Ai = {mu,sigma};
        //concer
        if (Ai1[0] < L_mu_Q1[i] && Ai1[1] < L_sigma_Q1[i]) {
            V = {L_mu_Q1[i], L_sigma_Q1[i]};
            dNN = Dbl1[i];
            lb1 += CornerDistE(Ai1, V, dNN);
        } else if (Ai1[0] > U_mu_Q1[i] && Ai1[1] < L_sigma_Q1[i]) {
            V = {U_mu_Q1[i], L_sigma_Q1[i]};
            dNN = Dbr1[i];
            lb1 += CornerDistE(Ai1, V, dNN);
        } else if (Ai1[0] < L_mu_Q1[i] && Ai1[1] > U_sigma_Q1[i]) {
            V = {L_mu_Q1[i], U_sigma_Q1[i]};
            dNN = Dtl1[i];
            lb1 += CornerDistE(Ai1, V, dNN);
        } else if (Ai1[0] > U_mu_Q1[i] && Ai1[1] > U_sigma_Q1[i]){
            V = {U_mu_Q1[i], U_sigma_Q1[i]};
            dNN = Dtr1[i];
            lb1 += CornerDistE(Ai1, V, dNN);
        }
        //edge
        else if (Ai1[0] <= L_mu_Q1[i] && L_sigma_Q1[i] <= Ai1[1] && Ai1[1] <= U_sigma_Q1[i]) {
            V_plus = {L_mu_Q1[i], L_sigma_Q1[i]};
            V_minus = {L_mu_Q1[i], U_sigma_Q1[i]};
            d_plus = Dbl1[i];
            d_minus = Dtl1[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle=false;
            double delta_plus = EdgeDistE(Ai1, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q1[i],U_mu_Q1[i],L_sigma_Q1[i],U_sigma_Q1[i],T1[i],intervalLength_mu1[i],intervalLength_sigma1[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai1, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q1[i],U_mu_Q1[i],L_sigma_Q1[i],U_sigma_Q1[i],T1[i],intervalLength_mu1[i],intervalLength_sigma1[i],lambda,num[i]);
            lb1 += MIN(delta_plus, delta_minus);
        } else if (Ai1[0] >= U_mu_Q1[i] && L_sigma_Q1[i] <= Ai1[1] && Ai1[1] <= U_sigma_Q1[i]) {
            V_plus = {U_mu_Q1[i], L_sigma_Q1[i]};
            V_minus = {U_mu_Q1[i], U_sigma_Q1[i]};
            d_plus = Dbr1[i];
            d_minus = Dtr1[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle=false;
            double delta_plus = EdgeDistE(Ai1, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q1[i],U_mu_Q1[i],L_sigma_Q1[i],U_sigma_Q1[i],T1[i],intervalLength_mu1[i],intervalLength_sigma1[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai1, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q1[i],U_mu_Q1[i],L_sigma_Q1[i],U_sigma_Q1[i],T1[i],intervalLength_mu1[i],intervalLength_sigma1[i],lambda,num[i]);
            lb1 += MIN(delta_plus, delta_minus);
        } else if (L_mu_Q1[i] <= Ai1[0] && Ai1[0] <= U_mu_Q1[i] && Ai1[1] <= L_sigma_Q1[i]) {
            V_plus = {L_mu_Q1[i], L_sigma_Q1[i]};
            V_minus = {U_mu_Q1[i], L_sigma_Q1[i]};
            d_plus = Dbl1[i];
            d_minus = Dbr1[i];
            MiddleLeftOrRight=false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai1, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q1[i],U_mu_Q1[i],L_sigma_Q1[i],U_sigma_Q1[i],T1[i],intervalLength_mu1[i],intervalLength_sigma1[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai1, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q1[i],U_mu_Q1[i],L_sigma_Q1[i],U_sigma_Q1[i],T1[i],intervalLength_mu1[i],intervalLength_sigma1[i],lambda,num[i]);
            lb1 += MIN(delta_plus, delta_minus);
        } else if (L_mu_Q1[i] <= Ai1[0] && Ai1[0] <= U_mu_Q1[i] && Ai1[1] >= U_sigma_Q1[i]){
            V_plus = {L_mu_Q1[i], U_sigma_Q1[i]};
            V_minus = {U_mu_Q1[i], U_sigma_Q1[i]};
            d_plus = Dtl1[i];
            d_minus = Dtr1[i];
            MiddleLeftOrRight=false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai1, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q1[i],U_mu_Q1[i],L_sigma_Q1[i],U_sigma_Q1[i],T1[i],intervalLength_mu1[i],intervalLength_sigma1[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai1, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q1[i],U_mu_Q1[i],L_sigma_Q1[i],U_sigma_Q1[i],T1[i],intervalLength_mu1[i],intervalLength_sigma1[i],lambda,num[i]);
            lb1 += MIN(delta_plus, delta_minus);
        }
         else {//center
            lb1 += CenterDist2D(Ai1[0],Ai1[1],L_mu_Q1[i],U_mu_Q1[i],L_sigma_Q1[i],U_sigma_Q1[i],T1[i],intervalLength_mu1[i],intervalLength_sigma1[i],lambda,num[i]);
        }


        if (Ai2[0] < L_mu_Q2[i] && Ai2[1] < L_sigma_Q2[i]) {
            V = {L_mu_Q2[i], L_sigma_Q2[i]};
            dNN = Dbl2[i];
            lb2 += CornerDistE(Ai2, V, dNN);
        } else if (Ai2[0] > U_mu_Q2[i] && Ai2[1] < L_sigma_Q2[i]) {
            V = {U_mu_Q2[i], L_sigma_Q2[i]};
            dNN = Dbr2[i];
            lb2 += CornerDistE(Ai2, V, dNN);
        } else if (Ai2[0] < L_mu_Q2[i] && Ai2[1] > U_sigma_Q2[i]) {
            V = {L_mu_Q2[i], U_sigma_Q2[i]};
            dNN = Dtl2[i];
            lb2 += CornerDistE(Ai2, V, dNN);
        } else if (Ai2[0] > U_mu_Q2[i] && Ai2[1] > U_sigma_Q2[i]) {
            V = {U_mu_Q2[i], U_sigma_Q2[i]};
            dNN = Dtr2[i];
            lb2 += CornerDistE(Ai2, V, dNN);
        }
        // edge
        else if (Ai2[0] <= L_mu_Q2[i] && L_sigma_Q2[i] <= Ai2[1] && Ai2[1] <= U_sigma_Q2[i]) {
            V_plus = {L_mu_Q2[i], L_sigma_Q2[i]};
            V_minus = {L_mu_Q2[i], U_sigma_Q2[i]};
            d_plus = Dbl2[i];
            d_minus = Dtl2[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle = false;
            double delta_plus = EdgeDistE(Ai2, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q2[i], U_mu_Q2[i], L_sigma_Q2[i], U_sigma_Q2[i], T2[i], intervalLength_mu2[i], intervalLength_sigma2[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai2, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q2[i], U_mu_Q2[i], L_sigma_Q2[i], U_sigma_Q2[i], T2[i], intervalLength_mu2[i], intervalLength_sigma2[i], lambda, num[i]);
            lb2 += MIN(delta_plus, delta_minus);
        } else if (Ai2[0] >= U_mu_Q2[i] && L_sigma_Q2[i] <= Ai2[1] && Ai2[1] <= U_sigma_Q2[i]) {
            V_plus = {U_mu_Q2[i], L_sigma_Q2[i]};
            V_minus = {U_mu_Q2[i], U_sigma_Q2[i]};
            d_plus = Dbr2[i];
            d_minus = Dtr2[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle = false;
            double delta_plus = EdgeDistE(Ai2, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q2[i], U_mu_Q2[i], L_sigma_Q2[i], U_sigma_Q2[i], T2[i], intervalLength_mu2[i], intervalLength_sigma2[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai2, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q2[i], U_mu_Q2[i], L_sigma_Q2[i], U_sigma_Q2[i], T2[i], intervalLength_mu2[i], intervalLength_sigma2[i], lambda, num[i]);
            lb2 += MIN(delta_plus, delta_minus);
        } else if (L_mu_Q2[i] <= Ai2[0] && Ai2[0] <= U_mu_Q2[i] && Ai2[1] <= L_sigma_Q2[i]) {
            V_plus = {L_mu_Q2[i], L_sigma_Q2[i]};
            V_minus = {U_mu_Q2[i], L_sigma_Q2[i]};
            d_plus = Dbl2[i];
            d_minus = Dbr2[i];
            MiddleLeftOrRight = false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai2, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q2[i], U_mu_Q2[i], L_sigma_Q2[i], U_sigma_Q2[i], T2[i], intervalLength_mu2[i], intervalLength_sigma2[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai2, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q2[i], U_mu_Q2[i], L_sigma_Q2[i], U_sigma_Q2[i], T2[i], intervalLength_mu2[i], intervalLength_sigma2[i], lambda, num[i]);
            lb2 += MIN(delta_plus, delta_minus);
        } else if (L_mu_Q2[i] <= Ai2[0] && Ai2[0] <= U_mu_Q2[i] && Ai2[1] >= U_sigma_Q2[i]) {
            V_plus = {L_mu_Q2[i], U_sigma_Q2[i]};
            V_minus = {U_mu_Q2[i], U_sigma_Q2[i]};
            d_plus = Dtl2[i];
            d_minus = Dtr2[i];
            MiddleLeftOrRight = false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai2, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q2[i], U_mu_Q2[i], L_sigma_Q2[i], U_sigma_Q2[i], T2[i], intervalLength_mu2[i], intervalLength_sigma2[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai2, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q2[i], U_mu_Q2[i], L_sigma_Q2[i], U_sigma_Q2[i], T2[i], intervalLength_mu2[i], intervalLength_sigma2[i], lambda, num[i]);
            lb2 += MIN(delta_plus, delta_minus);
        }
        else { // center
            lb2 += CenterDist2D(Ai2[0], Ai2[1], L_mu_Q2[i], U_mu_Q2[i], L_sigma_Q2[i], U_sigma_Q2[i], T2[i], intervalLength_mu2[i], intervalLength_sigma2[i], lambda, num[i]);
        }
    }

return sqrt(lb1+lb2);
}

double LB_MC3E_sort_strength(vector<vector<double>>& A,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<double>& Dbl, vector<double>& Dbr,
    vector<double>& Dtl, vector<double>& Dtr,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,
    vector<vector<int>> split_indices_sort)
{
    int m = A.size();
    int K = A[0].size();
    int lambda  = 3;
    double lb = 0.0;
    vector<double> Ai(2);
    vector<double> V(2);
    double dNN;
    vector<double> V_plus(2);
    vector<double> V_minus(2);
    double d_plus;
    double d_minus;
    bool MiddleLeftOrRight=false;
    bool TopOrBottomMiddle=false;
    vector<double> A_sort1;
    vector<double> A_sort2;
  double mu_A1,sigma_A1;
  double mu_A2,sigma_A2;
      //debug
      // vector<double> A_mu(m);
      // vector<double> A_sigma(m);
      // mu_sigma(A,K,m,A_mu,A_sigma);

    for (int i = 0; i < m; ++i) {
        //sort split
        splitA_sort_single(A[i],split_indices_sort[i],A_sort1,A_sort2);
        //debug----half
        // A_pca={A[i][0],A[i][1]};
        //mu and sigma
        mu_sigma_single(A_sort1,mu_A1,sigma_A1);
        mu_sigma_single(A_sort2,mu_A2,sigma_A2);
        Ai = {mu_A1,mu_A2};
        // double mu = (A_pca[0]+A_pca[1])/2;
        // double sigma = std::sqrt((std::pow(A_pca[0] - mu, 2) + std::pow(A_pca[1] - mu, 2)) / 2);
        // Ai = {mu,sigma};
        // Ai={A_mu[i],A_sigma[i]};
        //corner
        if (Ai[0] < L_mu_Q[i] && Ai[1] < L_sigma_Q[i]) {
            V = {L_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbl[i];
            lb += CornerDistE(Ai, V, dNN);
        } else if (Ai[0] > U_mu_Q[i] && Ai[1] < L_sigma_Q[i]) {
            V = {U_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbr[i];
            lb += CornerDistE(Ai, V, dNN);
        } else if (Ai[0] < L_mu_Q[i] && Ai[1] > U_sigma_Q[i]) {
            V = {L_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtl[i];
            lb += CornerDistE(Ai, V, dNN);
        } else if (Ai[0] > U_mu_Q[i] && Ai[1] > U_sigma_Q[i]){
            V = {U_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtr[i];
            lb += CornerDistE(Ai, V, dNN);
        }
        //edge
        else if (Ai[0] <= L_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i]) {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {L_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbl[i];
            d_minus = Dtl[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle=false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        } else if (Ai[0] >= U_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i]) {
            V_plus = {U_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbr[i];
            d_minus = Dtr[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle=false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        } else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] <= L_sigma_Q[i]) {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], L_sigma_Q[i]};
            d_plus = Dbl[i];
            d_minus = Dbr[i];
            MiddleLeftOrRight=false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        } else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] >= U_sigma_Q[i]){
            V_plus = {L_mu_Q[i], U_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dtl[i];
            d_minus = Dtr[i];
            MiddleLeftOrRight=false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            lb += MIN(delta_plus, delta_minus);
        }
         else {//center
            lb += CenterDist2D(Ai[0],Ai[1],L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
        }
    }
    return sqrt(K*lb);
}
