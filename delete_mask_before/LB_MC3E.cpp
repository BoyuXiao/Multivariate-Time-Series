
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <utility>
#include "funcs.h"
using namespace std;

// 优化版本
  pair<double,double> initializeT(std::vector<vector<Cell>>& T, std::pair<vector<double>,vector<double>>& Q, double L_mu_Q, double U_mu_Q, double L_sigma_Q, double U_sigma_Q, int lambda, int n,int start,int end) {
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
      for (int i = start; i <= end; ++i) {
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
          T[sigmaIndex][muIndex].dataList.emplace_back(mu,sigma,i);
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
                  T[i][j].left_bottom_t = T[row][col].left_bottom_t;
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
                      T[i][j].left_bottom_t = T[row][col].left_bottom_t;
                  }
              // }
              row = MAX(i-1,0);
              col = MAX(j-1,0);
              //遍历一个格子来给点赋值
              for (int k=0;k<T[row][col].dataList.size();k++) {
                  if (i<1||j<1) break;
                  mu = get<0>(T[row][col].dataList[k]);
                  sigma = get<1>(T[row][col].dataList[k]);
                  dist_X = sqrt(dist(T[i][j].Lleft_bottom.first,mu)+dist(T[i][j].Lleft_bottom.second,sigma));
                  if (dist_X < T[i][j].left_bottom_nearest_dist) {
                      T[i][j].left_bottom_nearest_dist = dist_X;
                      T[i][j].left_bottom_nearest_mu = mu;
                      T[i][j].left_bottom_nearest_sigma = sigma;
                      T[i][j].left_bottom_t = get<2>(T[row][col].dataList[k]);
                  }
              }
              //遍历一个格子来给left edge赋值
              row = i;
              col = MAX(0,j-1);
              if (T[row][col].dataList.size()!=0) {
                  for (int k=0;k<T[row][col].dataList.size();k++) {
                      if (j<1)break;
                      mu = get<0>(T[row][col].dataList[k]);
                      sigma = get<1>(T[row][col].dataList[k]);
                      dist_X = abs(T[i][j].Lleft_bottom.first - mu);
                      if (dist_X < T[i][j].minDistToLleft) {
                          T[i][j].minDistToLleft = dist_X;
                          T[i][j].left_nearest_mu = mu;
                          T[i][j].left_nearest_sigma = sigma;
                          T[i][j].left_t = get<2>(T[row][col].dataList[k]);
                      }
                  }
              }else {
                  // if (T[row][col].minDistToLleft!=MAX_VAL) {
                      mu = T[row][col].left_nearest_mu;
                      sigma = T[row][col].left_nearest_sigma;
                      dist_X = abs(T[i][j].Lleft_bottom.first - mu);
                      T[i][j].minDistToLleft = dist_X;
                      T[i][j].left_nearest_mu = mu;
                      T[i][j].left_nearest_sigma = sigma;
                      T[i][j].left_t = T[row][col].left_t;
                  // }
              }
              row = MAX(i-1,0);
              col = j;
              //遍历一个格子来给bottom edge赋值
              if (T[row][col].dataList.size()!=0) {
                  for (int k=0;k<T[row][col].dataList.size();k++) {
                      if (i<1)break;
                      mu = get<0>(T[row][col].dataList[k]);
                      sigma = get<1>(T[row][col].dataList[k]);
                      dist_X = abs(T[i][j].Lright_bottom.second - sigma);
                      if (dist_X < T[i][j].minDistToLbottom) {
                          T[i][j].minDistToLbottom = dist_X;
                          T[i][j].bottom_nearest_mu = mu;
                          T[i][j].bottom_nearest_sigma = sigma;
                          T[i][j].bottom_t = get<2>(T[row][col].dataList[k]);
                      }
                  }
              }else {
                  // if (T[row][col].minDistToLbottom!=MAX_VAL) {
                      mu = T[row][col].bottom_nearest_mu;
                      sigma = T[row][col].bottom_nearest_sigma;
                      dist_X = abs(T[i][j].Lright_bottom.second - sigma);
                      T[i][j].minDistToLbottom = dist_X;
                      T[i][j].bottom_nearest_mu = mu;
                      T[i][j].bottom_nearest_sigma = sigma;
                      T[i][j].bottom_t = T[row][col].bottom_t;
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
                  T[i][j].right_top_t = T[row][col].right_top_t;
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
                      T[i][j].right_top_t = T[row][col].right_top_t;
                  }
              // }
              row = MIN(i+1,numIntervals-1);
              col = MIN(j+1,numIntervals-1);
              //遍历一个格子来给点赋值
              for (int k=0;k<T[row][col].dataList.size();k++) {
                  if (i+2>numIntervals||j+2>numIntervals)break;
                  mu = get<0>(T[row][col].dataList[k]);
                  sigma = get<1>(T[row][col].dataList[k]);
                  dist_X = sqrt(dist(T[i][j].Lright_top.first,mu)+dist(T[i][j].Lright_top.second,sigma));
                  if (dist_X < T[i][j].right_top_nearest_dist) {
                      T[i][j].right_top_nearest_dist = dist_X;
                      T[i][j].right_top_nearest_mu = mu;
                      T[i][j].right_top_nearest_sigma = sigma;
                      T[i][j].right_top_t = get<2>(T[row][col].dataList[k]);
                  }
              }
              //遍历一个格子来给right edge赋值
              row = i;
              col = MIN(j+1,numIntervals - 1);
              if (T[row][col].dataList.size()!=0) {
                  for (int k=0;k<T[row][col].dataList.size();k++) {
                      if (j+2>numIntervals)break;
                      mu = get<0>(T[row][col].dataList[k]);
                      sigma = get<1>(T[row][col].dataList[k]);
                      dist_X = abs(mu - T[i][j].Lright_bottom.first);
                      if (dist_X < T[i][j].minDistToLright) {
                          T[i][j].minDistToLright = dist_X;
                          T[i][j].right_nearest_mu = mu;
                          T[i][j].right_nearest_sigma = sigma;
                          T[i][j].right_t = get<2>(T[row][col].dataList[k]);
                      }
                  }
              }else {
                  // if (T[row][col].minDistToLright!=MAX_VAL) {
                      mu = T[row][col].right_nearest_mu;
                      sigma = T[row][col].right_nearest_sigma;
                      dist_X = abs(mu - T[i][j].Lright_bottom.first);
                      T[i][j].minDistToLright = dist_X;
                      T[i][j].right_nearest_mu = mu;
                      T[i][j].right_nearest_sigma = sigma;
                      T[i][j].right_t = T[row][col].right_t;
                  // }
              }
              row = MIN(i+1,numIntervals-1);
              col = j;
              //遍历一个格子来给top edge赋值
              if (T[row][col].dataList.size()!=0) {
                  for (int k=0;k<T[row][col].dataList.size();k++) {
                      if (i+2>numIntervals)break;
                      mu = get<0>(T[row][col].dataList[k]);
                      sigma = get<1>(T[row][col].dataList[k]);
                      dist_X = abs(sigma - T[i][j].Lright_top.second);
                      if (dist_X < T[i][j].minDistToLtop) {
                          T[i][j].minDistToLtop = dist_X;
                          T[i][j].top_nearest_mu = mu;
                          T[i][j].top_nearest_sigma = sigma;
                          T[i][j].top_t = get<2>(T[row][col].dataList[k]);
                      }
                  }
              }else {
                  // if (T[row][col].minDistToLtop!=MAX_VAL) {
                      mu = T[row][col].top_nearest_mu;
                      sigma = T[row][col].top_nearest_sigma;
                      dist_X = abs(sigma - T[i][j].Lright_top.second);
                      T[i][j].minDistToLtop = dist_X;
                      T[i][j].top_nearest_mu = mu;
                      T[i][j].top_nearest_sigma = sigma;
                      T[i][j].top_t = T[row][col].top_t;
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
                   T[i][j].right_bottom_t = T[row][col].right_bottom_t;
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
                       T[i][j].right_bottom_t = T[row][col].right_bottom_t;
                   }
               // }
               row = MAX(i-1,0);
               col = MIN(j+1,numIntervals - 1);
               //遍历一个格子来给点赋值
               for (int k=0;k<T[row][col].dataList.size();k++) {
                   if (i<1||j+2>numIntervals)break;
                   mu = get<0>(T[row][col].dataList[k]);
                   sigma = get<1>(T[row][col].dataList[k]);
                   dist_X = sqrt(dist(T[i][j].Lright_bottom.first,mu)+dist(T[i][j].Lright_bottom.second,sigma));
                   if (dist_X < T[i][j].right_bottom_nearest_dist) {
                       T[i][j].right_bottom_nearest_dist = dist_X;
                       T[i][j].right_bottom_nearest_mu = mu;
                       T[i][j].right_bottom_nearest_sigma = sigma;
                       T[i][j].right_bottom_t = get<2>(T[row][col].dataList[k]);
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
                   T[i][j].left_top_t = T[row][col].left_top_t;
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
                       T[i][j].left_top_t = T[row][col].left_top_t;
                   }
               // }
               row = MIN(i+1,numIntervals-1);
               col = MAX(j-1,0);
               //遍历一个格子来给点赋值
               for (int k=0;k<T[row][col].dataList.size();k++) {
                   if (j<1||i+2>numIntervals)break;
                   mu = get<0>(T[row][col].dataList[k]);
                   sigma = get<1>(T[row][col].dataList[k]);
                   dist_X = sqrt(dist(T[i][j].Lleft_top.first,mu)+dist(T[i][j].Lleft_top.second,sigma));
                   if (dist_X < T[i][j].left_top_nearest_dist) {
                       T[i][j].left_top_nearest_dist = dist_X;
                       T[i][j].left_top_nearest_mu = mu;
                       T[i][j].left_top_nearest_sigma = sigma;
                       T[i][j].left_top_t = get<2>(T[row][col].dataList[k]);
                   }
               }
           }
       }

       //update
       // for (int i = 0; i < numIntervals; ++i) {
       //     for (int j = 0; j < numIntervals; ++j) {
       //         T[i][j].minDistToLleft = MIN(MIN(T[i][j].minDistToLleft,T[i][j].left_bottom_nearest_dist),T[i][j].left_top_nearest_dist);
       //         T[i][j].minDistToLright = MIN(MIN(T[i][j].minDistToLright,T[i][j].right_bottom_nearest_dist),T[i][j].right_top_nearest_dist);
       //         T[i][j].minDistToLbottom = MIN(MIN(T[i][j].minDistToLbottom,T[i][j].left_bottom_nearest_dist),T[i][j].right_bottom_nearest_dist);
       //         T[i][j].minDistToLtop = MIN(MIN(T[i][j].minDistToLtop,T[i][j].left_top_nearest_dist),T[i][j].right_top_nearest_dist);
       //     }
       // }
      for (int i = 0; i < numIntervals; ++i) {
          for (int j = 0; j < numIntervals; ++j) {
              //left
              if (T[i][j].left_bottom_nearest_dist<T[i][j].minDistToLleft) {
                  T[i][j].minDistToLleft = T[i][j].left_bottom_nearest_dist;
                  T[i][j].left_t = T[i][j].left_bottom_t;
              }
              if (T[i][j].left_top_nearest_dist<T[i][j].minDistToLleft) {
                  T[i][j].minDistToLleft = T[i][j].left_top_nearest_dist;
                  T[i][j].left_t = T[i][j].left_top_t;
              }
              //right
              if (T[i][j].right_bottom_nearest_dist<T[i][j].minDistToLright) {
                  T[i][j].minDistToLright = T[i][j].right_bottom_nearest_dist;
                  T[i][j].right_t = T[i][j].right_bottom_t;
              }
              if (T[i][j].right_top_nearest_dist<T[i][j].minDistToLright) {
                  T[i][j].minDistToLright = T[i][j].right_top_nearest_dist;
                  T[i][j].right_t = T[i][j].right_top_t;
              }
              //bottom
              if (T[i][j].left_bottom_nearest_dist<T[i][j].minDistToLbottom) {
                  T[i][j].minDistToLbottom = T[i][j].left_bottom_nearest_dist;
                  T[i][j].bottom_t = T[i][j].left_bottom_t;
              }
              if (T[i][j].right_bottom_nearest_dist<T[i][j].minDistToLbottom) {
                  T[i][j].minDistToLbottom = T[i][j].right_bottom_nearest_dist;
                  T[i][j].bottom_t = T[i][j].right_bottom_t;
              }
            //top
              if (T[i][j].left_top_nearest_dist<T[i][j].minDistToLtop) {
                  T[i][j].minDistToLtop = T[i][j].left_top_nearest_dist;
                  T[i][j].top_t = T[i][j].left_top_t;
              }
              if (T[i][j].right_top_nearest_dist<T[i][j].minDistToLtop) {
                  T[i][j].minDistToLtop = T[i][j].right_top_nearest_dist;
                  T[i][j].top_t = T[i][j].right_top_t;
              }
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



double CornerDistE(const vector<double>& X, const vector<double>& V, double dNN,vector<char>&mask) {
    double dMu = (X[0] - V[0])*(X[0] - V[0]);
    double dSigma = (X[1] - V[1])*(X[1] - V[1]);
    return dMu + dSigma + dNN*dNN + 2 * dNN * sqrt(MIN(dMu, dSigma));
}


double CenterDist2D(double mu_X, double sigma_X, double L_mu_Q, double U_mu_Q, double L_sigma_Q, double U_sigma_Q, std::vector<vector<Cell>>& T,double intervalLength_mu,double intervalLength_sigma,int lambda,int num,vector<char>&mask) {
      double dist_x = 0.0;
      int res_t=-1;
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
          double epsilon = MAX_VAL;
          if (leftDist<epsilon) {
              epsilon = leftDist;
              res_t = T[row][col].left_t;
          }
          if (rightDist<epsilon) {
              epsilon = rightDist;
              res_t = T[row][col].right_t;
          }
          if (bottomDist<epsilon) {
              epsilon = bottomDist;
              res_t = T[row][col].bottom_t;
          }
          if (topDist<epsilon) {
              epsilon = topDist;
              res_t = T[row][col].top_t;
          }
          if (T[row][col].dataList.size()== 1) {
              double dist_inter = dist(mu_X,get<0>(T[row][col].dataList[0]))+dist(sigma_X,get<1>(T[row][col].dataList[0]));
              if (dist_inter<epsilon) {
                  epsilon = dist_inter;
                  res_t = get<2>(T[row][col].dataList[0]);
              }
          }
          dist_x += epsilon;
        }
    mask[res_t] = 1;
    return dist_x;
}

double CenterDist2D_edge_no_mask(double mu_X, double sigma_X, double L_mu_Q, double U_mu_Q, double L_sigma_Q, double U_sigma_Q, std::vector<vector<Cell>>& T,double intervalLength_mu,double intervalLength_sigma,int lambda,int num) {
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
              epsilon = MIN(epsilon, dist(mu_X,get<0>(T[row][col].dataList[0]))+dist(sigma_X,get<1>(T[row][col].dataList[0])));
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
            return dMu + CenterDist2D_edge_no_mask(V[0],X[1],L_mu_Q,U_mu_Q,L_sigma_Q,U_sigma_Q,T,intervalLength_mu,intervalLength_sigma,lambda,num);
        } else {
            // return dNN*dNN + dMu + dSigma - 2 * dNN * sqrt(dSigma);
            return MAX(dMu + CenterDist2D_edge_no_mask(V[0],X[1],L_mu_Q,U_mu_Q,L_sigma_Q,U_sigma_Q,T,intervalLength_mu,intervalLength_sigma,lambda,num),dNN*dNN + dMu + dSigma - 2 * dNN * sqrt(dSigma));
        }
    } else if (TopOrBottomMiddle) {
        if (dNN*dNN <= dMu) {
            return dSigma + CenterDist2D_edge_no_mask(X[0],V[1],L_mu_Q,U_mu_Q,L_sigma_Q,U_sigma_Q,T,intervalLength_mu,intervalLength_sigma,lambda,num);
            // return dSigma+CenterDist(X[0],L_mu_Q,U_mu_Q,T_mu,intervalLength_mu,lambda,num);
        } else {
            // return dNN*dNN + dMu + dSigma - 2 * dNN * sqrt(dMu);
            return MAX(dSigma + CenterDist2D_edge_no_mask(X[0],V[1],L_mu_Q,U_mu_Q,L_sigma_Q,U_sigma_Q,T,intervalLength_mu,intervalLength_sigma,lambda,num),dNN*dNN + dMu + dSigma - 2 * dNN * sqrt(dMu));
        }
    }

    return -1;
}



vector<vector<pair<double,int>>> Preprocessing_MC3E(vector<vector<double>>& Q,
    vector<double>& mu_Q,vector<double>& sigma_Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,int r) {
    int lambda =3;
    int m = Q.size();
    lower_upper_lemire(mu_Q, m, r, L_mu_Q, U_mu_Q);
    lower_upper_lemire(sigma_Q, m, r, L_sigma_Q, U_sigma_Q);

    // vector<vector<double>> result(4, vector<double>(m));
      vector<vector<pair<double,int>>> result(4, vector<pair<double,int>>(m));


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
        int res_t = start[i];
        for (int t = start[i]; t <= end[i]; ++t) {
            Qt = {mu_Q[t], sigma_Q[t]};
            corner = {L_mu_Q[i], L_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            if (dist<minDbl) {
                minDbl = dist;
                res_t = t;
            }
        }
        result[0][i].first = minDbl;
        result[0][i].second = res_t;
    }
    for (int i = 0; i < m; ++i) {
        minDbr = MAX_VAL;
        int res_t = start[i];
        for (int t = start[i]; t <= end[i]; ++t) {
            Qt = {mu_Q[t], sigma_Q[t]};
            corner = {U_mu_Q[i], L_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            if (dist<minDbr) {
                minDbr = dist;
                res_t = t;
            }
        }
        result[1][i].first = minDbr;
        result[1][i].second = res_t;
    }
    for (int i = 0; i < m; ++i) {
        minDtl = MAX_VAL;
        int res_t = start[i];
        for (int t = start[i]; t <= end[i]; ++t) {
            Qt = {mu_Q[t], sigma_Q[t]};
            corner = {L_mu_Q[i], U_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            if (dist<minDtl) {
                minDtl = dist;
                res_t = t;
            }
        }
        result[2][i].first = minDtl;
        result[2][i].second = res_t;
    }
    for (int i = 0; i < m; ++i) {
        minDtr = MAX_VAL;
        int res_t = start[i];
        for (int t = start[i]; t <= end[i]; ++t) {
            Qt = {mu_Q[t], sigma_Q[t]};
            corner = {U_mu_Q[i], U_sigma_Q[i]};
            dist = euclideanDistance_2(Qt, corner,Qt.size());
            if (dist<minDtr) {
                minDtr = dist;
                res_t = t;
            }
        }
        result[3][i].first = minDtr;
        result[3][i].second = res_t;
    }



    std::pair<double, double> intervalLength;
    //计算Q的m个T表
    for (int i = 0; i < m; ++i) {
        num[i] = end[i] - start[i] + 1;
        pair<vector<double>,vector<double>> Q_part = {mu_Q,sigma_Q};
        intervalLength = initializeT(T[i], Q_part, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i],U_sigma_Q[i],lambda, num[i],start[i],end[i]);
        intervalLength_mu[i]=intervalLength.first;
        intervalLength_sigma[i]=intervalLength.second;
    }

    for (int i = 0; i < result.size(); ++i) {
        for (int j = 0; j < result[i].size(); ++j) {
            result[i][j].first = sqrt(result[i][j].first);
        }
    }

    return result;
}

double LB_GRAPH(vector<vector<double>>& A,vector<vector<double>>& Q,
    vector<double>& U_mu_Q,vector<double>& L_mu_Q,
    vector<double>& U_sigma_Q,vector<double>& L_sigma_Q,
    vector<pair<double,int>>& Dbl, vector<pair<double,int>>& Dbr,
    vector<pair<double,int>>& Dtl, vector<pair<double,int>>& Dtr,
    vector<vector<vector<Cell>>>& T,
    vector <double>& intervalLength_mu,vector <double>& intervalLength_sigma ,vector <int>& num,vector<double>&lb_record,vector<char>&mask,int r) {

    int m = A.size();
    int K = A[0].size();
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

        //corner
        if (Ai[0] < L_mu_Q[i] && Ai[1] < L_sigma_Q[i]) {
            V = {L_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbl[i].first;
            mask[Dbl[i].second] = 1;
            lb_record[i] = CornerDistE(Ai, V, dNN,mask);
            // lb += lb_record[i];
        } else if (Ai[0] > U_mu_Q[i] && Ai[1] < L_sigma_Q[i]) {
            V = {U_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbr[i].first;
            mask[Dbr[i].second] = 1;
            lb_record[i] = CornerDistE(Ai, V, dNN,mask);
            // lb += lb_record[i];
        } else if (Ai[0] < L_mu_Q[i] && Ai[1] > U_sigma_Q[i]) {
            V = {L_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtl[i].first;
            mask[Dtl[i].second] = 1;
            lb_record[i] = CornerDistE(Ai, V, dNN,mask);
            // lb += lb_record[i];
        } else if (Ai[0] > U_mu_Q[i] && Ai[1] > U_sigma_Q[i]){
            V = {U_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtr[i].first;
            mask[Dtr[i].second] = 1;
            lb_record[i] = CornerDistE(Ai, V, dNN,mask);
            // lb += lb_record[i];
        }
        //edge
        else if (Ai[0] <= L_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i]) {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {L_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbl[i].first;
            d_minus = Dtl[i].first;
            MiddleLeftOrRight = true;
            TopOrBottomMiddle=false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            if (delta_minus<delta_plus) {
                lb_record[i] = delta_minus;
                mask[Dtl[i].second] = 1;
            }else {
                lb_record[i] = delta_plus;
                mask[Dbl[i].second] = 1;
            }
            // lb += lb_record[i];
        } else if (Ai[0] >= U_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i]) {
            V_plus = {U_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbr[i].first;
            d_minus = Dtr[i].first;
            MiddleLeftOrRight = true;
            TopOrBottomMiddle=false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            if (delta_minus<delta_plus) {
                lb_record[i] = delta_minus;
                mask[Dtr[i].second] = 1;
            }else {
                lb_record[i] = delta_plus;
                mask[Dbr[i].second] = 1;
            }
            // lb += lb_record[i];
        } else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] <= L_sigma_Q[i]) {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], L_sigma_Q[i]};
            d_plus = Dbl[i].first;
            d_minus = Dbr[i].first;
            MiddleLeftOrRight=false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            if (delta_minus<delta_plus) {
                lb_record[i] = delta_minus;
                mask[Dbr[i].second] = 1;
            }else {
                lb_record[i] = delta_plus;
                mask[Dbl[i].second] = 1;
            }
            // lb += lb_record[i];
        } else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] >= U_sigma_Q[i]){
            V_plus = {L_mu_Q[i], U_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dtl[i].first;
            d_minus = Dtr[i].first;
            MiddleLeftOrRight=false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus,MiddleLeftOrRight,TopOrBottomMiddle,L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i]);
            if (delta_minus<delta_plus) {
                lb_record[i] = delta_minus;
                mask[Dtr[i].second] = 1;
            }else {
                lb_record[i] = delta_plus;
                mask[Dtl[i].second] = 1;
            }
            // lb += lb_record[i];
        }
        //center
         else {
             lb_record[i] = CenterDist2D(mu_A[i],sigma_A[i],L_mu_Q[i],U_mu_Q[i],L_sigma_Q[i],U_sigma_Q[i],T[i],intervalLength_mu[i],intervalLength_sigma[i],lambda,num[i],mask);
             // lb += lb_record[i];
        }
        lb += lb_record[i];
    }
      if (K!=2) {
          lb = K*lb;
          for (auto & x : lb_record) {
              x=K*x;
          }
      }


      // // //计算右侧空缺的点
      vector<double>max_val;
      upper_lemire(lb_record,m,r,max_val);
      for (int i=0;i<m;i++) {
          if (mask[i])continue;
          double delta = 0;
          int start = MAX(0,i-r);
          int end = MIN(m-1,i+r);
          vector<vector<double>> A_part(A.begin() + start, A.begin() + end + 1);
          delta = LB_MV_graph(Q[i],A_part,K);
          lb += MAX(0.0,delta-max_val[i]);
      }

      return sqrt(lb);
}
