
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <utility>
#include "funcs.h"
#include "lb_early.h"
using namespace std;

static bool mc3_lb_exceeds(double lb, int K, double best_dist)
{
    if (!lb_can_early_stop(best_dist))
    {
        return false;
    }
    return static_cast<double>(K) * lb >= lb_threshold_sq(best_dist);
}

/** LB_MC3 / LB_MC3E：统一为 sqrt(K * lb)，不再区分 K=2。 */
static double mc3_finalize(double lb, int K)
{
    return sqrt(static_cast<double>(K) * lb);
}

static void mc3_scale_accumulator(double &lb, std::vector<double> &lb_record, int K)
{
    lb = static_cast<double>(K) * lb;
    for (double &x : lb_record)
    {
        x = static_cast<double>(K) * x;
    }
}

/** lb 已乘 K 后与 best_dist 比较（用于 MC3E 后半段）。 */
static bool mc3_scaled_lb_exceeds(double lb_scaled, double best_dist)
{
    if (!lb_can_early_stop(best_dist))
    {
        return false;
    }
    return lb_scaled >= lb_threshold_sq(best_dist);
}

// 优化版本
pair<double, double> initializeT(std::vector<vector<Cell>> &T, std::pair<vector<double>, vector<double>> &Q, double L_mu_Q, double U_mu_Q, double L_sigma_Q, double U_sigma_Q, int lambda, int n, int start, int end)
{
    int numIntervals = lambda * static_cast<int>(std::sqrt(n));
    double intervalLength_mu = (U_mu_Q - L_mu_Q) / numIntervals;
    double intervalLength_sigma = (U_sigma_Q - L_sigma_Q) / numIntervals;
    T.resize(numIntervals + 1);
    for (int i = 0; i < numIntervals + 1; ++i)
    {
        double sigmaBottom = L_sigma_Q + i * intervalLength_sigma;
        double sigmaTop = sigmaBottom + intervalLength_sigma;
        for (int j = 0; j < numIntervals + 1; ++j)
        {
            double muLeft = L_mu_Q + j * intervalLength_mu;
            double muRight = muLeft + intervalLength_mu;
            T[i].emplace_back(muLeft, muRight, sigmaBottom, sigmaTop);
        }
    }
    T[0][0].numIntervals_sigma = numIntervals;
    T[0][0].numIntervals_mu = numIntervals;
    if (intervalLength_mu != 0)
        intervalLength_mu = 1 / intervalLength_mu;
    if (intervalLength_sigma != 0)
        intervalLength_sigma = 1 / intervalLength_sigma;
    // if (intervalLength_mu==0)intervalLength_mu = 1;
    // if (intervalLength_sigma==0)intervalLength_sigma = 1;
    for (int i = start; i <= end; ++i)
    {
        double mu = Q.first[i];
        double sigma = Q.second[i];
        int muIndex;
        int sigmaIndex;
        if (mu != U_mu_Q)
            muIndex = static_cast<int>((mu - L_mu_Q) * intervalLength_mu);
        else
            muIndex = numIntervals - 1;
        // if (mu!=U_mu_Q) muIndex = static_cast<int>((mu - L_mu_Q) / intervalLength_mu);

        if (sigma != U_sigma_Q)
            sigmaIndex = static_cast<int>((sigma - L_sigma_Q) * intervalLength_sigma);
        else
            sigmaIndex = numIntervals - 1;
        // if (sigma!=U_sigma_Q) sigmaIndex = static_cast<int>((sigma - L_sigma_Q) / intervalLength_sigma);
        T[sigmaIndex][muIndex].dataList.emplace_back(mu, sigma);
        T[sigmaIndex][muIndex].size++;
    }

    // leftbottom point and left,bottom edge
    for (int i = 0; i < numIntervals; ++i)
    {
        for (int j = 0; j < numIntervals; ++j)
        {
            double mu, sigma, dist_X;
            //继承两个点
            int row = i;
            int col = MAX(j - 1, 0);
            // if (T[row][col].left_bottom_nearest_dist!=MAX_VAL) {
            mu = T[row][col].left_bottom_nearest_mu;
            sigma = T[row][col].left_bottom_nearest_sigma;
            dist_X = sqrt(dist(T[i][j].Lleft_bottom.first, mu) + dist(T[i][j].Lleft_bottom.second, sigma));
            T[i][j].left_bottom_nearest_dist = dist_X;
            T[i][j].left_bottom_nearest_mu = mu;
            T[i][j].left_bottom_nearest_sigma = sigma;
            // }
            row = MAX(i - 1, 0);
            col = j;
            // if (T[row][col].left_bottom_nearest_dist!=MAX_VAL) {
            mu = T[row][col].left_bottom_nearest_mu;
            sigma = T[row][col].left_bottom_nearest_sigma;
            dist_X = sqrt(dist(T[i][j].Lleft_bottom.first, mu) + dist(T[i][j].Lleft_bottom.second, sigma));
            if (dist_X < T[i][j].left_bottom_nearest_dist)
            {
                T[i][j].left_bottom_nearest_dist = dist_X;
                T[i][j].left_bottom_nearest_mu = mu;
                T[i][j].left_bottom_nearest_sigma = sigma;
            }
            // }
            row = MAX(i - 1, 0);
            col = MAX(j - 1, 0);
            //遍历一个格子来给点赋值
            for (int k = 0; k < T[row][col].size; k++)
            {
                if (i < 1 || j < 1)
                    break;
                mu = T[row][col].dataList[k].first;
                sigma = T[row][col].dataList[k].second;
                dist_X = sqrt(dist(T[i][j].Lleft_bottom.first, mu) + dist(T[i][j].Lleft_bottom.second, sigma));
                if (dist_X < T[i][j].left_bottom_nearest_dist)
                {
                    T[i][j].left_bottom_nearest_dist = dist_X;
                    T[i][j].left_bottom_nearest_mu = mu;
                    T[i][j].left_bottom_nearest_sigma = sigma;
                }
            }
            //遍历一个格子来给left edge赋值
            row = i;
            col = MAX(0, j - 1);
            if (T[row][col].size != 0)
            {
                for (int k = 0; k < T[row][col].size; k++)
                {
                    if (j < 1)
                        break;
                    mu = T[row][col].dataList[k].first;
                    sigma = T[row][col].dataList[k].second;
                    dist_X = abs(T[i][j].Lleft_bottom.first - mu);
                    if (dist_X < T[i][j].minDistToLleft)
                    {
                        T[i][j].minDistToLleft = dist_X;
                        T[i][j].left_nearest_mu = mu;
                        T[i][j].left_nearest_sigma = sigma;
                    }
                }
            }
            else
            {
                // if (T[row][col].minDistToLleft!=MAX_VAL) {
                mu = T[row][col].left_nearest_mu;
                sigma = T[row][col].left_nearest_sigma;
                dist_X = abs(T[i][j].Lleft_bottom.first - mu);
                T[i][j].minDistToLleft = dist_X;
                T[i][j].left_nearest_mu = mu;
                T[i][j].left_nearest_sigma = sigma;
                // }
            }
            row = MAX(i - 1, 0);
            col = j;
            //遍历一个格子来给bottom edge赋值
            if (T[row][col].size != 0)
            {
                for (int k = 0; k < T[row][col].size; k++)
                {
                    if (i < 1)
                        break;
                    mu = T[row][col].dataList[k].first;
                    sigma = T[row][col].dataList[k].second;
                    dist_X = abs(T[i][j].Lright_bottom.second - sigma);
                    if (dist_X < T[i][j].minDistToLbottom)
                    {
                        T[i][j].minDistToLbottom = dist_X;
                        T[i][j].bottom_nearest_mu = mu;
                        T[i][j].bottom_nearest_sigma = sigma;
                    }
                }
            }
            else
            {
                // if (T[row][col].minDistToLbottom!=MAX_VAL) {
                mu = T[row][col].bottom_nearest_mu;
                sigma = T[row][col].bottom_nearest_sigma;
                dist_X = abs(T[i][j].Lright_bottom.second - sigma);
                T[i][j].minDistToLbottom = dist_X;
                T[i][j].bottom_nearest_mu = mu;
                T[i][j].bottom_nearest_sigma = sigma;
                // }
            }
        }
    }

    // righttop point and right,top edge
    for (int i = numIntervals - 1; i >= 0; --i)
    {
        for (int j = numIntervals - 1; j >= 0; --j)
        {
            double mu, sigma, dist_X;
            //继承两个点
            int row = MIN(i + 1, numIntervals - 1);
            int col = j;
            // if (T[row][col].right_top_nearest_dist!=MAX_VAL) {
            mu = T[row][col].right_top_nearest_mu;
            sigma = T[row][col].right_top_nearest_sigma;
            dist_X = sqrt(dist(T[i][j].Lright_top.first, mu) + dist(T[i][j].Lright_top.second, sigma));
            T[i][j].right_top_nearest_dist = dist_X;
            T[i][j].right_top_nearest_mu = mu;
            T[i][j].right_top_nearest_sigma = sigma;
            // }
            row = i;
            col = MIN(j + 1, numIntervals - 1);
            // if (T[row][col].right_top_nearest_dist!=MAX_VAL) {
            mu = T[row][col].right_top_nearest_mu;
            sigma = T[row][col].right_top_nearest_sigma;
            dist_X = sqrt(dist(T[i][j].Lright_top.first, mu) + dist(T[i][j].Lright_top.second, sigma));
            if (dist_X < T[i][j].right_top_nearest_dist)
            {
                T[i][j].right_top_nearest_dist = dist_X;
                T[i][j].right_top_nearest_mu = mu;
                T[i][j].right_top_nearest_sigma = sigma;
            }
            // }
            row = MIN(i + 1, numIntervals - 1);
            col = MIN(j + 1, numIntervals - 1);
            //遍历一个格子来给点赋值
            for (int k = 0; k < T[row][col].size; k++)
            {
                if (i + 2 > numIntervals || j + 2 > numIntervals)
                    break;
                mu = T[row][col].dataList[k].first;
                sigma = T[row][col].dataList[k].second;
                dist_X = sqrt(dist(T[i][j].Lright_top.first, mu) + dist(T[i][j].Lright_top.second, sigma));
                if (dist_X < T[i][j].right_top_nearest_dist)
                {
                    T[i][j].right_top_nearest_dist = dist_X;
                    T[i][j].right_top_nearest_mu = mu;
                    T[i][j].right_top_nearest_sigma = sigma;
                }
            }
            //遍历一个格子来给right edge赋值
            row = i;
            col = MIN(j + 1, numIntervals - 1);
            if (T[row][col].size != 0)
            {
                for (int k = 0; k < T[row][col].size; k++)
                {
                    if (j + 2 > numIntervals)
                        break;
                    mu = T[row][col].dataList[k].first;
                    sigma = T[row][col].dataList[k].second;
                    dist_X = abs(mu - T[i][j].Lright_bottom.first);
                    if (dist_X < T[i][j].minDistToLright)
                    {
                        T[i][j].minDistToLright = dist_X;
                        T[i][j].right_nearest_mu = mu;
                        T[i][j].right_nearest_sigma = sigma;
                    }
                }
            }
            else
            {
                // if (T[row][col].minDistToLright!=MAX_VAL) {
                mu = T[row][col].right_nearest_mu;
                sigma = T[row][col].right_nearest_sigma;
                dist_X = abs(mu - T[i][j].Lright_bottom.first);
                T[i][j].minDistToLright = dist_X;
                T[i][j].right_nearest_mu = mu;
                T[i][j].right_nearest_sigma = sigma;
                // }
            }
            row = MIN(i + 1, numIntervals - 1);
            col = j;
            //遍历一个格子来给top edge赋值
            if (T[row][col].size != 0)
            {
                for (int k = 0; k < T[row][col].size; k++)
                {
                    if (i + 2 > numIntervals)
                        break;
                    mu = T[row][col].dataList[k].first;
                    sigma = T[row][col].dataList[k].second;
                    dist_X = abs(sigma - T[i][j].Lright_top.second);
                    if (dist_X < T[i][j].minDistToLtop)
                    {
                        T[i][j].minDistToLtop = dist_X;
                        T[i][j].top_nearest_mu = mu;
                        T[i][j].top_nearest_sigma = sigma;
                    }
                }
            }
            else
            {
                // if (T[row][col].minDistToLtop!=MAX_VAL) {
                mu = T[row][col].top_nearest_mu;
                sigma = T[row][col].top_nearest_sigma;
                dist_X = abs(sigma - T[i][j].Lright_top.second);
                T[i][j].minDistToLtop = dist_X;
                T[i][j].top_nearest_mu = mu;
                T[i][j].top_nearest_sigma = sigma;
                // }
            }
        }
    }

    // rightbottom point
    for (int i = 0; i < numIntervals; ++i)
    {
        for (int j = numIntervals - 1; j >= 0; --j)
        {
            double mu, sigma, dist_X;
            //继承两个点
            int row = i;
            int col = MIN(j + 1, numIntervals - 1);
            // if (T[row][col].right_bottom_nearest_dist!=MAX_VAL) {
            mu = T[row][col].right_bottom_nearest_mu;
            sigma = T[row][col].right_bottom_nearest_sigma;
            dist_X = sqrt(dist(T[i][j].Lright_bottom.first, mu) + dist(T[i][j].Lright_bottom.second, sigma));
            T[i][j].right_bottom_nearest_dist = dist_X;
            T[i][j].right_bottom_nearest_mu = mu;
            T[i][j].right_bottom_nearest_sigma = sigma;
            // }
            row = MAX(i - 1, 0);
            col = j;
            // if (T[row][col].right_bottom_nearest_dist!=MAX_VAL) {
            mu = T[row][col].right_bottom_nearest_mu;
            sigma = T[row][col].right_bottom_nearest_sigma;
            dist_X = sqrt(dist(T[i][j].Lright_bottom.first, mu) + dist(T[i][j].Lright_bottom.second, sigma));
            if (dist_X < T[i][j].right_bottom_nearest_dist)
            {
                T[i][j].right_bottom_nearest_dist = dist_X;
                T[i][j].right_bottom_nearest_mu = mu;
                T[i][j].right_bottom_nearest_sigma = sigma;
            }
            // }
            row = MAX(i - 1, 0);
            col = MIN(j + 1, numIntervals - 1);
            //遍历一个格子来给点赋值
            for (int k = 0; k < T[row][col].size; k++)
            {
                if (i < 1 || j + 2 > numIntervals)
                    break;
                mu = T[row][col].dataList[k].first;
                sigma = T[row][col].dataList[k].second;
                dist_X = sqrt(dist(T[i][j].Lright_bottom.first, mu) + dist(T[i][j].Lright_bottom.second, sigma));
                if (dist_X < T[i][j].right_bottom_nearest_dist)
                {
                    T[i][j].right_bottom_nearest_dist = dist_X;
                    T[i][j].right_bottom_nearest_mu = mu;
                    T[i][j].right_bottom_nearest_sigma = sigma;
                }
            }
        }
    }

    // lefttop point
    for (int i = numIntervals - 1; i >= 0; --i)
    {
        for (int j = 0; j < numIntervals; ++j)
        {
            double mu, sigma, dist_X;
            //继承两个点
            int row = i;
            int col = MAX(j - 1, 0);
            // if (T[row][col].left_top_nearest_dist!=MAX_VAL) {
            mu = T[row][col].left_top_nearest_mu;
            sigma = T[row][col].left_top_nearest_sigma;
            dist_X = sqrt(dist(T[i][j].Lleft_top.first, mu) + dist(T[i][j].Lleft_top.second, sigma));
            T[i][j].left_top_nearest_dist = dist_X;
            T[i][j].left_top_nearest_mu = mu;
            T[i][j].left_top_nearest_sigma = sigma;
            // }
            row = MIN(i + 1, numIntervals - 1);
            col = j;
            // if (T[row][col].left_top_nearest_dist!=MAX_VAL) {
            mu = T[row][col].left_top_nearest_mu;
            sigma = T[row][col].left_top_nearest_sigma;
            dist_X = sqrt(dist(T[i][j].Lleft_top.first, mu) + dist(T[i][j].Lleft_top.second, sigma));
            if (dist_X < T[i][j].left_top_nearest_dist)
            {
                T[i][j].left_top_nearest_dist = dist_X;
                T[i][j].left_top_nearest_mu = mu;
                T[i][j].left_top_nearest_sigma = sigma;
            }
            // }
            row = MIN(i + 1, numIntervals - 1);
            col = MAX(j - 1, 0);
            //遍历一个格子来给点赋值
            for (int k = 0; k < T[row][col].size; k++)
            {
                if (j < 1 || i + 2 > numIntervals)
                    break;
                mu = T[row][col].dataList[k].first;
                sigma = T[row][col].dataList[k].second;
                dist_X = sqrt(dist(T[i][j].Lleft_top.first, mu) + dist(T[i][j].Lleft_top.second, sigma));
                if (dist_X < T[i][j].left_top_nearest_dist)
                {
                    T[i][j].left_top_nearest_dist = dist_X;
                    T[i][j].left_top_nearest_mu = mu;
                    T[i][j].left_top_nearest_sigma = sigma;
                }
            }
        }
    }

    // update
    //  for (int i = 0; i < numIntervals; ++i) {
    //      for (int j = 0; j < numIntervals; ++j) {
    //          T[i][j].minDistToLleft = MIN(MIN(T[i][j].minDistToLleft,T[i][j].left_bottom_nearest_dist),T[i][j].left_top_nearest_dist);
    //          T[i][j].minDistToLright = MIN(MIN(T[i][j].minDistToLright,T[i][j].right_bottom_nearest_dist),T[i][j].right_top_nearest_dist);
    //          T[i][j].minDistToLbottom = MIN(MIN(T[i][j].minDistToLbottom,T[i][j].left_bottom_nearest_dist),T[i][j].right_bottom_nearest_dist);
    //          T[i][j].minDistToLtop = MIN(MIN(T[i][j].minDistToLtop,T[i][j].left_top_nearest_dist),T[i][j].right_top_nearest_dist);
    //      }
    //  }
    for (int i = 0; i < numIntervals; ++i)
    {
        for (int j = 0; j < numIntervals; ++j)
        {
            // left
            if (T[i][j].left_bottom_nearest_dist < T[i][j].minDistToLleft)
            {
                T[i][j].minDistToLleft = T[i][j].left_bottom_nearest_dist;
            }
            if (T[i][j].left_top_nearest_dist < T[i][j].minDistToLleft)
            {
                T[i][j].minDistToLleft = T[i][j].left_top_nearest_dist;
            }
            // right
            if (T[i][j].right_bottom_nearest_dist < T[i][j].minDistToLright)
            {
                T[i][j].minDistToLright = T[i][j].right_bottom_nearest_dist;
            }
            if (T[i][j].right_top_nearest_dist < T[i][j].minDistToLright)
            {
                T[i][j].minDistToLright = T[i][j].right_top_nearest_dist;
            }
            // bottom
            if (T[i][j].left_bottom_nearest_dist < T[i][j].minDistToLbottom)
            {
                T[i][j].minDistToLbottom = T[i][j].left_bottom_nearest_dist;
            }
            if (T[i][j].right_bottom_nearest_dist < T[i][j].minDistToLbottom)
            {
                T[i][j].minDistToLbottom = T[i][j].right_bottom_nearest_dist;
            }
            // top
            if (T[i][j].left_top_nearest_dist < T[i][j].minDistToLtop)
            {
                T[i][j].minDistToLtop = T[i][j].left_top_nearest_dist;
            }
            if (T[i][j].right_top_nearest_dist < T[i][j].minDistToLtop)
            {
                T[i][j].minDistToLtop = T[i][j].right_top_nearest_dist;
            }
        }
    }

    // copy
    for (int i = 0; i < numIntervals; ++i)
    {
        T[i][numIntervals] = T[i][numIntervals - 1];
    }
    for (int j = 0; j < numIntervals; ++j)
    {
        T[numIntervals][j] = T[numIntervals - 1][j];
    }
    T[numIntervals][numIntervals] = T[numIntervals - 1][numIntervals - 1];

    return {intervalLength_mu, intervalLength_sigma};
}

double CornerDistE(const vector<double> &X, const vector<double> &V, double dNN)
{
    double dMu = (X[0] - V[0]) * (X[0] - V[0]);
    double dSigma = (X[1] - V[1]) * (X[1] - V[1]);
    return dMu + dSigma + dNN * dNN + 2 * dNN * sqrt(MIN(dMu, dSigma));
}

double CenterDist2D(double mu_X, double sigma_X, double L_mu_Q, double U_mu_Q, double L_sigma_Q, double U_sigma_Q, std::vector<vector<Cell>> &T, double intervalLength_mu, double intervalLength_sigma, int lambda, int num)
{
    double dist_x = 0.0;
    int col = round((mu_X - L_mu_Q) * intervalLength_mu);
    int row = round((sigma_X - L_sigma_Q) * intervalLength_sigma);
    // int col = round((mu_X - L_mu_Q) / intervalLength_mu);
    // int row = round((sigma_X - L_sigma_Q) / intervalLength_sigma);
    if (T[row][col].size >= 2)
    {
        return 0;
    }
    else
    {
        double leftDist = dist(T[row][col].minDistToLleft + mu_X, T[row][col].Lleft_bottom.first);
        double rightDist = dist(T[row][col].minDistToLright + T[row][col].Lright_bottom.first, mu_X);
        double bottomDist = dist(T[row][col].minDistToLbottom + sigma_X, T[row][col].Lleft_bottom.second);
        double topDist = dist(T[row][col].minDistToLtop + T[row][col].Lleft_top.second, sigma_X);
        dist_x = MIN(MIN(MIN(leftDist, rightDist), bottomDist), topDist);
        if (T[row][col].size == 1)
        {
            double dist_inter = dist(mu_X, T[row][col].dataList[0].first) + dist(sigma_X, T[row][col].dataList[0].second);
            dist_x = MIN(dist_x, dist_inter);
        }
    }

    return dist_x;
}

double EdgeDistE(const vector<double> &X, const vector<double> &V, double dNN, bool MiddleLeftOrRight, bool TopOrBottomMiddle, double L_mu_Q, double U_mu_Q,
                 double L_sigma_Q, double U_sigma_Q, vector<vector<Cell>> &T, double intervalLength_mu, double intervalLength_sigma, int lambda, int num)
{
    double dMu = (X[0] - V[0]) * (X[0] - V[0]);
    double dSigma = (X[1] - V[1]) * (X[1] - V[1]);
    if (MiddleLeftOrRight)
    {
        if (dNN * dNN <= dSigma)
        {
            return dMu + CenterDist2D(V[0], X[1], L_mu_Q, U_mu_Q, L_sigma_Q, U_sigma_Q, T, intervalLength_mu, intervalLength_sigma, lambda, num);
        }
        else
        {
            // return dNN*dNN + dMu + dSigma - 2 * dNN * sqrt(dSigma);
            return MAX(dMu + CenterDist2D(V[0], X[1], L_mu_Q, U_mu_Q, L_sigma_Q, U_sigma_Q, T, intervalLength_mu, intervalLength_sigma, lambda, num), dNN * dNN + dMu + dSigma - 2 * dNN * sqrt(dSigma));
        }
    }
    else if (TopOrBottomMiddle)
    {
        if (dNN * dNN <= dMu)
        {
            return dSigma + CenterDist2D(X[0], V[1], L_mu_Q, U_mu_Q, L_sigma_Q, U_sigma_Q, T, intervalLength_mu, intervalLength_sigma, lambda, num);
            // return dSigma+CenterDist(X[0],L_mu_Q,U_mu_Q,T_mu,intervalLength_mu,lambda,num);
        }
        else
        {
            // return dNN*dNN + dMu + dSigma - 2 * dNN * sqrt(dMu);
            return MAX(dSigma + CenterDist2D(X[0], V[1], L_mu_Q, U_mu_Q, L_sigma_Q, U_sigma_Q, T, intervalLength_mu, intervalLength_sigma, lambda, num), dNN * dNN + dMu + dSigma - 2 * dNN * sqrt(dMu));
        }
    }

    return -1;
}

vector<vector<double>> Preprocessing_MC3(vector<vector<double>> &Q,
                                         vector<double> &mu_Q, vector<double> &sigma_Q,
                                         vector<double> &U_mu_Q, vector<double> &L_mu_Q,
                                         vector<double> &U_sigma_Q, vector<double> &L_sigma_Q,
                                         vector<vector<vector<Cell>>> &T,
                                         vector<double> &intervalLength_mu, vector<double> &intervalLength_sigma, vector<int> &num, int r)
{
    int lambda = 3;
    int m = Q.size();
    lower_upper_lemire(mu_Q, m, r, L_mu_Q, U_mu_Q);
    lower_upper_lemire(sigma_Q, m, r, L_sigma_Q, U_sigma_Q);

    vector<vector<double>> result(4, vector<double>(m));

    double minDbl, minDbr, minDtl, minDtr;
    vector<int> start(m);
    vector<int> end(m);
    vector<double> Qt(2);
    vector<double> corner_bl(2);
    vector<double> corner_br(2);
    vector<double> corner_tl(2);
    vector<double> corner_tr(2);
    double dist_bl;
    double dist_br;
    double dist_tl;
    double dist_tr;

    for (int i = 0; i < m; ++i)
    {
        start[i] = MAX(0, i - r);
        end[i] = MIN(m - 1, i + r);
    }
    //计算四个dNN预选值
    for (int i = 0; i < m; ++i)
    {
        minDbl = MAX_VAL;
        minDbr = MAX_VAL;
        minDtl = MAX_VAL;
        minDtr = MAX_VAL;
        for (int t = start[i]; t <= end[i]; ++t)
        {
            Qt = {mu_Q[t], sigma_Q[t]};
            corner_bl = {L_mu_Q[i], L_sigma_Q[i]};
            corner_br = {U_mu_Q[i], L_sigma_Q[i]};
            corner_tl = {L_mu_Q[i], U_sigma_Q[i]};
            corner_tr = {U_mu_Q[i], U_sigma_Q[i]};
            dist_bl = euclideanDistance_2(Qt, corner_bl);
            dist_br = euclideanDistance_2(Qt, corner_br);
            dist_tl = euclideanDistance_2(Qt, corner_tl);
            dist_tr = euclideanDistance_2(Qt, corner_tr);
            minDbl = MIN(dist_bl, minDbl);
            minDbr = MIN(dist_br, minDbr);
            minDtl = MIN(dist_tl, minDtl);
            minDtr = MIN(dist_tr, minDtr);
        }
        result[0][i] = minDbl;
        result[1][i] = minDbr;
        result[2][i] = minDtl;
        result[3][i] = minDtr;
    }

    std::pair<double, double> intervalLength;
    //计算Q的m个T表
    for (int i = 0; i < m; ++i)
    {
        num[i] = end[i] - start[i] + 1;
        pair<vector<double>, vector<double>> Q_part = {mu_Q, sigma_Q};
        intervalLength = initializeT(T[i], Q_part, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], lambda, num[i], start[i], end[i]);
        intervalLength_mu[i] = intervalLength.first;
        intervalLength_sigma[i] = intervalLength.second;
    }

    for (int i = 0; i < result.size(); ++i)
    {
        for (int j = 0; j < result[i].size(); ++j)
        {
            result[i][j] = sqrt(result[i][j]);
        }
    }

    return result;
}

double LB_MC3(vector<vector<double>> &A,
              vector<double> &U_mu_Q, vector<double> &L_mu_Q,
              vector<double> &U_sigma_Q, vector<double> &L_sigma_Q,
              vector<double> &mu_A, vector<double> &sigma_A,
              vector<double> &Dbl, vector<double> &Dbr,
              vector<double> &Dtl, vector<double> &Dtr,
              vector<vector<vector<Cell>>> &T,
              vector<double> &intervalLength_mu, vector<double> &intervalLength_sigma, vector<int> &num)
{
    int m = A.size();
    int K = A[0].size();
    int lambda = 3;
    double lb = 0.0;

    vector<double> Ai(2);
    vector<double> V(2);
    double dNN;
    vector<double> V_plus(2);
    vector<double> V_minus(2);
    double d_plus;
    double d_minus;
    bool MiddleLeftOrRight = false;
    bool TopOrBottomMiddle = false;

    for (int i = 0; i < m; ++i)
    {
        Ai = {mu_A[i], sigma_A[i]};

        // concer
        if (Ai[0] < L_mu_Q[i] && Ai[1] < L_sigma_Q[i])
        {
            V = {L_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbl[i];
            lb += CornerDistE(Ai, V, dNN);
        }
        else if (Ai[0] > U_mu_Q[i] && Ai[1] < L_sigma_Q[i])
        {
            V = {U_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbr[i];
            lb += CornerDistE(Ai, V, dNN);
        }
        else if (Ai[0] < L_mu_Q[i] && Ai[1] > U_sigma_Q[i])
        {
            V = {L_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtl[i];
            lb += CornerDistE(Ai, V, dNN);
        }
        else if (Ai[0] > U_mu_Q[i] && Ai[1] > U_sigma_Q[i])
        {
            V = {U_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtr[i];
            lb += CornerDistE(Ai, V, dNN);
        }
        // edge
        else if (Ai[0] <= L_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i])
        {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {L_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbl[i];
            d_minus = Dtl[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle = false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            lb += MIN(delta_plus, delta_minus);
        }
        else if (Ai[0] >= U_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i])
        {
            V_plus = {U_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbr[i];
            d_minus = Dtr[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle = false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            lb += MIN(delta_plus, delta_minus);
        }
        else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] <= L_sigma_Q[i])
        {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], L_sigma_Q[i]};
            d_plus = Dbl[i];
            d_minus = Dbr[i];
            MiddleLeftOrRight = false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            lb += MIN(delta_plus, delta_minus);
        }
        else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] >= U_sigma_Q[i])
        {
            V_plus = {L_mu_Q[i], U_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dtl[i];
            d_minus = Dtr[i];
            MiddleLeftOrRight = false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            lb += MIN(delta_plus, delta_minus);
        }
        // center
        else
        {
            lb += CenterDist2D(mu_A[i], sigma_A[i], L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
        }
    }
    return mc3_finalize(lb, K);
}

double LB_MC3E(vector<vector<double>> &A, vector<vector<double>> &Q, vector<vector<double>> &A_T, vector<vector<double>> &Q_T,
               vector<double> &U_mu_Q, vector<double> &L_mu_Q,
               vector<double> &U_sigma_Q, vector<double> &L_sigma_Q,
               vector<double> &mu_A, vector<double> &sigma_A,
               vector<double> &Dbl, vector<double> &Dbr,
               vector<double> &Dtl, vector<double> &Dtr,
               vector<vector<vector<Cell>>> &T,
               vector<double> &intervalLength_mu, vector<double> &intervalLength_sigma, vector<int> &num, vector<double> &lb_record, int r)
{

    int m = A.size();
    int K = A[0].size();
    int lambda = 3;
    double lb = 0.0;

    vector<double> Ai(2);
    vector<double> V(2);
    double dNN;
    vector<double> V_plus(2);
    vector<double> V_minus(2);
    double d_plus;
    double d_minus;
    bool MiddleLeftOrRight = false;
    bool TopOrBottomMiddle = false;
    vector<double> A_L(m), A_U(m);

    for (int i = 0; i < m; ++i)
    {
        Ai = {mu_A[i], sigma_A[i]};

        // corner
        if (Ai[0] < L_mu_Q[i] && Ai[1] < L_sigma_Q[i])
        {
            V = {L_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbl[i];
            lb_record[i] = CornerDistE(Ai, V, dNN);
        }
        else if (Ai[0] > U_mu_Q[i] && Ai[1] < L_sigma_Q[i])
        {
            V = {U_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbr[i];
            lb_record[i] = CornerDistE(Ai, V, dNN);
        }
        else if (Ai[0] < L_mu_Q[i] && Ai[1] > U_sigma_Q[i])
        {
            V = {L_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtl[i];
            lb_record[i] = CornerDistE(Ai, V, dNN);
        }
        else if (Ai[0] > U_mu_Q[i] && Ai[1] > U_sigma_Q[i])
        {
            V = {U_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtr[i];
            lb_record[i] = CornerDistE(Ai, V, dNN);
        }
        // edge
        else if (Ai[0] <= L_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i])
        {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {L_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbl[i];
            d_minus = Dtl[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle = false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            if (delta_minus < delta_plus)
            {
                lb_record[i] = delta_minus;
            }
            else
            {
                lb_record[i] = delta_plus;
            }
        }
        else if (Ai[0] >= U_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i])
        {
            V_plus = {U_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbr[i];
            d_minus = Dtr[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle = false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            if (delta_minus < delta_plus)
            {
                lb_record[i] = delta_minus;
            }
            else
            {
                lb_record[i] = delta_plus;
            }
        }
        else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] <= L_sigma_Q[i])
        {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], L_sigma_Q[i]};
            d_plus = Dbl[i];
            d_minus = Dbr[i];
            MiddleLeftOrRight = false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            if (delta_minus < delta_plus)
            {
                lb_record[i] = delta_minus;
            }
            else
            {
                lb_record[i] = delta_plus;
            }
        }
        else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] >= U_sigma_Q[i])
        {
            V_plus = {L_mu_Q[i], U_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dtl[i];
            d_minus = Dtr[i];
            MiddleLeftOrRight = false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle, L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            if (delta_minus < delta_plus)
            {
                lb_record[i] = delta_minus;
            }
            else
            {
                lb_record[i] = delta_plus;
            }
        }
        // center
        else
        {
            lb_record[i] = CenterDist2D(mu_A[i], sigma_A[i], L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
        }
        lb += lb_record[i];
    }
    mc3_scale_accumulator(lb, lb_record, K);

    // // //计算右侧空缺的点
    vector<double> max_val(m);
    upper_lemire(lb_record, m, r, max_val);

    //计算右侧的全部MV
    vector<double> mv_right(m, 0);
    for (int k = 0; k < K; k++)
    {
        lower_upper_lemire(A_T[k], m, r, A_L, A_U);
        for (int i = 0; i < m; i++)
        {
            if (Q_T[k][i] < A_L[i])
            {
                mv_right[i] += dist(Q_T[k][i], A_L[i]);
            }
            else if (Q_T[k][i] > A_U[i])
            {
                mv_right[i] += dist(Q_T[k][i], A_U[i]);
            }
        }
    }

    for (int i = 0; i < m; i++)
    {
        // cout<<mv_right[i]<<endl;
        lb += MAX(0.0, mv_right[i] - max_val[i]);
    }

    return sqrt(lb);
}

double LB_MC3_early(vector<vector<double>> &A,
                    vector<double> &U_mu_Q, vector<double> &L_mu_Q,
                    vector<double> &U_sigma_Q, vector<double> &L_sigma_Q,
                    vector<double> &mu_A, vector<double> &sigma_A,
                    vector<double> &Dbl, vector<double> &Dbr,
                    vector<double> &Dtl, vector<double> &Dtr,
                    vector<vector<vector<Cell>>> &T,
                    vector<double> &intervalLength_mu, vector<double> &intervalLength_sigma, vector<int> &num,
                    double best_dist)
{
    if (!lb_can_early_stop(best_dist))
    {
        return LB_MC3(A, U_mu_Q, L_mu_Q, U_sigma_Q, L_sigma_Q, mu_A, sigma_A,
                      Dbl, Dbr, Dtl, Dtr, T, intervalLength_mu, intervalLength_sigma, num);
    }

    int m = static_cast<int>(A.size());
    int K = static_cast<int>(A[0].size());
    int lambda = 3;
    double lb = 0.0;
    vector<double> Ai(2);
    vector<double> V(2);
    double dNN;
    vector<double> V_plus(2);
    vector<double> V_minus(2);
    double d_plus;
    double d_minus;
    bool MiddleLeftOrRight = false;
    bool TopOrBottomMiddle = false;

    for (int i = 0; i < m; ++i)
    {
        Ai = {mu_A[i], sigma_A[i]};
        if (Ai[0] < L_mu_Q[i] && Ai[1] < L_sigma_Q[i])
        {
            V = {L_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbl[i];
            lb += CornerDistE(Ai, V, dNN);
        }
        else if (Ai[0] > U_mu_Q[i] && Ai[1] < L_sigma_Q[i])
        {
            V = {U_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbr[i];
            lb += CornerDistE(Ai, V, dNN);
        }
        else if (Ai[0] < L_mu_Q[i] && Ai[1] > U_sigma_Q[i])
        {
            V = {L_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtl[i];
            lb += CornerDistE(Ai, V, dNN);
        }
        else if (Ai[0] > U_mu_Q[i] && Ai[1] > U_sigma_Q[i])
        {
            V = {U_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtr[i];
            lb += CornerDistE(Ai, V, dNN);
        }
        else if (Ai[0] <= L_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i])
        {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {L_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbl[i];
            d_minus = Dtl[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle = false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle,
                                          L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i],
                                          intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle,
                                           L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i],
                                           intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            lb += MIN(delta_plus, delta_minus);
        }
        else if (Ai[0] >= U_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i])
        {
            V_plus = {U_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbr[i];
            d_minus = Dtr[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle = false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle,
                                          L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i],
                                          intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle,
                                           L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i],
                                           intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            lb += MIN(delta_plus, delta_minus);
        }
        else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] <= L_sigma_Q[i])
        {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], L_sigma_Q[i]};
            d_plus = Dbl[i];
            d_minus = Dbr[i];
            MiddleLeftOrRight = false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle,
                                          L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i],
                                          intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle,
                                           L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i],
                                           intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            lb += MIN(delta_plus, delta_minus);
        }
        else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] >= U_sigma_Q[i])
        {
            V_plus = {L_mu_Q[i], U_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dtl[i];
            d_minus = Dtr[i];
            MiddleLeftOrRight = false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle,
                                          L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i],
                                          intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle,
                                           L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i],
                                           intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            lb += MIN(delta_plus, delta_minus);
        }
        else
        {
            lb += CenterDist2D(mu_A[i], sigma_A[i], L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i],
                               T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
        }
        if (mc3_lb_exceeds(lb, K, best_dist))
        {
            return mc3_finalize(lb, K);
        }
    }
    return mc3_finalize(lb, K);
}

double LB_MC3E_early(vector<vector<double>> &A, vector<vector<double>> &Q, vector<vector<double>> &A_T,
                     vector<vector<double>> &Q_T,
                     vector<double> &U_mu_Q, vector<double> &L_mu_Q,
                     vector<double> &U_sigma_Q, vector<double> &L_sigma_Q,
                     vector<double> &mu_A, vector<double> &sigma_A,
                     vector<double> &Dbl, vector<double> &Dbr,
                     vector<double> &Dtl, vector<double> &Dtr,
                     vector<vector<vector<Cell>>> &T,
                     vector<double> &intervalLength_mu, vector<double> &intervalLength_sigma, vector<int> &num,
                     vector<double> &lb_record, int r, double best_dist)
{
    if (!lb_can_early_stop(best_dist))
    {
        return LB_MC3E(A, Q, A_T, Q_T, U_mu_Q, L_mu_Q, U_sigma_Q, L_sigma_Q, mu_A, sigma_A,
                       Dbl, Dbr, Dtl, Dtr, T, intervalLength_mu, intervalLength_sigma, num, lb_record, r);
    }

    int m = static_cast<int>(A.size());
    int K = static_cast<int>(A[0].size());
    int lambda = 3;
    double lb = 0.0;
    vector<double> Ai(2);
    vector<double> V(2);
    double dNN;
    vector<double> V_plus(2);
    vector<double> V_minus(2);
    double d_plus;
    double d_minus;
    bool MiddleLeftOrRight = false;
    bool TopOrBottomMiddle = false;
    vector<double> A_L(m), A_U(m);

    for (int i = 0; i < m; ++i)
    {
        Ai = {mu_A[i], sigma_A[i]};
        if (Ai[0] < L_mu_Q[i] && Ai[1] < L_sigma_Q[i])
        {
            V = {L_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbl[i];
            lb_record[i] = CornerDistE(Ai, V, dNN);
        }
        else if (Ai[0] > U_mu_Q[i] && Ai[1] < L_sigma_Q[i])
        {
            V = {U_mu_Q[i], L_sigma_Q[i]};
            dNN = Dbr[i];
            lb_record[i] = CornerDistE(Ai, V, dNN);
        }
        else if (Ai[0] < L_mu_Q[i] && Ai[1] > U_sigma_Q[i])
        {
            V = {L_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtl[i];
            lb_record[i] = CornerDistE(Ai, V, dNN);
        }
        else if (Ai[0] > U_mu_Q[i] && Ai[1] > U_sigma_Q[i])
        {
            V = {U_mu_Q[i], U_sigma_Q[i]};
            dNN = Dtr[i];
            lb_record[i] = CornerDistE(Ai, V, dNN);
        }
        else if (Ai[0] <= L_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i])
        {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {L_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbl[i];
            d_minus = Dtl[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle = false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle,
                                          L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i],
                                          intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle,
                                           L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i],
                                           intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            lb_record[i] = MIN(delta_plus, delta_minus);
        }
        else if (Ai[0] >= U_mu_Q[i] && L_sigma_Q[i] <= Ai[1] && Ai[1] <= U_sigma_Q[i])
        {
            V_plus = {U_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dbr[i];
            d_minus = Dtr[i];
            MiddleLeftOrRight = true;
            TopOrBottomMiddle = false;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle,
                                          L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i],
                                          intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle,
                                           L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i],
                                           intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            lb_record[i] = MIN(delta_plus, delta_minus);
        }
        else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] <= L_sigma_Q[i])
        {
            V_plus = {L_mu_Q[i], L_sigma_Q[i]};
            V_minus = {U_mu_Q[i], L_sigma_Q[i]};
            d_plus = Dbl[i];
            d_minus = Dbr[i];
            MiddleLeftOrRight = false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle,
                                          L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i],
                                          intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle,
                                           L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i],
                                           intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            lb_record[i] = MIN(delta_plus, delta_minus);
        }
        else if (L_mu_Q[i] <= Ai[0] && Ai[0] <= U_mu_Q[i] && Ai[1] >= U_sigma_Q[i])
        {
            V_plus = {L_mu_Q[i], U_sigma_Q[i]};
            V_minus = {U_mu_Q[i], U_sigma_Q[i]};
            d_plus = Dtl[i];
            d_minus = Dtr[i];
            MiddleLeftOrRight = false;
            TopOrBottomMiddle = true;
            double delta_plus = EdgeDistE(Ai, V_plus, d_plus, MiddleLeftOrRight, TopOrBottomMiddle,
                                          L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i],
                                          intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            double delta_minus = EdgeDistE(Ai, V_minus, d_minus, MiddleLeftOrRight, TopOrBottomMiddle,
                                           L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i], T[i],
                                           intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
            lb_record[i] = MIN(delta_plus, delta_minus);
        }
        else
        {
            lb_record[i] = CenterDist2D(mu_A[i], sigma_A[i], L_mu_Q[i], U_mu_Q[i], L_sigma_Q[i], U_sigma_Q[i],
                                        T[i], intervalLength_mu[i], intervalLength_sigma[i], lambda, num[i]);
        }
        lb += lb_record[i];
        if (mc3_lb_exceeds(lb, K, best_dist))
        {
            return mc3_finalize(lb, K);
        }
    }

    mc3_scale_accumulator(lb, lb_record, K);

    vector<double> max_val(m);
    upper_lemire(lb_record, m, r, max_val);
    vector<double> mv_right(m, 0);
    for (int k = 0; k < K; k++)
    {
        lower_upper_lemire(A_T[k], m, r, A_L, A_U);
        for (int i = 0; i < m; i++)
        {
            if (Q_T[k][i] < A_L[i])
            {
                mv_right[i] += dist(Q_T[k][i], A_L[i]);
            }
            else if (Q_T[k][i] > A_U[i])
            {
                mv_right[i] += dist(Q_T[k][i], A_U[i]);
            }
        }
    }

    for (int i = 0; i < m; i++)
    {
        lb += MAX(0.0, mv_right[i] - max_val[i]);
        if (mc3_scaled_lb_exceeds(lb, best_dist))
        {
            return sqrt(lb);
        }
    }
    return sqrt(lb);
}
