#include <vector>
#include "funcs.h"
using namespace std;

void mu_sigma(vector<vector<double>>& A, vector<double>& mu_A, vector<double>& sigma_A) {
    int m = A.size();
    int K = A[0].size();
    mu_A.resize(m);
    sigma_A.resize(m);
    for (int i = 0; i < m; ++i) {
        double sum_A = 0.0;
        for (int k = 0; k < K; ++k) {
            sum_A += A[i][k];
        }
        mu_A[i] = sum_A / K;

        double variance_A = 0.0;
        for (int k = 0; k < K; ++k) {
            variance_A += (A[i][k] - mu_A[i])*(A[i][k] - mu_A[i]);
        }
        variance_A /= K;
        sigma_A[i] = sqrt(variance_A);
    }
}

void mu_sigma_single(vector<double>& A,double& mu_A, double& sigma_A) {
    int  K = A.size();
    double sum_A = 0.0;
    for (double val : A) {
        sum_A += val;
    }
    mu_A = sum_A / K;

    double variance_A = 0.0;
    for (double val : A) {
        variance_A += (val - mu_A) * (val - mu_A);
    }
    variance_A /= K;

    sigma_A = std::sqrt(variance_A);
}






