#include <vector>
#include "funcs.h"
using namespace std;

void mu_sigma(vector<vector<double>>& A, int K, int m, vector<double>& mu_A, vector<double>& sigma_A) {
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






