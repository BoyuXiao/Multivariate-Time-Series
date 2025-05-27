#include <vector>
#include "funcs.h"
using namespace std;

void only_mu(vector<vector<double>>& A, int K, int m, vector<double>& mu_A) {
    for (int i = 0; i < m; ++i) {
        double sum_A = 0.0;
        for (int k = 0; k < K; ++k) {
            sum_A += A[i][k];
        }
        mu_A[i] = sum_A / K;
    }
}







