#include "funcs.h"

void lower_upper_lemire(vector<double>& Q, int m, int r, vector<double>& lowerBound, vector<double>& upperBound)
{
    int k = 2 * r + 1;
    std::vector<int> q_u(m);
    std::vector<int> q_l(m);
    int head_l = 0, tail_l = -1;
    int head_u = 0, tail_u = -1;

    for (int i = 0; i < k - 1 - r; i++) {
        while (head_l <= tail_l && Q[q_l[tail_l]] >= Q[i]) tail_l--;
        q_l[++tail_l] = i;
        while (head_u <= tail_u && Q[q_u[tail_u]] <= Q[i]) tail_u--;
        q_u[++tail_u] = i;
    }

    for (int i = k - 1 - r; i < m; i++) {
        while (head_l <= tail_l && Q[q_l[tail_l]] >= Q[i]) tail_l--;
        q_l[++tail_l] = i;
        while (q_l[head_l] <= i - k) head_l++;
        lowerBound.push_back(Q[q_l[head_l]]);

        while (head_u <= tail_u && Q[q_u[tail_u]] <= Q[i]) tail_u--;
        q_u[++tail_u] = i;
        while (q_u[head_u] <= i - k) head_u++;
        upperBound.push_back(Q[q_u[head_u]]);
    }

    for (int i = m; i < m + r; i++) {
        while (q_u[head_u] <= i - k) head_u++;
        upperBound.push_back(Q[q_u[head_u]]);
        while (q_l[head_l] <= i - k) head_l++;
        lowerBound.push_back(Q[q_l[head_l]]);
    }
}

void lower_upper_lemire_index(vector<double>& Q, int m, int r, vector<pair<int,double>>& lowerBound, vector<pair<int,double>>& upperBound)
{
    int k = 2 * r + 1;
    vector<int> q_u(m);
    vector<int> q_l(m);
    int head_l = 0, tail_l = -1;
    int head_u = 0, tail_u = -1;

    for (int i = 0; i < k - 1 - r; i++) {
        while (head_l <= tail_l && Q[q_l[tail_l]] >= Q[i]) tail_l--;
        q_l[++tail_l] = i;
        while (head_u <= tail_u && Q[q_u[tail_u]] <= Q[i]) tail_u--;
        q_u[++tail_u] = i;
    }

    for (int i = k - 1 - r; i < m; i++) {
        while (head_l <= tail_l && Q[q_l[tail_l]] >= Q[i]) tail_l--;
        q_l[++tail_l] = i;
        while (q_l[head_l] <= i - k) head_l++;
        lowerBound.emplace_back(q_l[head_l],Q[q_l[head_l]]);

        while (head_u <= tail_u && Q[q_u[tail_u]] <= Q[i]) tail_u--;
        q_u[++tail_u] = i;
        while (q_u[head_u] <= i - k) head_u++;
        upperBound.emplace_back(q_u[head_u],Q[q_u[head_u]]);
    }

    for (int i = m; i < m + r; i++) {
        while (q_u[head_u] <= i - k) head_u++;
        upperBound.emplace_back(q_u[head_u],Q[q_u[head_u]]);
        while (q_l[head_l] <= i - k) head_l++;
        lowerBound.emplace_back(q_l[head_l],Q[q_l[head_l]]);
    }
}


void lower_upper_lemire_pca(int m, vector<double>& mu_A, vector<double>& sigma_A,
                  double& U_mu_Q, double& L_mu_Q, double& U_sigma_Q, double& L_sigma_Q)
{
    double u_mu,l_mu,u_sigma,l_sigma;
    u_mu = mu_A[0];
    l_mu = mu_A[0];
    u_sigma = sigma_A[0];
    l_sigma = sigma_A[0];

    for (int i = 1; i < m; ++i) {
        if (mu_A[i] > u_mu) {
            u_mu = mu_A[i];
        }
        if (mu_A[i] < l_mu) {
            l_mu = mu_A[i];
        }
    }

    for (int i = 1; i < m; ++i) {
        if (sigma_A[i] > u_sigma) {
            u_sigma = sigma_A[i];
        }
        if (sigma_A[i] < l_sigma) {
            l_sigma = sigma_A[i];
        }
    }
    U_mu_Q=u_mu;
    L_mu_Q=l_mu;
    U_sigma_Q=u_sigma;
    L_sigma_Q=l_sigma;
}

void upper_lemire(vector<double>& Q, int m, int r, vector<double>& upperBound)
{
    int k = 2 * r + 1;
    std::vector<int> q_u(m);
    int head_u = 0, tail_u = -1;

    for (int i = 0; i < k - 1 - r; i++) {
        while (head_u <= tail_u && Q[q_u[tail_u]] <= Q[i]) tail_u--;
        q_u[++tail_u] = i;
    }

    for (int i = k - 1 - r; i < m; i++) {
        while (head_u <= tail_u && Q[q_u[tail_u]] <= Q[i]) tail_u--;
        q_u[++tail_u] = i;
        while (q_u[head_u] <= i - k) head_u++;
        upperBound.push_back(Q[q_u[head_u]]);
    }

    for (int i = m; i < m + r; i++) {
        while (q_u[head_u] <= i - k) head_u++;
        upperBound.push_back(Q[q_u[head_u]]);
    }
}
//
void lower_lemire(vector<double>& Q, int m, int r, vector<double>& lowerBound)
{
    int k = 2 * r + 1;
    std::vector<int> q_l(m);
    int head_l = 0, tail_l = -1;

    for (int i = 0; i < k - 1 - r; i++) {
        while (head_l <= tail_l && Q[q_l[tail_l]] >= Q[i]) tail_l--;
        q_l[++tail_l] = i;
    }

    for (int i = k - 1 - r; i < m; i++) {
        while (head_l <= tail_l && Q[q_l[tail_l]] >= Q[i]) tail_l--;
        q_l[++tail_l] = i;
        while (q_l[head_l] <= i - k) head_l++;
        lowerBound.push_back(Q[q_l[head_l]]);
    }

    for (int i = m; i < m + r; i++) {
        while (q_l[head_l] <= i - k) head_l++;
        lowerBound.push_back(Q[q_l[head_l]]);
    }
}