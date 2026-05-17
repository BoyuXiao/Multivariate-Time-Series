#include <vector>
#include "funcs.h"
#include "lb_early.h"
using namespace std;

namespace {

void fill_train_mu_sigma(const vector<vector<double>>& A, int m, vector<double>& mu_A,
                         vector<double>& sigma_A) {
    mu_A.resize(m);
    sigma_A.resize(m);
    mu_sigma(const_cast<vector<vector<double>>&>(A), mu_A, sigma_A);
}

/** LB_MC2 = sqrt(K * (LB_Keogh(mu)^2 + LB_Keogh(sigma)^2))；LB_Keogh_2 已返回平方和。 */
double mc2_inner_from_keogh_sq(double keogh_mu_sq, double keogh_sigma_sq, int K) {
    return static_cast<double>(K) * (keogh_mu_sq + keogh_sigma_sq);
}

}  // namespace

double LB_MC2_early_cached(const vector<double>& mu_A,
                          const vector<double>& sigma_A,
                          const vector<double>& mu_Q,
                          const vector<double>& sigma_Q,
                          const vector<double>& mu_q_lower,
                          const vector<double>& mu_q_upper,
                          const vector<double>& sigma_q_lower,
                          const vector<double>& sigma_q_upper,
                          int K,
                          int m,
                          double best_dist) {
    const double thresh_sq = lb_threshold_sq(best_dist);
    const double result = mc2_inner_from_keogh_sq(
        LB_Keogh_2_env_early(const_cast<vector<double>&>(mu_A), mu_q_lower, mu_q_upper, m, best_dist),
        LB_Keogh_2_env_early(const_cast<vector<double>&>(sigma_A), sigma_q_lower, sigma_q_upper, m,
                             best_dist),
        K);
    if (lb_can_early_stop(best_dist) && result >= thresh_sq) {
        return sqrt(result);
    }
    return sqrt(result);
}

double LB_MC2_early(vector<vector<double>>& Q, vector<vector<double>>& A, vector<double>& mu_Q,
                    vector<double>& sigma_Q, int K, int m, int r, double best_dist) {
    vector<double> mu_A(m);
    vector<double> sigma_A(m);
    fill_train_mu_sigma(A, m, mu_A, sigma_A);
    const double thresh_sq = lb_threshold_sq(best_dist);
    const double result = mc2_inner_from_keogh_sq(
        LB_Keogh_2_early(mu_A, mu_Q, m, r, best_dist),
        LB_Keogh_2_early(sigma_A, sigma_Q, m, r, best_dist), K);
    if (lb_can_early_stop(best_dist) && result >= thresh_sq) {
        return sqrt(result);
    }
    return sqrt(result);
}

double LB_MC2(vector<vector<double>>& Q, vector<vector<double>>& A, vector<double>& mu_Q,
              vector<double>& sigma_Q, int K, int m, int r) {
    vector<double> mu_A(m);
    vector<double> sigma_A(m);
    fill_train_mu_sigma(A, m, mu_A, sigma_A);
    const double result = mc2_inner_from_keogh_sq(LB_Keogh_2(mu_A, mu_Q, m, r),
                                                 LB_Keogh_2(sigma_A, sigma_Q, m, r), K);
    return sqrt(result);
}
