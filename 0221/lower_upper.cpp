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