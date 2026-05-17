#include "funcs.h"

void lower_upper_lemire(vector<double> &Q, int m, int r, vector<double> &lowerBound, vector<double> &upperBound)
{
    int k = 2 * r + 1;
    int l_index = 0;
    int u_index = 0;
    std::vector<int> q_u(m);
    std::vector<int> q_l(m);
    int head_l = 0, tail_l = -1;
    int head_u = 0, tail_u = -1;
    // lowerBound.resize(m);
    // upperBound.resize(m);

    for (int i = 0; i < k - 1 - r; i++)
    {
        while (head_l <= tail_l && Q[q_l[tail_l]] >= Q[i])
            tail_l--;
        q_l[++tail_l] = i;
        while (head_u <= tail_u && Q[q_u[tail_u]] <= Q[i])
            tail_u--;
        q_u[++tail_u] = i;
    }

    for (int i = k - 1 - r; i < m; i++)
    {
        while (head_l <= tail_l && Q[q_l[tail_l]] >= Q[i])
            tail_l--;
        q_l[++tail_l] = i;
        while (q_l[head_l] <= i - k)
            head_l++;
        lowerBound[l_index++] = Q[q_l[head_l]];

        while (head_u <= tail_u && Q[q_u[tail_u]] <= Q[i])
            tail_u--;
        q_u[++tail_u] = i;
        while (q_u[head_u] <= i - k)
            head_u++;
        upperBound[u_index++] = Q[q_u[head_u]];
    }

    for (int i = m; i < m + r; i++)
    {
        while (q_u[head_u] <= i - k)
            head_u++;
        upperBound[u_index++] = Q[q_u[head_u]];
        while (q_l[head_l] <= i - k)
            head_l++;
        lowerBound[l_index++] = Q[q_l[head_l]];
    }
}

void lower_upper_lemire_index(vector<double> &Q, int m, int r, vector<pair<int, double>> &lowerBound, vector<pair<int, double>> &upperBound)
{
    int k = 2 * r + 1;
    int l_index = 0;
    int u_index = 0;
    vector<int> q_u(m);
    vector<int> q_l(m);
    int head_l = 0, tail_l = -1;
    int head_u = 0, tail_u = -1;
    // lowerBound.resize(m);
    // upperBound.resize(m);

    for (int i = 0; i < k - 1 - r; i++)
    {
        while (head_l <= tail_l && Q[q_l[tail_l]] >= Q[i])
            tail_l--;
        q_l[++tail_l] = i;
        while (head_u <= tail_u && Q[q_u[tail_u]] <= Q[i])
            tail_u--;
        q_u[++tail_u] = i;
    }

    for (int i = k - 1 - r; i < m; i++)
    {
        while (head_l <= tail_l && Q[q_l[tail_l]] >= Q[i])
            tail_l--;
        q_l[++tail_l] = i;
        while (q_l[head_l] <= i - k)
            head_l++;
        // lowerBound[l_index].first = q_l[head_l];
        // lowerBound[l_index++].second = Q[q_l[head_l]];
        lowerBound[l_index++] = {q_l[head_l], Q[q_l[head_l]]};

        while (head_u <= tail_u && Q[q_u[tail_u]] <= Q[i])
            tail_u--;
        q_u[++tail_u] = i;
        while (q_u[head_u] <= i - k)
            head_u++;
        // upperBound[u_index].first = q_u[head_u];
        // upperBound[u_index++].second = Q[q_u[head_u]];
        upperBound[u_index++] = {q_u[head_u], Q[q_u[head_u]]};
    }

    for (int i = m; i < m + r; i++)
    {
        while (q_u[head_u] <= i - k)
            head_u++;
        // upperBound[u_index].first = q_u[head_u];
        // upperBound[u_index++].second = Q[q_u[head_u]];
        upperBound[u_index++] = {q_u[head_u], Q[q_u[head_u]]};
        while (q_l[head_l] <= i - k)
            head_l++;
        // lowerBound[l_index].first = q_l[head_l];
        // lowerBound[l_index++].second = Q[q_l[head_l]];
        lowerBound[l_index++] = {q_l[head_l], Q[q_l[head_l]]};
    }
}

void upper_lemire(vector<double> &Q, int m, int r, vector<double> &upperBound)
{
    int k = 2 * r + 1;
    int u_index = 0;
    std::vector<int> q_u(m);
    int head_u = 0, tail_u = -1;
    // upperBound.resize(m);

    for (int i = 0; i < k - 1 - r; i++)
    {
        while (head_u <= tail_u && Q[q_u[tail_u]] <= Q[i])
            tail_u--;
        q_u[++tail_u] = i;
    }

    for (int i = k - 1 - r; i < m; i++)
    {
        while (head_u <= tail_u && Q[q_u[tail_u]] <= Q[i])
            tail_u--;
        q_u[++tail_u] = i;
        while (q_u[head_u] <= i - k)
            head_u++;
        upperBound[u_index++] = Q[q_u[head_u]];
    }

    for (int i = m; i < m + r; i++)
    {
        while (q_u[head_u] <= i - k)
            head_u++;
        upperBound[u_index++] = Q[q_u[head_u]];
    }
}
//
void lower_lemire(vector<double> &Q, int m, int r, vector<double> &lowerBound)
{
    int k = 2 * r + 1;
    int l_index = 0;
    std::vector<int> q_l(m);
    int head_l = 0, tail_l = -1;
    // lowerBound.resize(m);

    for (int i = 0; i < k - 1 - r; i++)
    {
        while (head_l <= tail_l && Q[q_l[tail_l]] >= Q[i])
            tail_l--;
        q_l[++tail_l] = i;
    }

    for (int i = k - 1 - r; i < m; i++)
    {
        while (head_l <= tail_l && Q[q_l[tail_l]] >= Q[i])
            tail_l--;
        q_l[++tail_l] = i;
        while (q_l[head_l] <= i - k)
            head_l++;
        lowerBound[l_index++] = Q[q_l[head_l]];
    }

    for (int i = m; i < m + r; i++)
    {
        while (q_l[head_l] <= i - k)
            head_l++;
        lowerBound[l_index++] = Q[q_l[head_l]];
    }
}