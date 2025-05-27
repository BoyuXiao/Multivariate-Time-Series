#include "funcs.h"

double dtw(vector<vector<double>>& A, vector<vector<double>>& B, int K, int m, int r)
{
    vector<double> cost(2 * r + 1, INF);
    vector<double> cost_prev(2 * r + 1, INF);
    vector<double> cost_tmp;
    int i,j,k;
    double x,y,z,min_cost;

    for (i=0; i<m; i++)
    {
        k = max(0,r-i);
        min_cost = INF;

        for(j=max(0,i-r); j<=min(m-1,i+r); j++, k++)
        {
            if ((i==0)&&(j==0))
            {
                cost[k]=euclideanDistance_2(A[0],B[0]);
                min_cost = cost[k];
                continue;
            }

            if ((j-1<0)||(k-1<0))     y = INF;
            else                      y = cost[k-1];
            if ((i-1<0)||(k+1>2*r))   x = INF;
            else                      x = cost_prev[k+1];
            if ((i-1<0)||(j-1<0))     z = INF;
            else                      z = cost_prev[k];

            cost[k] = min( min( x, y) , z) + euclideanDistance_2(A[i],B[j]);

            if (cost[k] < min_cost)
                min_cost = cost[k];
        }

        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--;

    double final_dtw = cost_prev[k];
    return sqrt(final_dtw);
}

