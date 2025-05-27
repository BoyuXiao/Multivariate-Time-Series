#ifndef CELL_H
#define CELL_H
struct Cell {
    std::pair<double, double> Lleft_top;
    std::pair<double, double> Lleft_bottom;
    std::pair<double, double> Lright_top;
    std::pair<double, double> Lright_bottom;

    std::vector<tuple<double, double,int>> dataList;

    double minDistToLleft;
    double minDistToLright;
    double minDistToLbottom;
    double minDistToLtop;

    double left_nearest_mu;
    double left_nearest_sigma;

    double right_nearest_mu;
    double right_nearest_sigma;

    double bottom_nearest_mu;
    double bottom_nearest_sigma;

    double top_nearest_mu;
    double top_nearest_sigma;

    double left_bottom_nearest_mu;
    double left_bottom_nearest_sigma;
    double left_bottom_nearest_dist;
    double left_top_nearest_mu;
    double left_top_nearest_sigma;
    double left_top_nearest_dist;

    double right_bottom_nearest_mu;
    double right_bottom_nearest_sigma;
    double right_bottom_nearest_dist;
    double right_top_nearest_mu;
    double right_top_nearest_sigma;
    double right_top_nearest_dist;

    int numIntervals_mu;
    int numIntervals_sigma;

    int left_t;
    int right_t;
    int bottom_t;
    int top_t;

    int left_bottom_t;
    int right_bottom_t;
    int right_top_t;
    int left_top_t;


    Cell(double muLeft, double muRight, double sigmaBottom, double sigmaTop)
        : minDistToLleft(MAX_VAL), minDistToLright(MAX_VAL), minDistToLbottom(MAX_VAL), minDistToLtop(MAX_VAL),
        left_nearest_mu(MAX_VAL),left_nearest_sigma(MAX_VAL),right_nearest_mu(MAX_VAL),right_nearest_sigma(MAX_VAL),
        bottom_nearest_mu(MAX_VAL),bottom_nearest_sigma(MAX_VAL),top_nearest_mu(MAX_VAL),top_nearest_sigma(MAX_VAL),
        left_bottom_nearest_mu(MAX_VAL),left_bottom_nearest_sigma(MAX_VAL),left_top_nearest_mu(MAX_VAL),left_top_nearest_sigma(MAX_VAL),
        right_bottom_nearest_mu(MAX_VAL),right_bottom_nearest_sigma(MAX_VAL),right_top_nearest_mu(MAX_VAL),right_top_nearest_sigma(MAX_VAL)
        {
        Lleft_top = {muLeft, sigmaTop};
        Lleft_bottom = {muLeft, sigmaBottom};
        Lright_top = {muRight, sigmaTop};
        Lright_bottom = {muRight, sigmaBottom};
    }


};
#endif //CELL_H
