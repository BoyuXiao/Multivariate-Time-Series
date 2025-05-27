#ifndef CELL_H
#define CELL_H
struct Cell {
    // 线段Lleft的上下顶点坐标
    std::pair<double, double> Lleft_top;
    std::pair<double, double> Lleft_bottom;
    // 线段Lright的上下顶点坐标
    std::pair<double, double> Lright_top;
    std::pair<double, double> Lright_bottom;
    // 线段Lbottom的左右顶点坐标
    std::pair<double, double> Lbottom_left;
    std::pair<double, double> Lbottom_right;
    // 线段Ltop的左右顶点坐标
    std::pair<double, double> Ltop_left;
    std::pair<double, double> Ltop_right;

    std::vector<pair<double, double>> dataList;

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



    Cell(double muLeft, double muRight, double sigmaBottom, double sigmaTop)
        : minDistToLleft(MAX_VAL), minDistToLright(MAX_VAL), minDistToLbottom(MAX_VAL), minDistToLtop(MAX_VAL),
        left_nearest_mu(MAX_VAL),left_nearest_sigma(MAX_VAL),right_nearest_mu(MAX_VAL),right_nearest_sigma(MAX_VAL),
        bottom_nearest_mu(MAX_VAL),bottom_nearest_sigma(MAX_VAL),top_nearest_mu(MAX_VAL),top_nearest_sigma(MAX_VAL)
        {
        Lleft_top = {muLeft, sigmaTop};
        Lleft_bottom = {muLeft, sigmaBottom};
        Lright_top = {muRight, sigmaTop};
        Lright_bottom = {muRight, sigmaBottom};
        Lbottom_left = {muLeft, sigmaBottom};
        Lbottom_right = {muRight, sigmaBottom};
        Ltop_left = {muLeft, sigmaTop};
        Ltop_right = {muRight, sigmaTop};
    }
};
#endif //CELL_H
