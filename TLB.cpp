#include "funcs.h"
#include "nn1_lb.h"

/** 1NN + 下界剪枝 benchmark（由 main.cpp / TLB.exe 调用）。经典 TLB 比值见 TLB_ratio.exe。 */
void TLB(std::vector<std::vector<std::vector<double>>>& test_set,
         std::vector<std::vector<std::vector<double>>>& train_set,
         const std::string& name,
         const Nn1Config& config) {
    run_nn1_lb_test(test_set, train_set, name, config);
}
