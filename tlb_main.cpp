/** 独立 TLB 指标程序：mean(LB/DTW)，与 1NN 剪枝 benchmark（TLB.exe / main.cpp）解耦。 */

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "def.h"
#include "funcs.h"
#include "tlb_ratio.h"

namespace fs = std::filesystem;

static std::string make_result_filename() {
    const auto now = std::chrono::system_clock::now();
    const time_t tt = std::chrono::system_clock::to_time_t(now);
    tm tm_buf{};
#ifdef _WIN32
    localtime_s(&tm_buf, &tt);
#else
    localtime_r(&tt, &tm_buf);
#endif
    std::ostringstream oss;
    oss << "tlb_" << std::put_time(&tm_buf, "%Y%m%d_%H%M%S") << ".txt";
    return oss.str();
}

static fs::path detect_res_dir() {
    if (fs::exists("CMakeLists.txt")) {
        return fs::path("res");
    }
    if (fs::exists("../../CMakeLists.txt")) {
        return fs::path("../../res");
    }
    return fs::path("res");
}

static std::string init_result_file(TlbRatioConfig& config) {
    if (!config.result_file.empty()) {
        return config.result_file;
    }
    const fs::path res_dir = detect_res_dir();
    fs::create_directories(res_dir);
    config.result_file = (res_dir / make_result_filename()).string();
    return config.result_file;
}

int main(int argc, char* argv[]) {
    const std::vector<std::string> test_paths = {
        data1_test,  data2_test,  data3_test,  data4_test,  data5_test,
        data7_test,  data8_test,  data9_test,  data11_test, data12_test,
        data13_test, data14_test, data15_test, data16_test,
        data18_test, data20_test, data22_test, data23_test, data24_test, data25_test, data26_test,
    };
    const std::vector<std::string> train_paths = {
        data1_train,  data2_train,  data3_train,  data4_train,  data5_train,
        data7_train,  data8_train,  data9_train,  data11_train, data12_train,
        data13_train, data14_train, data15_train, data16_train,
        data18_train, data20_train, data22_train, data23_train, data24_train, data25_train,
        data26_train,
    };
    const std::vector<std::string> names = {
        data1_name,  data2_name,  data3_name,  data4_name,  data5_name,
        data7_name,  data8_name,  data9_name,  data11_name, data12_name,
        data13_name, data14_name, data15_name, data16_name,
        data18_name, data20_name, data22_name, data23_name, data24_name, data25_name, data26_name,
    };

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " start_index end_index [r_ratio] [ti_p] [pc_q] [result_file] [methods]\n"
                  << "  Computes mean(LB/DTW) over all query-train pairs (full LB, no 1NN early stop).\n"
                  << "  methods: comma-separated, e.g. MC1 or MV,MC2 (default: all)\n"
                  << "  Defaults: r_ratio=0.01, ti_p=5, pc_q=2.0\n"
                  << "  Example: " << argv[0] << " 0 3 0.01 5 2.0 res/tlb_sweep.txt MC1\n";
        return 1;
    }

    const int start_index = std::atoi(argv[1]);
    const int end_index = std::atoi(argv[2]);

    TlbRatioConfig config;
    if (argc >= 4) {
        config.r_ratio = std::atof(argv[3]);
    }
    if (argc >= 5) {
        config.ti_p = std::atoi(argv[4]);
    }
    if (argc >= 6) {
        config.pc_q = std::atof(argv[5]);
    }
    if (argc >= 7) {
        config.result_file = argv[6];
    }
    if (argc >= 8) {
        config.methods_filter.push_back(argv[7]);
    }

    if (start_index < 0 || start_index >= static_cast<int>(test_paths.size()) || end_index < 0 ||
        end_index >= static_cast<int>(test_paths.size()) || start_index > end_index) {
        std::cerr << "Invalid index (0-" << test_paths.size() - 1 << ").\n";
        return 1;
    }

    const std::string result_path = init_result_file(config);
    std::cout << "TLB ratio test  r_ratio=" << config.r_ratio << " ti_p=" << config.ti_p
              << " pc_q=" << config.pc_q;
    if (!config.methods_filter.empty()) {
        std::cout << " methods=" << config.methods_filter[0];
    } else {
        std::cout << " methods=all";
    }
    std::cout << "\nResults -> " << result_path << std::endl;

    for (int i = start_index; i <= end_index; ++i) {
        auto test_set = readCSV(test_paths[i]);
        auto train_set = readCSV(train_paths[i]);
        run_tlb_ratio_test(test_set, train_set, names[i], config);
    }
    return 0;
}
