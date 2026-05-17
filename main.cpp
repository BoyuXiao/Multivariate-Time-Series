#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include "funcs.h"
#include "nn1_lb.h"
#include "def.h"
using namespace std;
namespace fs = std::filesystem;

static string make_result_filename() {
    const auto now = chrono::system_clock::now();
    const time_t tt = chrono::system_clock::to_time_t(now);
    tm tm_buf{};
#ifdef _WIN32
    localtime_s(&tm_buf, &tt);
#else
    localtime_r(&tt, &tm_buf);
#endif
    ostringstream oss;
    oss << "nn1_" << put_time(&tm_buf, "%Y%m%d_%H%M%S") << ".txt";
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

static string init_result_file(Nn1Config& config) {
    if (!config.result_file.empty()) {
        return config.result_file;
    }
    const fs::path res_dir = detect_res_dir();
    fs::create_directories(res_dir);
    config.result_file = (res_dir / make_result_filename()).string();
    return config.result_file;
}

// 1NN + 下界剪枝 benchmark（输出 res/nn1_*.txt）。经典 TLB 比值请用 TLB_ratio.exe。
// Table 2: 21 multidimensional UEA datasets
int main(int argc, char* argv[]) {
    vector<string> testFilePaths = {
        data1_test, data2_test, data3_test, data4_test, data5_test,
        data7_test, data8_test, data9_test, data11_test, data12_test,
        data13_test, data14_test, data15_test, data16_test,
        data18_test, data20_test, data22_test, data23_test, data24_test, data25_test, data26_test,
    };

    vector<string> trainFilePaths = {
        data1_train, data2_train, data3_train, data4_train, data5_train,
        data7_train, data8_train, data9_train, data11_train, data12_train,
        data13_train, data14_train, data15_train, data16_train,
        data18_train, data20_train, data22_train, data23_train, data24_train, data25_train, data26_train,
    };

    vector<string> dataNames = {
        data1_name, data2_name, data3_name, data4_name, data5_name,
        data7_name, data8_name, data9_name, data11_name, data12_name,
        data13_name, data14_name, data15_name, data16_name,
        data18_name, data20_name, data22_name, data23_name, data24_name, data25_name, data26_name,
    };

    if (argc < 3) {
        cerr << "Usage: " << argv[0]
             << " start_index end_index [r_ratio] [early_stop:0|1] [result_file] [methods]" << endl;
        cerr << "  Defaults: r_ratio=0.01, early_stop=1" << endl;
        cerr << "  methods: comma-separated, e.g. MC1,MC2,MC3 (default: all)" << endl;
        cerr << "  Example: " << argv[0] << " 0 0" << endl;
        cerr << "  Example: " << argv[0] << " 0 20 0.01 1 res/nn1_patch.txt MC1,MC2,MC3" << endl;
        return 1;
    }

    int start_index = atoi(argv[1]);
    int end_index = atoi(argv[2]);

    Nn1Config config;
    if (argc >= 4) {
        config.r_ratio = atof(argv[3]);
    }
    if (argc >= 5) {
        config.use_early_stopping = (atoi(argv[4]) != 0);
    }
    if (argc >= 6) {
        config.result_file = argv[5];
    }
    if (argc >= 7) {
        config.methods_filter.push_back(argv[6]);
    }

    if (start_index < 0 || start_index >= static_cast<int>(testFilePaths.size()) ||
        end_index < 0 || end_index >= static_cast<int>(testFilePaths.size()) ||
        start_index > end_index) {
        cerr << "Invalid start or end index (valid range: 0-" << testFilePaths.size() - 1 << ")." << endl;
        return 1;
    }

    cout << "r_ratio=" << config.r_ratio << " early_stop=" << (config.use_early_stopping ? "on" : "off");
    if (!config.methods_filter.empty()) {
        cout << " methods=" << config.methods_filter[0];
    }
    cout << endl;
    const string result_path = init_result_file(config);
    cout << "Results -> " << result_path << endl;

    for (int i = start_index; i <= end_index; i++) {
        vector<vector<vector<double>>> data_test = readCSV(testFilePaths[i]);
        vector<vector<vector<double>>> data_train = readCSV(trainFilePaths[i]);
        TLB(data_test, data_train, dataNames[i], config);
    }
    return 0;
}
