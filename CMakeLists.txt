cmake_minimum_required(VERSION 3.20)
project(TLB)

set(CMAKE_CXX_STANDARD 20)
if(MSVC)
    add_compile_options(/utf-8)
endif()
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")

# 下界 / DTW 核心（不含 1NN 与 TLB 比值逻辑）
set(LB_CORE_SOURCES
        dtw.cpp
        dist.cpp
        LB_MV.cpp
        LB_MC1.cpp
        LB_Keogh.cpp
        LB_MC2.cpp
        LB_MC_NAIVE1.cpp
        LB_MC_NAIVE2.cpp
        mu_sigma.cpp
        read_data.cpp
        lower_upper.cpp
        only_mu.cpp
        transpose.cpp
        LB_MC3E.cpp
        split_two_dim.cpp
        LB_TI.cpp
        LB_PC.cpp
        LB_Petitjean.cpp
)

# TLB.exe：1NN + 剪枝加速 benchmark（原 main）
add_executable(TLB main.cpp TLB.cpp nn1_lb.cpp ${LB_CORE_SOURCES})

# TLB_ratio.exe：经典 mean(LB/DTW)，与 1NN 解耦
add_executable(TLB_ratio tlb_main.cpp tlb_ratio.cpp ${LB_CORE_SOURCES})

# 单元 / 冒烟测试
add_executable(TLB_test tests/test_tlb.cpp TLB.cpp nn1_lb.cpp ${LB_CORE_SOURCES})

target_include_directories(TLB PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(TLB_ratio PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(TLB_test PRIVATE ${CMAKE_SOURCE_DIR})
