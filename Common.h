#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <map>
#include <omp.h>
#include <vector>
#include <tuple>
#include <thread>
#include <unordered_map>
#include <utility>
#include "Equipment.h"

using cf = std::complex<DeviceArgType>;

constexpr auto PI = 3.14159265358979323846;
constexpr auto epsilon = 1e-6;
constexpr auto Imaginer = std::complex<float>{0, 1};
const cf alpha = { cosf(PI * 2 / 3), sinf(PI * 2 / 3) };
const auto maxThreadsNum = std::thread::hardware_concurrency();
const Eigen::MatrixXcf St = (Eigen::MatrixXcf(3, 3) << 1, 1, 1, alpha * alpha, alpha, 1, alpha, alpha * alpha, 1).finished();//�������(1 / 3 * S)����

using ThreeSequenceData = struct ThreeSequenceData
{
    Eigen::MatrixXf lineData1;
    Eigen::MatrixXf lineData2;
    Eigen::MatrixXf lineData0;
};
