#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <float.h>
#include <map>
#include <omp.h>
#include <tuple>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>

//constexpr auto M_PI = 3.14159265358979323846;
constexpr auto epsilon = 1e-6;
const auto Imaginer = std::complex<float>{ 0, 1 };
const std::complex<float> alpha = { cosf(M_PI * 2 / 3), sinf(M_PI * 2 / 3) };
const auto maxThreadsNum = std::thread::hardware_concurrency();
const Eigen::MatrixXcf St = (Eigen::MatrixXcf(3, 3) << 1, 1, 1, alpha * alpha, alpha, 1, alpha, alpha * alpha, 1).finished();//ÆäÄæ¾ØÕó¼´(1 / 3 * S)¾ØÕó
