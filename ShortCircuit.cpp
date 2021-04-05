#define EIGEN_USE_MKL_ALL
#include <iostream>
#include "Common.h"

int main() {

	Eigen::MatrixXcf mat(3, 3);
	mat << -1, 2, -3, 4, -5, 6, 7, 8, -9;
	std::cout << mat.cwiseAbs() << std::endl;
	return 0;
}