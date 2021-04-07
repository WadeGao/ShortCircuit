#define EIGEN_USE_MKL_ALL
#include <iostream>
#include "Common.h"
#include "Grid.h"

int main()
{
    
    Eigen::MatrixXf mat(5, 6);
    mat <<
        1,  2,  0.000, 0.105, 0.000,  1.05,
        2,  3,  0.024, 0.065, 0.032,  1.00,
        2,  4,  0.030, 0.080, 0.040,  1.00,
        4,  3,  0.018, 0.050, 0.026,  1.00,
        5,  4,  0.000, 0.184, 0.000,  0.96;
    Grid my_grid({ mat, mat, mat }, { {5,4,0.96},{1,2,1.05} });

    std::cout << my_grid.getAdmittanceMatrix(SEQUENCE::POSITIVE) << std::endl;
    return 0;
}
