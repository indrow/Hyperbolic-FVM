//
// Created by indro on 10/03/21.
//

#ifndef HERTZ_DCELL_H
#define HERTZ_DCELL_H

#include <Eigen/Dense>
#include <vector>

namespace hertz::fvm {
    template<typename Real, typename Derived>
    class DCell {
    public:
        typedef Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> RealMatrix;
        typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> RealVector;

        DCell();

        bool calculate;

        Derived u, maxWaveSpeed, residual;

        std::vector<Derived> urk;

        RealMatrix AP3;

        RealVector CP3;
    };

    template<typename Real, typename Derived>
    DCell<Real, Derived>::DCell() : calculate(false) {
        urk = std::vector<Derived>(3);
        AP3 = RealMatrix(10, 10);
        CP3 = RealVector(10);
    }
}

#endif //HERTZ_DCELL_H
