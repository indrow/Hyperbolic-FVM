//
// Created by indro on 10/03/21.
//

#ifndef HERTZ_DFACE_H
#define HERTZ_DFACE_H

#include <Eigen/Dense>
#include <vector>

namespace hertz::fvm {
    template<typename Real, typename Derived>
    class DFace {
    public:
        typedef Eigen::Matrix<Real, Eigen::Dynamic, 1> VectorXr;

        DFace();

        bool calculated;

        Derived maxWaveSpeed, flux;

        std::vector<VectorXr> quadrature;
    };

    template<typename Real, typename Derived>
    DFace<Real, Derived>::DFace() : calculated(false) {

    }
}

#endif //HERTZ_DFACE_H
