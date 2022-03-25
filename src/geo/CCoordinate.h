//
// Created by indro on 10/03/21.
//

#ifndef HERTZ_CCOORDINATE_H
#define HERTZ_CCOORDINATE_H

#include <Eigen/Dense>

namespace hertz::geo {
    template<typename Real>
    class CCoordinate : public Eigen::Matrix<Real, 2, 1> {
    public:
        typedef Eigen::Matrix<Real, 2, 1> Base;

        inline explicit CCoordinate(Real val = (Real) 0) { this->setConstant(val); }

        CCoordinate(Real x, Real y) : Base(x, y) {}

        template<typename Derived>
        inline explicit CCoordinate(const Eigen::MatrixBase<Derived> &other) : Base(other) {}

        template<typename Derived>
        CCoordinate &operator=(const Eigen::MatrixBase<Derived> &other) {
            this->Base::operator=(other);
            return *this;
        }

    };

}

#endif //HERTZ_CCOORDINATE_H
