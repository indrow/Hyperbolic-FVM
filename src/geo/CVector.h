//
// Created by indro on 10/03/21.
//

#ifndef HERTZ_CVECTOR_H
#define HERTZ_CVECTOR_H

namespace hertz::geo {
    template<typename Real>
    class CVector : public Eigen::Matrix<Real, 2, 1> {
    public:
        typedef Eigen::Matrix<Real, 2, 1> Base;

        inline explicit CVector(Real val = (Real) 0) { this->setConstant(val); }

        CVector(Real x, Real y) : Base(x, y) {}

        template<typename Derived>
        inline explicit CVector(const Eigen::MatrixBase<Derived> &other) : Base(other) {}

        template<typename Derived>
        CVector &operator=(const Eigen::MatrixBase<Derived> &other) {
            this->Base::operator=(other);
            return *this;
        }

    };
}

#endif //HERTZ_CVECTOR_H
