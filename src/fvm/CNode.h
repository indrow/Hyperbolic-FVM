//
// Created by indro on 10/03/21.
//

#ifndef HERTZ_CNODE_H
#define HERTZ_CNODE_H

namespace hertz::fvm {
    using hertz::geo::CCoordinate;

    template<typename Real>
    class CNode : public hertz::geo::CCoordinate<Real> {
    public:
        typedef CCoordinate <Real> Base;

        inline explicit CNode(Real value = (Real) 0) { this->setConstant(value); }

        inline CNode(const Real x, const Real y, const Real z) : Base(x, y, z) {}

        template<typename Derived>
        inline explicit CNode(const Eigen::MatrixBase<Derived> &p) : Base(p) {}

        template<typename Derived>
        CNode &operator=(const Eigen::MatrixBase<Derived> &p) {
            this->Base::operator=(p);
            return *this;
        }
    };
}

#endif //HERTZ_CNODE_H
