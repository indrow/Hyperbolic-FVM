//
// Created by indro on 10/03/21.
//

#ifndef HERTZ_PHYSICS_SCALAR_H
#define HERTZ_PHYSICS_SCALAR_H

#include "common/CPi.h"
#include "geo/CCoordinate.h"
#include "geo/CVector.h"

namespace hertz::physics {
    using hertz::geo::CCoordinate;
    using hertz::geo::CVector;

    template<typename Real>
    class Scalar {
    public:
        static Real linearAdvection(const Real &u);

        static Real Burgers(const Real &u);

        static Real dLinearAdvection(const Real &u);

        static Real dBurgers(const Real &u);

        static Real sine(const CCoordinate<Real> &coordinate);

        static Real Rusanov(const Real &ul, const Real &ur, const Real &a_norm, Real (*fn)(const Real &),
                            Real (*gn)(const Real &), const CVector<Real> &norm);
    };

    template<typename Real>
    Real Scalar<Real>::linearAdvection(const Real &u) {
        return u;
    }

    template<typename Real>
    Real Scalar<Real>::dLinearAdvection(const Real &u) {
        return 1.0;
    }

    template<typename Real>
    Real Scalar<Real>::sine(const CCoordinate<Real> &coordinate) {
        return std::sin(CPi<Real>() * (coordinate.x()));
    }

    template<typename Real>
    Real Scalar<Real>::Rusanov(const Real &ul, const Real &ur, const Real &a_norm, Real (*fn)(const Real &),
                               Real (*gn)(const Real &), const CVector<Real> &norm) {

        return 0.5 * (CVector<Real>(fn(ul) + fn(ur), gn(ul) + gn(ur)).derived().template dot(norm) -
                a_norm * (ur - ul));
    }

    template<typename Real>
    Real Scalar<Real>::Burgers(const Real &u) {
        return u * u;
    }

    template<typename Real>
    Real Scalar<Real>::dBurgers(const Real &u) {
        return 2.0 * u;
    }
}

#endif //HERTZ_PHYSICS_SCALAR_H
