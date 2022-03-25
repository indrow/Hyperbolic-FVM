//
// Created by indro on 10/03/21.
//

#ifndef HERTZ_CLINE_H
#define HERTZ_CLINE_H

#include "CCoordinate.h"
#include "CVector.h"

namespace hertz::geo {
    template<typename Real>
    class CLine {
    private:
        CCoordinate<Real> _startPoint;
        CCoordinate<Real> _endPoint;
        CCoordinate<Real> _dummyPoint;
        CVector<Real> _vector;
        Real _length;

    public:
        CLine();

        void setStartPoint(const CCoordinate<Real> &aPoint);

        void setEndPoint(const CCoordinate<Real> &aPoint);

        void setPoints(const CCoordinate<Real> &aStartPoint, const CCoordinate<Real> &anEndPoint);

        const CCoordinate<Real> &startPoint() const;

        const CCoordinate<Real> &endPoint() const;

        const CVector<Real> &vector();

        const Real &length();

        CCoordinate<Real> midPoint(Real ratio);
    };

    /*============================ Implementation =======================*/

    template<typename Real>
    CLine<Real>::CLine() {
        _startPoint(0);
        _endPoint(0);
        _dummyPoint(0);
        _vector(0);
        _length = (Real) 0;
    }

    template<typename Real>
    void CLine<Real>::setStartPoint(const CCoordinate<Real> &aPoint) {
        _startPoint = aPoint;
    }

    template<typename Real>
    void CLine<Real>::setEndPoint(const CCoordinate<Real> &aPoint) {
        _endPoint = aPoint;
    }

    template<typename Real>
    void CLine<Real>::setPoints(const CCoordinate<Real> &aStartPoint, const CCoordinate<Real> &anEndPoint) {
        _startPoint = aStartPoint;
        _endPoint = anEndPoint;
    }

    template<typename Real>
    const CCoordinate<Real> &CLine<Real>::startPoint() const {
        return _startPoint;
    }

    template<typename Real>
    const CCoordinate<Real> &CLine<Real>::endPoint() const {
        return _endPoint;
    }

    template<typename Real>
    const CVector<Real> &CLine<Real>::vector() {
        _vector = _endPoint - _startPoint;
        return _vector;
    }

    template<typename Real>
    const Real &CLine<Real>::length() {
        _length = vector().norm();
        return _length;
    }

    template<typename Real>
    CCoordinate<Real> CLine<Real>::midPoint(Real ratio) {
        _dummyPoint = _startPoint + vector() * ratio;
        return _dummyPoint;
    }
}

#endif //HERTZ_CLINE_H
