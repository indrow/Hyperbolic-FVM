//
// Created by indro on 10/03/21.
//

#ifndef HERTZ_CFACE_H
#define HERTZ_CFACE_H

#include <vector>

#include "common/CType.h"
#include "geo/CCoordinate.h"
#include "geo/CVector.h"

namespace hertz::fvm {
    using hertz::geo::CCoordinate;
    using hertz::geo::CVector;

    enum CFaceType {
        FT_NONE,
        FT_POINT,
        FT_LINE,
    };

    enum BoundaryType {
        BT_NONE,
        BT_INTERIOR,
        BT_PERIODIC,
    };

    class CFaceBoundaryType {
    private:
        BoundaryType _boundaryType;

    public:

        CFaceBoundaryType() : _boundaryType(BT_NONE) {}

        explicit CFaceBoundaryType(const BoundaryType boundaryType) : _boundaryType(boundaryType) {}

        CFaceBoundaryType &operator=(const BoundaryType &boundaryType) {
            this->_boundaryType = boundaryType;
            return *this;
        }

        BoundaryType operator()() const {
            return this->_boundaryType;
        }

        bool operator==(const BoundaryType &boundaryType) const {
            return this->_boundaryType == boundaryType;
        }

        bool operator!=(const BoundaryType &boundaryType) const {
            return this->_boundaryType != boundaryType;
        }
    };

    template<typename Real>
    class CFace {
    private:
        CFaceType _faceType;
        CFaceBoundaryType _boundaryType;

        std::vector<Integer> _nodeIds;
        Integer _nbNode;
        Integer _pairFaceId;
        Integer _cellId, _neighborId;
        bool _hasPair;

        Real _area;
        CVector<Real> _normal;
        CCoordinate<Real> _barycenter;
        std::vector<CCoordinate<Real>> _quadrature;
        Integer _nbQuadrature;

    public:
        CFace();

        void setFaceType(CFaceType aFaceType);

        [[nodiscard]] const CFaceType &faceType() const;

        void setFaceBoundaryType(BoundaryType aBoundaryType);

        [[nodiscard]] const CFaceBoundaryType &boundaryType() const;

        void addNodeId(Integer aNodeId);

        void replaceNode(Integer nodeIndex, Integer aNodeId);

        [[nodiscard]] const Integer &nodeId(const Integer &nodeIndex) const;

        [[nodiscard]] const std::vector<Integer> &nodeIds() const;

        [[nodiscard]] const Integer &nbNode() const;

        void setPairFaceId(Integer aPairFaceId);

        [[nodiscard]] const Integer &pairFaceId() const;

        [[nodiscard]] bool hasPair() const;

        void setCellId(Integer aCellId);

        [[nodiscard]] const Integer &cellId() const;

        void setNeighborId(Integer aNeighborId);

        [[nodiscard]] const Integer &neighborId() const;

        void setArea(const Real &anArea);

        [[nodiscard]] const Real &area() const;

        void setNormal(const CVector<Real> &aNormalVector);

        [[nodiscard]] const CVector<Real> &normal() const;

        void setBarycenter(const CCoordinate<Real> &aPoint);

        [[nodiscard]] const CCoordinate<Real> &barycenter() const;

        void addQuadraturePt(const CCoordinate<Real> &aPoint);

        [[nodiscard]] const CCoordinate<Real> &quadrature(Integer quadratureIndex) const;

        [[nodiscard]] const Integer &nbQuadrature();

        void reset();
    };

    /*============================ Implementation =======================*/

    template<typename Real>
    CFace<Real>::CFace() : _faceType(FT_NONE), _boundaryType(BT_NONE), _nbNode(0), _pairFaceId(0), _cellId(0),
                           _neighborId(0), _hasPair(false), _area(0), _normal(0), _nbQuadrature(0) {
    }

    template<typename Real>
    void CFace<Real>::setFaceType(CFaceType aFaceType) {
        _faceType = aFaceType;
    }

    template<typename Real>
    const CFaceType &CFace<Real>::faceType() const {
        return _faceType;
    }

    template<typename Real>
    void CFace<Real>::setFaceBoundaryType(BoundaryType aBoundaryType) {
        _boundaryType = aBoundaryType;
    }

    template<typename Real>
    const CFaceBoundaryType &CFace<Real>::boundaryType() const {
        return _boundaryType;
    }

    template<typename Real>
    void CFace<Real>::addNodeId(Integer aNodeId) {
        _nodeIds.push_back(aNodeId);
        _nbNode++;
    }

    template<typename Real>
    void CFace<Real>::replaceNode(Integer nodeIndex, Integer aNodeId) {
        _nodeIds.at(nodeIndex) = aNodeId;
    }

    template<typename Real>
    const Integer &CFace<Real>::nodeId(const Integer &nodeIndex) const {
        return _nodeIds.at(nodeIndex);
    }

    template<typename Real>
    const std::vector<Integer> &CFace<Real>::nodeIds() const {
        return _nodeIds;
    }

    template<typename Real>
    const Integer &CFace<Real>::nbNode() const {
        return _nbNode;
    }

    template<typename Real>
    void CFace<Real>::setPairFaceId(Integer aPairFaceId) {
        _pairFaceId = aPairFaceId;
        _hasPair = true;
    }

    template<typename Real>
    const Integer &CFace<Real>::pairFaceId() const {
        return _pairFaceId;
    }

    template<typename Real>
    bool CFace<Real>::hasPair() const {
        return _hasPair;
    }

    template<typename Real>
    void CFace<Real>::setCellId(Integer aCellId) {
        _cellId = aCellId;
    }

    template<typename Real>
    const Integer &CFace<Real>::cellId() const {
        return _cellId;
    }

    template<typename Real>
    void CFace<Real>::setNeighborId(Integer aNeighborId) {
        _neighborId = aNeighborId;
    }

    template<typename Real>
    const Integer &CFace<Real>::neighborId() const {
        return _neighborId;
    }

    template<typename Real>
    void CFace<Real>::setArea(const Real &anArea) {
        _area = anArea;
    }

    template<typename Real>
    const Real &CFace<Real>::area() const {
        return _area;
    }

    template<typename Real>
    void CFace<Real>::setNormal(const CVector<Real> &aNormalVector) {
        _normal = aNormalVector;
    }

    template<typename Real>
    const CVector<Real> &CFace<Real>::normal() const {
        return _normal;
    }

    template<typename Real>
    void CFace<Real>::setBarycenter(const CCoordinate<Real> &aPoint) {
        _barycenter = aPoint;
    }

    template<typename Real>
    const CCoordinate<Real> &CFace<Real>::barycenter() const {
        return _barycenter;
    }

    template<typename Real>
    void CFace<Real>::addQuadraturePt(const CCoordinate<Real> &aPoint) {
        _quadrature.push_back(aPoint);
        _nbQuadrature++;
    }

    template<typename Real>
    const CCoordinate<Real> &CFace<Real>::quadrature(Integer quadratureIndex) const {
        return _quadrature.at(quadratureIndex);
    }

    template<typename Real>
    const Integer &CFace<Real>::nbQuadrature() {
        return _nbQuadrature;
    }

    template<typename Real>
    void CFace<Real>::reset() {
        _faceType = FT_NONE;
        _boundaryType = BT_NONE;

        _nodeIds.clear();
        _nbNode = 0;

        _pairFaceId = 0;
        _cellId = 0;
        _neighborId = 0;
        _hasPair = false;

        _area = 0.0;
        _normal.setZero();
        _nbQuadrature = 0;
    }
}

#endif //HERTZ_CFACE_H
