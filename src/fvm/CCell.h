//
// Created by indro on 10/03/21.
//

#ifndef HERTZ_CCELL_H
#define HERTZ_CCELL_H

#include <unordered_map>

#include "CFace.h"

namespace hertz::fvm {
    using hertz::fvm::BoundaryType;

    enum CCellType {
        CT_NONE,
        CT_LINE,
        CT_TRIANGLE,
        CT_QUADRILATERAL
    };

    class CCellBoundaryType {
    private:
        BoundaryType _boundaryType;

    public:
        CCellBoundaryType() : _boundaryType(BT_NONE) {}

        explicit CCellBoundaryType(const BoundaryType boundaryType) : _boundaryType(boundaryType) {}

        CCellBoundaryType &operator=(const BoundaryType &boundaryType) {
            _boundaryType = boundaryType;
            return *this;
        }

        BoundaryType operator()() const {
            return _boundaryType;
        }

        bool operator==(const BoundaryType &boundaryType) const {
            return _boundaryType == boundaryType;
        }

        bool operator!=(const BoundaryType &boundaryType) const {
            return _boundaryType != boundaryType;
        }
    };

    template<typename Real>
    class CCell {
    private:
        CCellType _cellType;
        CCellBoundaryType _boundaryType;

        std::vector<Integer> _nodeIds;
        std::vector<Integer> _faceIds;
        Integer _nbNode, _nbFace;

        std::unordered_map<Integer, std::vector<Integer>> _stencil;
        std::unordered_map<Integer, Integer> _nbCellInStencil;
        Integer _nbStencil;

        std::vector<CCoordinate < Real>> _quadrature;
        CCoordinate <Real> _aPoint;
        CCoordinate <Real> _barycenter, _newPoint;
        Integer _nbQuadrature;

        Integer _linkedCellId;

        Real _volume, _h, _radii;

    public:
        CCell();

        void setCellType(CCellType aCellType);

        [[nodiscard]] const CCellType &cellType() const;

        void setCellBoundaryType(BoundaryType aBoundaryType);

        [[nodiscard]] const CCellBoundaryType &boundaryType() const;

        void addNodeId(const Integer &aNodeId);

        [[nodiscard]] const Integer &nodeId(const Integer &aNodeIndex) const;

        [[nodiscard]] const Integer &nbNode() const;

        void addFaceId(const Integer &aFaceId);

        [[nodiscard]] const Integer &faceId(const Integer &aFaceIndex) const;

        [[nodiscard]] const Integer &nbFace() const;

        void addStencil(const Integer &anIndex, const std::vector<Integer> &aStencil);

        void replaceStencil(const Integer &anIndex, const std::vector<Integer> &aStencil);

        [[nodiscard]] const std::vector<Integer> &stencil(const Integer &anIndex) const;

        [[nodiscard]] const Integer &nbCellInStencil(const Integer &anIndex) const;

        [[nodiscard]] const Integer &nbStencil() const;

        void addQuadrature(const CCoordinate <Real> &aQuadraturePoint);

        [[nodiscard]] const CCoordinate <Real> &quadrature(const Integer &anIndex);

        [[nodiscard]] const Integer &nbQuadrature();

        void setBarycenter(const CCoordinate <Real> &aPoint);

        [[nodiscard]] const CCoordinate <Real> &barycenter() const;

        void setLinkedCellId(const Integer &aCellId);

        const Integer &linkedCellId();

        void setVolume(const Real &theVolume);

        const Real &volume() const;

        void setH(const Real &theH);

        const Real &h() const;

        void setRadii(const Real &theRadii);

        const Real &radii() const;

        CCoordinate <Real> scale(const CCoordinate <Real> &origin);

        void reset();
    };

/*============================ Implementation =======================*/

    template<typename Real>
    CCell<Real>::CCell() : _cellType(CT_NONE), _boundaryType(BT_NONE), _nbNode(0), _nbFace(0), _nbStencil(0),
                           _nbQuadrature(0), _linkedCellId(0) {
    }

    template<typename Real>
    void CCell<Real>::setCellType(CCellType aCellType) {
        _cellType = aCellType;
    }

    template<typename Real>
    const CCellType &CCell<Real>::cellType() const {
        return _cellType;
    }

    template<typename Real>
    void CCell<Real>::setCellBoundaryType(BoundaryType aBoundaryType) {
        _boundaryType = aBoundaryType;
    }

    template<typename Real>
    const CCellBoundaryType &CCell<Real>::boundaryType() const {
        return _boundaryType;
    }

    template<typename Real>
    void CCell<Real>::addNodeId(const Integer &aNodeId) {
        _nodeIds.push_back(aNodeId);
        _nbNode++;
    }

    template<typename Real>
    const Integer &CCell<Real>::nodeId(const Integer &aNodeIndex) const {
        return _nodeIds.at(aNodeIndex);
    }

    template<typename Real>
    const Integer &CCell<Real>::nbNode() const {
        return _nbNode;
    }

    template<typename Real>
    void CCell<Real>::addFaceId(const Integer &aFaceId) {
        _faceIds.push_back(aFaceId);
        _nbFace++;
    }

    template<typename Real>
    const Integer &CCell<Real>::faceId(const Integer &aFaceIndex) const {
        return _faceIds.at(aFaceIndex);
    }

    template<typename Real>
    const Integer &CCell<Real>::nbFace() const {
        return _nbFace;
    }

    template<typename Real>
    void CCell<Real>::addStencil(const Integer &anIndex, const std::vector<Integer> &aStencil) {
        _stencil[anIndex] = aStencil;
        _nbCellInStencil[anIndex] = aStencil.size();
        _nbStencil++;
    }

    template<typename Real>
    void CCell<Real>::replaceStencil(const Integer &anIndex, const std::vector<Integer> &aStencil) {
        _stencil.at(anIndex) = aStencil;
        _nbCellInStencil.at(anIndex) = aStencil.size();
    }

    template<typename Real>
    const std::vector<Integer> &CCell<Real>::stencil(const Integer &anIndex) const {
        return _stencil.at(anIndex);
    }

    template<typename Real>
    const Integer &CCell<Real>::nbCellInStencil(const Integer &anIndex) const {
        return _nbCellInStencil.at(anIndex);
    }

    template<typename Real>
    const Integer &CCell<Real>::nbStencil() const {
        return _nbStencil;
    }

    template<typename Real>
    void CCell<Real>::addQuadrature(const CCoordinate <Real> &aQuadraturePoint) {
        _quadrature.push_back(aQuadraturePoint);
        _nbQuadrature++;
    }

    template<typename Real>
    const CCoordinate <Real> &CCell<Real>::quadrature(const Integer &anIndex) {
        return _quadrature.at(anIndex);
    }

    template<typename Real>
    const Integer &CCell<Real>::nbQuadrature() {
        return _nbQuadrature;
    }

    template<typename Real>
    void CCell<Real>::setBarycenter(const CCoordinate <Real> &aPoint) {
        _barycenter = aPoint;
    }

    template<typename Real>
    const CCoordinate <Real> &CCell<Real>::barycenter() const {
        return _barycenter;
    }

    template<typename Real>
    void CCell<Real>::setLinkedCellId(const Integer &aCellId) {
        _linkedCellId = aCellId;
    }

    template<typename Real>
    const Integer &CCell<Real>::linkedCellId() {
        return _linkedCellId;
    }

    template<typename Real>
    void CCell<Real>::reset() {
        _cellType = CT_NONE;
        _boundaryType = BT_NONE;

        _nodeIds.clear();
        _faceIds.clear();
        _nbNode = 0, _nbFace = 0;

        _stencil.clear();
        _nbCellInStencil.clear();

        _quadrature.clear();
        _barycenter = CCoordinate<Real>(0);
        _nbQuadrature = 0;

        _linkedCellId = 0;
    }

    template<typename Real>
    void CCell<Real>::setVolume(const Real &theVolume) {
        _volume = theVolume;
    }

    template<typename Real>
    void CCell<Real>::setH(const Real &theH) {
        _h = theH;
    }

    template<typename Real>
    void CCell<Real>::setRadii(const Real &theRadii) {
        _radii = theRadii;
    }

    template<typename Real>
    const Real &CCell<Real>::volume() const {
        return _volume;
    }

    template<typename Real>
    const Real &CCell<Real>::h() const {
        return _h;
    }

    template<typename Real>
    const Real &CCell<Real>::radii() const {
        return _radii;
    }

    template<typename Real>
    CCoordinate <Real> CCell<Real>::scale(const CCoordinate <Real> &origin) {
        _newPoint = (origin - _barycenter) / _h;
        return _newPoint;
    }
}

#endif //HERTZ_CCELL_H
