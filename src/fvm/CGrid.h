//
// Created by indro on 10/03/21.
//

#ifndef HERTZ_CGRID_H
#define HERTZ_CGRID_H

#include "CCell.h"
#include "CFace.h"
#include "CNode.h"

namespace hertz::fvm {
    using hertz::fvm::CNode;
    using hertz::fvm::CFace;
    using hertz::fvm::CCell;

    template<typename Real>
    class CGrid {
    private:
        std::vector<Integer> _nodeIds;
        std::vector<Integer> _faceIds;
        std::vector<Integer> _cellIds;

        std::unordered_map<Integer, Integer> _nodeIndices;
        std::unordered_map<Integer, Integer> _faceIndices;
        std::unordered_map<Integer, Integer> _cellIndices;

        std::unordered_map<Integer, CNode<Real>> _nodes;
        std::unordered_map<Integer, CFace<Real>> _faces;
        std::unordered_map<Integer, CCell<Real>> _cells;

        Integer _nodeIndex;
        Integer _faceIndex;
        Integer _cellIndex;

        Integer _internalNodeIndex;
        Integer _internalFaceIndex;
        Integer _internalCellIndex;

        Integer _lastNodeId;
        Integer _lastFaceId;
        Integer _lastCellId;

        std::vector<Real> _faceQuadratureWeights, _cellQuadratureWeights;
        Integer _nbFaceQuadrature, _nbCellQuadrature;

    public:
        CGrid();

        void addNode(const Integer &aNodeId, const CNode<Real> &aNode);

        void addFace(const Integer &aFaceId, const CFace<Real> &aFace);

        void addCell(const Integer &aCellId, const CCell<Real> &aCell);

        void removeNode(const Integer &aNodeId);

        void removeFace(const Integer &aFaceId);

        void removeCell(const Integer &aCellId);

        const Integer &nodeId(const Integer &anIndex) const;

        const Integer &faceId(const Integer &anIndex) const;

        const Integer &cellId(const Integer &anIndex) const;

        CNode<Real> &node(const Integer &aNodeId);

        CFace<Real> &face(const Integer &aFaceId);

        CCell<Real> &cell(const Integer &aCellId);

        CNode<Real> &nodeByIndices(const Integer &aNodeIndices);

        CFace<Real> &faceByIndices(const Integer &aFaceIndices);

        CCell<Real> &cellByIndices(const Integer &aCellIndices);

        void setInternalNode(Integer anIndex);

        void setInternalFace(Integer anIndex);

        void setInternalCell(Integer anIndex);

        const Integer &nbInternalNode() const;

        const Integer &nbInternalFace() const;

        const Integer &nbInternalCell() const;

        const Integer &nbNodes() const;

        const Integer &nbFaces() const;

        const Integer &nbCells() const;

        const Integer &lastNodeId() const;

        const Integer &lastFaceId() const;

        const Integer &lastCellId() const;

        void setFaceQuadratureWeights(std::vector<Real> weights);

        const std::vector<Real> &faceQuadratureWeights() const;

        const Integer &nbFaceQuadrature() const;

        void setCellQuadratureWeights(std::vector<Real> weights);

        const std::vector<Real> &cellQuadratureWeights() const;

        const Integer &nbCellQuadrature() const;


        void reset();


    };

    template<typename Real>
    CGrid<Real>::CGrid() : _nodeIndex(0), _faceIndex(0), _cellIndex(0), _internalNodeIndex(0), _internalFaceIndex(0),
                           _internalCellIndex(0), _lastNodeId(0), _lastFaceId(0), _lastCellId(0) {
    }

    template<typename Real>
    void CGrid<Real>::addNode(const Integer &aNodeId, const CNode<Real> &aNode) {
        _nodeIds.push_back(aNodeId);

        _nodeIndices[aNodeId] = _nodeIndex;
        _nodeIndex++;

        _nodes[aNodeId] = aNode;

        _lastNodeId = aNodeId > _lastNodeId ? aNodeId : _lastNodeId;
    }

    template<typename Real>
    void CGrid<Real>::addFace(const Integer &aFaceId, const CFace<Real> &aFace) {
        _faceIds.push_back(aFaceId);

        _faceIndices[aFaceId] = _faceIndex;
        _faceIndex++;

        _faces[aFaceId] = aFace;

        _lastFaceId = aFaceId > _lastFaceId ? aFaceId : _lastFaceId;
    }

    template<typename Real>
    void CGrid<Real>::addCell(const Integer &aCellId, const CCell<Real> &aCell) {
        _cellIds.push_back(aCellId);

        _cellIndices[aCellId] = _cellIndex;
        _cellIndex++;

        _cells[aCellId] = aCell;

        _lastCellId = aCellId > _lastCellId ? aCellId : _lastCellId;
    }

    template<typename Real>
    void CGrid<Real>::removeNode(const Integer &aNodeId) {
        _nodeIds.erase(std::find(_nodeIds.begin(), _nodeIds.end(), aNodeId));
        _nodeIndices.erase(aNodeId);
        _nodes.erase(aNodeId);

        _nodeIndex--;
    }

    template<typename Real>
    void CGrid<Real>::removeFace(const Integer &aFaceId) {
        _faceIds.erase(std::find(_faceIds.begin(), _faceIds.end(), aFaceId));
        _faceIndices.erase(aFaceId);
        _faces.erase(aFaceId);

        _faceIndex--;
    }

    template<typename Real>
    void CGrid<Real>::removeCell(const Integer &aCellId) {
        _cellIds.erase(std::find(_cellIds.begin(), _cellIds.end(), aCellId));
        _cellIndices.erase(aCellId);
        _cells.erase(aCellId);

        _cellIndex--;
    }

    template<typename Real>
    const Integer &CGrid<Real>::nodeId(const Integer &anIndex) const {
        return _nodeIds.at(anIndex);
    }

    template<typename Real>
    const Integer &CGrid<Real>::faceId(const Integer &anIndex) const {
        return _faceIds.at(anIndex);
    }

    template<typename Real>
    const Integer &CGrid<Real>::cellId(const Integer &anIndex) const {
        return _cellIds.at(anIndex);
    }

    template<typename Real>
    CNode<Real> &CGrid<Real>::node(const Integer &aNodeId) {
        return _nodes.at(aNodeId);
    }

    template<typename Real>
    CFace<Real> &CGrid<Real>::face(const Integer &aFaceId) {
        return _faces.at(aFaceId);
    }

    template<typename Real>
    CCell<Real> &CGrid<Real>::cell(const Integer &aCellId) {
        return _cells.at(aCellId);
    }

    template<typename Real>
    CNode<Real> &CGrid<Real>::nodeByIndices(const Integer &aNodeIndices) {
        return _nodes.at(_nodeIds.at(aNodeIndices));
    }

    template<typename Real>
    CFace<Real> &CGrid<Real>::faceByIndices(const Integer &aFaceIndices) {
        return _faces.at(_faceIds.at(aFaceIndices));
    }

    template<typename Real>
    CCell<Real> &CGrid<Real>::cellByIndices(const Integer &aCellIndices) {
        return _cells.at(_cellIds.at(aCellIndices));
    }

    template<typename Real>
    void CGrid<Real>::setInternalNode(Integer anIndex) {
        _internalNodeIndex = anIndex;
    }

    template<typename Real>
    void CGrid<Real>::setInternalFace(Integer anIndex) {
        _internalFaceIndex = anIndex;
    }

    template<typename Real>
    void CGrid<Real>::setInternalCell(Integer anIndex) {
        _internalCellIndex = anIndex;
    }

    template<typename Real>
    const Integer &CGrid<Real>::nbInternalNode() const {
        return _internalNodeIndex;
    }

    template<typename Real>
    const Integer &CGrid<Real>::nbInternalFace() const {
        return _internalFaceIndex;
    }

    template<typename Real>
    const Integer &CGrid<Real>::nbInternalCell() const {
        return _internalCellIndex;
    }

    template<typename Real>
    const Integer &CGrid<Real>::nbNodes() const {
        return _nodeIndex;
    }

    template<typename Real>
    const Integer &CGrid<Real>::nbFaces() const {
        return _faceIndex;
    }

    template<typename Real>
    const Integer &CGrid<Real>::nbCells() const {
        return _cellIndex;
    }

    template<typename Real>
    const Integer &CGrid<Real>::lastNodeId() const {
        return _lastNodeId;
    }

    template<typename Real>
    const Integer &CGrid<Real>::lastFaceId() const {
        return _lastFaceId;
    }


    template<typename Real>
    const Integer &CGrid<Real>::lastCellId() const {
        return _lastCellId;
    }

    template<typename Real>
    void CGrid<Real>::setFaceQuadratureWeights(std::vector<Real> weights) {
        _faceQuadratureWeights = weights;
        _nbFaceQuadrature = weights.size();
    }

    template<typename Real>
    const std::vector<Real> &CGrid<Real>::faceQuadratureWeights() const {
        return _faceQuadratureWeights;
    }

    template<typename Real>
    void CGrid<Real>::setCellQuadratureWeights(std::vector<Real> weights) {
        _cellQuadratureWeights = weights;
        _nbCellQuadrature = weights.size();
    }

    template<typename Real>
    const std::vector<Real> &CGrid<Real>::cellQuadratureWeights() const {
        return _cellQuadratureWeights;
    }

    template<typename Real>
    void CGrid<Real>::reset() {
        _nodeIds.clear();
        _faceIds.clear();
        _cellIds.clear();

        _nodeIndices.clear();
        _faceIndices.clear();
        _cellIndices.clear();

        _nodes.clear();
        _faces.clear();
        _cells.clear();

        _nodeIndex = 0;
        _faceIndex = 0;
        _cellIndex = 0;

        _internalNodeIndex = 0;
        _internalFaceIndex = 0;
        _internalCellIndex = 0;

        _lastNodeId = 0;
        _lastFaceId = 0;
        _lastCellId = 0;

        _cellQuadratureWeights.clear();
    }

    template<typename Real>
    const Integer &CGrid<Real>::nbCellQuadrature() const {
        return _nbCellQuadrature;
    }

    template<typename Real>
    const Integer &CGrid<Real>::nbFaceQuadrature() const {
        return _nbFaceQuadrature;
    }
}


#endif //HERTZ_CGRID_H
