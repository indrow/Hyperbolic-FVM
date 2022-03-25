//
// Created by indro on 10/03/21.
//

#ifndef HERTZ_CGRIDCREATE_H
#define HERTZ_CGRIDCREATE_H

#include "geo/CLine.h"
#include "CGrid.h"
#include "msh/CGmsh.h"

namespace hertz::fvm {
    using hertz::geo::CCoordinate;
    using hertz::geo::CLine;

    using hertz::msh::CGmsh;

    using hertz::fvm::CGrid;
    using hertz::fvm::BoundaryType;

    template<typename Real>
    class CGridCreate {
    private:
        CGmsh<Real> _cgmsh;
        CGrid<Real> _cgrid;

        std::map<std::vector<Integer>, BoundaryType> _faceBoundaryType;

        CFaceType _faceType;
        CCellType _cellType;

        CCoordinate<Real> _coordinate;
        CNode<Real> _node;
        CFace<Real> _face;
        CCell<Real> _cell;

        CFace<Real> *_facePtr, *_anotherFacePtr, *_localFacePtr;
        CCell<Real> *_cellPtr, *_localCellPtr;

        void populateNode();

        void domainBoundaries();

        void populateCells();

        void ghostCells();

        void translateCells(const Integer &toFace, const Integer &fromFace, const std::vector<Integer> &cells,
                            BoundaryType aType);

        void cellStencils();

        void geometry();

    public:
        CGridCreate();

        [[maybe_unused]] explicit CGridCreate(const CGmsh<Real> &aGmshMesh);

        void process(const CGmsh<Real> &aGmshMesh);

        const CGrid<Real> &grid();
    };

    template<typename Real>
    CGridCreate<Real>::CGridCreate() : _cgmsh(CGmsh<Real>()), _cgrid(CGrid<Real>()), _faceType(FT_NONE),
                                       _cellType(CT_NONE) {
    }

    template<typename Real>
    [[maybe_unused]] CGridCreate<Real>::CGridCreate(const CGmsh<Real> &aGmshMesh) : _faceType(FT_NONE),
                                                                                       _cellType(CT_NONE) {
        _cgmsh = aGmshMesh;
    }

    template<typename Real>
    void CGridCreate<Real>::process(const CGmsh<Real> &aGmshMesh) {
        _cgmsh = aGmshMesh;

        populateNode();
        domainBoundaries();
        populateCells();
        ghostCells();
        cellStencils();
        geometry();
    }

    template<typename Real>
    void CGridCreate<Real>::populateNode() {
        std::vector<Integer> nodeIds = _cgmsh.getNodeIds();
        std::vector<Real> nodesCoord = _cgmsh.getNodesCoord();

        /* collecting node */
        auto node = [this, &nodesCoord](const Integer &i) {
            _node = CNode<Real>::Map((std::vector<Real>(nodesCoord.begin() + i * 3,
                                                       nodesCoord.begin() + (i * 3) + 2)).data());
        };

        for (Integer i = 0; i < nodeIds.size(); ++i) {
            node(i);
            _cgrid.addNode(nodeIds.at(i), _node);
        }

        _cgrid.setInternalNode(_cgrid.nbNodes());
    }

    template<typename Real>
    void CGridCreate<Real>::domainBoundaries() {
        std::unordered_map<int, std::string> physicalGroupNames = _cgmsh.getPhysicalGroupNames();
        std::map<std::vector<Integer>, int> facesBoundaryType = _cgmsh.getFacesBoundaryType();

        std::unordered_map<int, BoundaryType> boundaryType;

        switch (_cgmsh.getFaceType()) {
            case 0:
                _faceType = FT_POINT;
                break;

            case 1:
                _faceType = FT_LINE;
                break;

            default:
                _faceType = FT_NONE;
        }

        switch (_cgmsh.getCellType()) {
            case 1:
                _cellType = CT_LINE;
                break;

            case 2:
                _cellType = CT_TRIANGLE;
                break;

            case 3:
                _cellType = CT_QUADRILATERAL;
                break;

            default:
                _cellType = CT_NONE;
        }

        /* all boundaries list */
        for (auto &pgn : physicalGroupNames) {
            std::transform(pgn.second.begin(), pgn.second.end(), pgn.second.begin(),
                           [](unsigned char c) { return std::tolower(c); });

            if (pgn.second == "periodic") {
                boundaryType[pgn.first] = BT_PERIODIC;
            }
        }

        /* faces at boundaries */
        for (const auto &faceBT : facesBoundaryType) {
            _faceBoundaryType[faceBT.first] = boundaryType.at(faceBT.second);
        }
    }

    template<typename Real>
    void CGridCreate<Real>::populateCells() {
        std::vector<Integer> cellIds = _cgmsh.getCellIds();
        std::vector<Integer> cellNodeIds = _cgmsh.getCellNodeIds();
        std::vector<Integer> cellFaces = _cgmsh.getCellsFace();

        std::vector<std::unordered_map<Integer, Integer>> periodicBCMap = _cgmsh.getPeriodicBCMap();

        Integer nbFaceNode = _faceType, nbCellNode = _cellType + 1, nbCellFaceNode = nbFaceNode * nbCellNode;

        auto collectCellNodeIds = [&cellNodeIds, &nbCellNode](const Integer &i) {
            return std::vector<Integer>(cellNodeIds.begin() + i * nbCellNode,
                                        cellNodeIds.begin() + (i + (Integer)1) * nbCellNode);
        };

        auto collectCellFaceIds = [&cellFaces, &nbCellFaceNode](const Integer &i) {
            return std::vector<Integer>(cellFaces.begin() + i * nbCellFaceNode,
                                        cellFaces.begin() + (i + (Integer)1) * nbCellFaceNode);
        };

        std::vector<Integer> sliceCellNodeIds, sliceCellFaces;

        auto aFace = [this, &cellIds, &sliceCellFaces, &nbFaceNode](const Integer &i, const Integer &j) {
            _face.setCellId(cellIds[i]);
            _face.setFaceType(_faceType);
            _face.addNodeId(sliceCellFaces[j * nbFaceNode]);
            _face.addNodeId(sliceCellFaces[j * nbFaceNode + 1]);
        };

        std::vector<Integer> faceNodeIds;

        std::map<std::vector<Integer>, std::vector<Integer>> neighboringMap;
        std::map<std::vector<Integer>, std::vector<Integer>> facePairIdMap;

        Integer nbQuadrature = _cgmsh.getNbQuadrature();
        std::vector<Real> quadrature = _cgmsh.getQuadrature();

        _cgrid.setCellQuadratureWeights(_cgmsh.getQuadratureWeights());

        auto collectCellQuadrature = [this, &quadrature](const Integer &counter) {
            _coordinate = CCoordinate<Real>::Map(std::vector<Real>(quadrature.begin() + counter * 3,
                                                                  quadrature.begin() + (counter * 3) + 2).data());
        };

        Integer quadratureCount = 0;

        /* main loop:
         *      - face nodes
         *      - cell nodes
         *      - cell quadrature / cubature
         * */
        for (Integer i = 0, faceId = 1, facePairId; i < cellIds.size(); ++i) {
            sliceCellNodeIds = collectCellNodeIds(i);
            sliceCellFaces = collectCellFaceIds(i);

            for (Integer j = 0; j < nbCellNode; ++j) {
                aFace(i, j);
                _cgrid.addFace(faceId, _face);
                _facePtr = &_cgrid.face(faceId);

                faceNodeIds = _face.nodeIds();

                if (_faceBoundaryType.find(faceNodeIds) != _faceBoundaryType.end()) {
                    _facePtr->setFaceBoundaryType(_faceBoundaryType[faceNodeIds]);
                } else {
                    _facePtr->setFaceBoundaryType(BT_INTERIOR);
                }

                std::sort(faceNodeIds.begin(), faceNodeIds.end());

                neighboringMap[faceNodeIds].push_back(cellIds[i]);
                facePairIdMap[faceNodeIds].push_back(faceId);

                if (neighboringMap[faceNodeIds].size() == 2) {
                    facePairId = facePairIdMap.at(faceNodeIds).front();
                    _anotherFacePtr = &_cgrid.face(facePairId);

                    _facePtr->setPairFaceId(facePairId);
                    _facePtr->setNeighborId(neighboringMap.at(faceNodeIds).front());

                    _anotherFacePtr->setPairFaceId(faceId);
                    _anotherFacePtr->setNeighborId(cellIds[i]);
                }

                _cell.addFaceId(faceId);
                _face.reset();

                faceId++;
            }

            std::for_each(sliceCellNodeIds.begin(), sliceCellNodeIds.end(),
                          [this](const Integer &node) { _cell.addNodeId(node); });

            for (Integer j = 0; j < nbQuadrature; ++j) {
                collectCellQuadrature(quadratureCount);
                _cell.addQuadrature(_coordinate);
                quadratureCount++;
            }

            _cell.setCellType(_cellType);
            _cell.setCellBoundaryType(BT_INTERIOR);
            _cgrid.addCell(cellIds.at(i), _cell);

            _cell.reset();
        }

        _cgrid.setInternalFace(_cgrid.nbFaces());
        _cgrid.setInternalCell(_cgrid.nbCells());

        Integer faceId, pairFaceId;
        std::vector<Integer> periodicPairIds;

        /* periodic boundary face mapping (necessary if periodic bc was detected) */
        if (!periodicBCMap.empty()) {
            goto search;
        } else {
            goto skip;
        }

        search:
        for (const auto &faceBT : _faceBoundaryType) {
            if (faceBT.second == BT_PERIODIC) {
                faceNodeIds = faceBT.first;

                std::sort(faceNodeIds.begin(), faceNodeIds.end());
                faceId = facePairIdMap.at(faceNodeIds).front();
                _facePtr = &_cgrid.face(faceId);

                for (const auto &periodicNode : periodicBCMap) {
                    for (const auto &aNode : _facePtr->nodeIds()) {
                        if (periodicNode.find(aNode) != periodicNode.end()) {
                            periodicPairIds.push_back(periodicNode.at(aNode));
                        }
                    }

                    if (periodicPairIds.size() == faceNodeIds.size()) {
                        std::sort(periodicPairIds.begin(), periodicPairIds.end());
                        pairFaceId = facePairIdMap.at(periodicPairIds).front();
                        _anotherFacePtr = &_cgrid.face(pairFaceId);

                        _facePtr->setPairFaceId(pairFaceId);
                        _facePtr->setNeighborId(_anotherFacePtr->cellId());

                        _anotherFacePtr->setPairFaceId(faceId);
                        _anotherFacePtr->setNeighborId(_facePtr->cellId());

                        periodicPairIds.clear();

                    } else {
                        periodicPairIds.clear();
                        continue;
                    }
                }
            }
        }

        skip:;
    }

    template<typename Real>
    void CGridCreate<Real>::ghostCells() {
        Integer target, next;
        std::vector<Integer> neighbors;

        CFace<Real> aFace;

        /* if cell face is a boundary face then create some ghost cells */
        auto collectNeighborOf = [this, &aFace, &neighbors](const Integer &id, const Integer &exceptionId) {

            _localCellPtr = &_cgrid.cell(id);

            for (Integer j = 0; j < _localCellPtr->nbFace(); ++j) {
                aFace = _cgrid.face(_localCellPtr->faceId(j));

                if (std::find(neighbors.begin(), neighbors.end(), aFace.neighborId()) == neighbors.end() &&
                    aFace.neighborId() != exceptionId && aFace.boundaryType() != BT_PERIODIC) {

                    neighbors.push_back(aFace.neighborId());
                }

                if (neighbors.size() > 6) break;
            }
        };

        /* for periodic bc use a translate cell function, else use a mirror cell function */
        for (Integer i = 0; i < _cgrid.nbInternalCell(); ++i) {
            _cellPtr = &_cgrid.cellByIndices(i);

            for (Integer j = 0; j < _cellPtr->nbFace(); ++j) {
                _face = _cgrid.face(_cellPtr->faceId(j));

                if (_face.boundaryType() == BT_PERIODIC) {
                    neighbors.push_back(_face.neighborId());
                    target = neighbors.front();
                    next = 1;

                    while (neighbors.size() < 7) {
                        collectNeighborOf(target, _cgrid.cellId(i));
                        target = neighbors.at(next);
                        next++;
                    }

                    translateCells(_cellPtr->faceId(j), _face.pairFaceId(), neighbors, BT_PERIODIC);
                    neighbors.clear();
                }
            }
        }
    }

    template<typename Real>
    void CGridCreate<Real>::translateCells(const Integer &toFace, const Integer &fromFace,
                                           const std::vector<Integer> &cells, BoundaryType aType) {
        /* note: this code is used specifically for fourth-order WENO schemes, stencil is in a fixed size */

        const Integer targetFaceId = toFace;
        const Integer originFaceId = fromFace;

        CFace<Real> face, pairFace;

        CCoordinate<Real> distance;
        CLine<Real> target, origin;

        auto createLine = [this](const Integer &aFaceId, CLine<Real> &aLine) {
            aLine.setStartPoint(_cgrid.node(_cgrid.face(aFaceId).nodeId(0)));
            aLine.setEndPoint(_cgrid.node(_cgrid.face(aFaceId).nodeId(1)));
        };

        createLine(originFaceId, origin);
        createLine(targetFaceId, target);

        distance = target.midPoint(0.5) - origin.midPoint(0.5);

        Integer lastNodeId, lastFaceId, lastCellId;
        std::unordered_map<Integer, Integer> nodeIds, faceIds, cellIds;

        lastNodeId = _cgrid.lastNodeId();
        lastFaceId = _cgrid.lastFaceId();
        lastCellId = _cgrid.lastCellId();

        /* list of required nodes, faces, and cells */
        std::for_each(cells.begin(), cells.end(), [this, &originFaceId, &targetFaceId, &lastNodeId, &lastFaceId,
                &lastCellId, &nodeIds, &faceIds, &cellIds](const Integer &aCellId) {

            if (cellIds.find(aCellId) == cellIds.end()) {
                lastCellId++;
                cellIds[aCellId] = lastCellId;
            }

            _localCellPtr = &_cgrid.cell(aCellId);

            for (Integer i = 0; i < _localCellPtr->nbFace(); ++i) {
                if (faceIds.find(_localCellPtr->faceId(i)) == faceIds.end()) {
                    lastFaceId++;
                    faceIds[_localCellPtr->faceId(i)] = lastFaceId;
                }

                if (nodeIds.find(_localCellPtr->nodeId(i)) == nodeIds.end()) {
                    lastNodeId++;
                    nodeIds[_localCellPtr->nodeId(i)] = lastNodeId;
                }
            }
        });

        std::vector<Integer> targetNodeIds = _cgrid.face(targetFaceId).nodeIds();

        /* some new nodes */
        std::for_each(nodeIds.begin(), nodeIds.end(), [this, &targetNodeIds, &distance](const auto &nodePair) {
            if (std::find(targetNodeIds.begin(), targetNodeIds.end(), nodePair.second) == targetNodeIds.end()) {
                _node = _cgrid.node(nodePair.first) + distance;
                _cgrid.addNode(nodePair.second, _node);
            }
        });

        /* some new faces */
        std::for_each(faceIds.begin(), faceIds.end(), [this, &originFaceId, &targetFaceId, &nodeIds, &faceIds, &cellIds]
                (const auto &facePair) {

            _face.reset();
            _localFacePtr = &_cgrid.face(facePair.first);

            for (Integer i = 0; i < _localFacePtr->nbNode(); ++i) {
                _face.addNodeId(nodeIds.at(_localFacePtr->nodeId(i)));
            }

            _face.setCellId(cellIds.at(_localFacePtr->cellId()));

            if (faceIds.find(_localFacePtr->pairFaceId()) != faceIds.end()) {
                _face.setPairFaceId(faceIds.at(_localFacePtr->pairFaceId()));

                if (cellIds.find(_cgrid.face(facePair.first).neighborId()) != cellIds.end()) {
                    _face.setNeighborId(cellIds.at(_localFacePtr->neighborId()));
                }
            }

            if (facePair.first == originFaceId) {
                _face.setNeighborId(_cgrid.face(targetFaceId).cellId());
                _cgrid.face(targetFaceId).setNeighborId(_face.cellId());
            }

            _cgrid.addFace(facePair.second, _face);
        });

        _cell.reset();

        /* some new cells */
        std::for_each(cellIds.begin(), cellIds.end(), [this, &distance, &nodeIds, &faceIds, aType](const auto &cellPair) {
            _localCellPtr = &_cgrid.cell(cellPair.first);

            for (Integer i = 0; i < _localCellPtr->nbNode(); ++i) {
                _cell.addNodeId(nodeIds.at(_localCellPtr->nodeId(i)));
                _cell.addFaceId(faceIds.at(_localCellPtr->faceId(i)));
            }

            for (Integer i = 0; i < _localCellPtr->nbQuadrature(); ++i) {
                _coordinate = _localCellPtr->quadrature(i) + distance;
                _cell.addQuadrature(_coordinate);
            }

            _cell.setCellBoundaryType(aType);
            _cell.setLinkedCellId(cellPair.first);
            _cgrid.addCell(cellPair.second, _cell);

            _cell.reset();
        });

        nodeIds.clear();
        faceIds.clear();
        cellIds.clear();

    }

    template<typename Real>
    void CGridCreate<Real>::cellStencils() {
        CFace<Real> aFace;

        auto collectNeighborOf = [this, &aFace](const Integer &id, const Integer &exceptionId,
                                                std::vector<Integer> &neighbors) {
            for (Integer j = 0; j < _cgrid.cell(id).nbFace(); ++j) {
                aFace = _cgrid.face(_cgrid.cell(id).faceId(j));

                if (aFace.neighborId() != exceptionId && aFace.neighborId() > 0) {

                    neighbors.push_back(aFace.neighborId());
                }
            }
        };

        Integer stencilNum;
        std::vector<Integer> stencil, bigStencil;
        std::vector<Integer> targetNeighbor, nextLayerNeighbor;

        stencil.reserve(4);
        bigStencil.reserve(10);

        /* list of all stencil in each cell id */
        auto stencilOf = [&, this](const Integer &cellId) {
            stencilNum = 0;

            stencil.push_back(cellId);
            bigStencil.push_back(cellId);

            _localCellPtr = &_cgrid.cell(cellId);

            collectNeighborOf(cellId, 0, targetNeighbor);

            for (const Integer &tn : targetNeighbor) {
                if (std::find(bigStencil.begin(), bigStencil.end(), tn) == bigStencil.end()) {
                    stencil.push_back(tn);
                    bigStencil.push_back(tn);

                    collectNeighborOf(tn, cellId, nextLayerNeighbor);

                    for (const Integer &sln : nextLayerNeighbor) {
                        if (std::find(bigStencil.begin(), bigStencil.end(), sln) == bigStencil.end()) {
                            stencil.push_back(sln);
                            bigStencil.push_back(sln);

                            _localCellPtr->addStencil(stencilNum, stencil);

                            stencil.erase(stencil.begin() + 2);
                            stencilNum++;
                        }
                    }
                    stencil.erase(stencil.begin() + 1);
                    nextLayerNeighbor.clear();
                }
            }

            if (bigStencil.size() > 10) {
                _localCellPtr->addStencil(stencilNum, bigStencil);
                goto clearance;

            } else {
                for (Integer n = 0; n < stencilNum; ++n) {
                    if (bigStencil.size() > 9) {
                        _localCellPtr->addStencil(stencilNum, bigStencil);
                        goto clearance;
                    } else {
                        stencil = _localCellPtr->stencil(n);
                        collectNeighborOf(stencil.back(), stencil.at(1), nextLayerNeighbor);

                        if (std::find(bigStencil.begin(), bigStencil.end(), nextLayerNeighbor.front()) ==
                            bigStencil.end()) {
                            stencil.push_back(nextLayerNeighbor.front());
                            bigStencil.push_back(nextLayerNeighbor.front());
                        }

                        _localCellPtr->addStencil(stencilNum, stencil);
                        stencilNum++;
                        nextLayerNeighbor.clear();
                    }
                }
            }

            clearance:
            stencil.clear();
            bigStencil.clear();
            targetNeighbor.clear();
        };

        /* main loop for stencil collection */
        for (Integer i = 0; i < _cgrid.nbInternalCell(); ++i) {
            stencilOf(_cgrid.cellId(i));
            _cellPtr = &_cgrid.cellByIndices(i);

            for (Integer j = 0; j < _cellPtr->nbFace(); ++j) {
                _facePtr = &_cgrid.face(_cellPtr->faceId(j));

                if (_facePtr->boundaryType() != BT_INTERIOR) {
                    stencilOf(_facePtr->neighborId());
                }
            }
        }
    }

    template<typename Real>
    void CGridCreate<Real>::geometry() {
        Integer nbNode;
        Real volume, semi_perimeter;
        CCoordinate<Real> barycenter;

        /* Calculating cell barycenter, volume, inner circle radius */
        for (Integer i = 0; i < _cgrid.nbCells(); ++i) {
            barycenter.setZero();
            volume = 0.0;
            semi_perimeter = 0.0;

            _cellPtr = &_cgrid.cellByIndices(i);
            nbNode = _cellPtr->nbNode();

            for (Integer j = 0; j < nbNode; ++j) {
                barycenter += _cgrid.node(_cellPtr->nodeId(j));

                volume += _cgrid.node(_cellPtr->nodeId(j)).x() * _cgrid.node(_cellPtr->nodeId((j + 1) % nbNode)).y() -
                          _cgrid.node(_cellPtr->nodeId((j + 1) % nbNode)).x() * _cgrid.node(_cellPtr->nodeId(j)).y();

                semi_perimeter += (_cgrid.node(_cellPtr->nodeId(j)) - _cgrid.node(_cellPtr->nodeId((j + 1) % nbNode))).norm();

            }

            barycenter /= (Real) nbNode;
            volume = fabs(volume) / 2.0;

            _cellPtr->setBarycenter(barycenter);
            _cellPtr->setVolume(volume);
            _cellPtr->setH(sqrt(volume));
            _cellPtr->setRadii(volume / semi_perimeter / 0.5);
        }

        /* Face quadrature points */
        Integer faceId;
        CLine<Real> line;
        CVector<Real> geoVector;
        std::unordered_map<Integer, bool> isCalculated;

        _cgrid.setFaceQuadratureWeights({0.5, 0.5});
        Real r1 = 0.5 - std::sqrt(3.0) / 6.0, r2 = 0.5 + std::sqrt(3.0) / 6.0;

        for (Integer i = 0; i < _cgrid.nbFaces(); ++i) {
            isCalculated[_cgrid.faceId(i)] = false;
        }

        for (Integer i = 0; i < _cgrid.nbFaces(); ++i) {
            faceId = _cgrid.faceId(i);
            _facePtr = &_cgrid.face(faceId);

            if (!isCalculated[_cgrid.faceId(i)]) {
                if (_facePtr->pairFaceId() > 0) {
                    _anotherFacePtr = &_cgrid.face(_facePtr->pairFaceId());

                } else
                    _anotherFacePtr = nullptr;

                line.setStartPoint(_cgrid.node(_facePtr->nodeId(0)));
                line.setEndPoint(_cgrid.node(_facePtr->nodeId(1)));

                _facePtr->addQuadraturePt(line.midPoint(r1));
                _facePtr->addQuadraturePt(line.midPoint(r2));
                _facePtr->setArea(line.length());

                geoVector = line.vector();
                geoVector = {geoVector.y(), -geoVector.x()};
                geoVector /= geoVector.norm();

                _facePtr->setNormal(geoVector);

                isCalculated[faceId] = true;

                if (_anotherFacePtr != nullptr && !isCalculated[_facePtr->pairFaceId()]) {
                    _anotherFacePtr->addQuadraturePt(_facePtr->quadrature(0));
                    _anotherFacePtr->addQuadraturePt(_facePtr->quadrature(1));
                    _anotherFacePtr->setArea(_facePtr->area());

                    geoVector = -geoVector;
                    _anotherFacePtr->setNormal(geoVector);

                    isCalculated[_facePtr->pairFaceId()] = true;
                }
            }
        }
    }

    template<typename Real>
    const CGrid<Real> &CGridCreate<Real>::grid() {
        return _cgrid;
    }
}


#endif //HERTZ_CGRIDCREATE_H
