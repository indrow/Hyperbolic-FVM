//
// Created by indro on 10/03/21.
//

#ifndef HERTZ_CGMSH_H
#define HERTZ_CGMSH_H

#include <fstream>
#include <gmsh.h>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>

#include "common/CType.h"

namespace hertz::msh {
    template<typename Real>
    class CGmsh {
    private:
        int _dimension;

        std::vector<Integer> _nodeIds;
        std::vector<Real> _nodeCoords;

        std::vector<std::unordered_map<Integer, Integer>> _periodicBCMap;
        Integer _nbPeriodicBC;

        std::unordered_map<int, std::string> _physicalGroupNames;
        std::map<std::vector<Integer>, int> _facesBoundaryType;

        std::vector<Real> _quadratureWeights;
        std::vector<Real> _quadratureLocalCoord;
        std::vector<Real> _quadratureCoord;

        std::vector<Integer> _cellIds;
        std::vector<Integer> _cellsFace;
        std::vector<Integer> _cellsNodeIds;

        int _faceType;
        int _cellType;
        Integer _nbQuadrature;

        void loadPeriodicity();

        void loadBoundary();

        void loadCellsData();

    public:
        CGmsh();

        void load(const std::string &strFilename);

        [[nodiscard]] const int &dimension() const;

        [[nodiscard]] const std::vector<Integer> &getNodeIds() const;

        [[nodiscard]] const std::vector<Real> &getNodesCoord() const;

        [[nodiscard]] const std::vector<std::unordered_map<Integer, Integer>> &getPeriodicBCMap() const;

        [[nodiscard]] const Integer &getNbPeriodicBC() const;

        [[nodiscard]] const std::unordered_map<int, std::string> &getPhysicalGroupNames() const;

        [[nodiscard]] const std::map<std::vector<Integer>, int> &getFacesBoundaryType() const;

        [[nodiscard]] const int &getFaceType() const;

        [[nodiscard]] const std::vector<Integer> &getCellIds() const;

        [[nodiscard]] const std::vector<Integer> &getCellsFace() const;

        [[nodiscard]] const std::vector<Integer> &getCellNodeIds() const;

        [[nodiscard]] const int &getCellType() const;

        [[nodiscard]] const Integer &getNbQuadrature() const;

        [[nodiscard]] const std::vector<Real> &getQuadratureWeights() const;

        [[nodiscard]] const std::vector<Real> &getQuadrature() const;

        void reset();

        ~CGmsh();

    };

    template<typename Real>
    CGmsh<Real>::CGmsh() : _dimension(0), _nbPeriodicBC(0), _faceType(0), _cellType(0), _nbQuadrature(0) {
    }

    template<typename Real>
    void CGmsh<Real>::load(const std::string &strFilename) {
        gmsh::initialize();
        gmsh::open(strFilename);

        _dimension = gmsh::model::getDimension();

        if (_dimension > 2) {
            gmsh::logger::write("Sorry, the current solver support only 2D geometry!", "error");
            gmsh::finalize();
            exit(-1);
        }

        std::vector<Real> parametricCoords;

        gmsh::model::mesh::getNodes(_nodeIds, _nodeCoords, parametricCoords, -1, -1, false, false);

        if (_nodeIds.empty()) {
            gmsh::logger::write("Something wrong with your mesh!", "error");
            gmsh::finalize();
            exit(-1);
        }

        loadPeriodicity();
        loadBoundary();
        loadCellsData();

        gmsh::finalize();
    }

    template<typename Real>
    void CGmsh<Real>::loadPeriodicity() {
        int idMaster;
        std::vector<Integer> nodeIds, nodeIdsMaster;
        std::vector<Real> affineTransform;

        std::unordered_map<Integer, Integer> periodicBCMap;

        for (int i = 1, j = 2; i < j; ++i) {
            gmsh::model::mesh::getPeriodicNodes(1, i, idMaster, nodeIds, nodeIdsMaster, affineTransform);

            if (nodeIds.empty()) {
                break;

            } else {
                for (Integer k = 0; k < nodeIds.size(); ++k) {
                    periodicBCMap[nodeIds[k]] = nodeIdsMaster[k];
                    periodicBCMap[nodeIdsMaster[k]] = nodeIds[k];
                }

                _periodicBCMap.push_back(periodicBCMap);
                _nbPeriodicBC++;

                j++;
                nodeIds.clear();
                nodeIdsMaster.clear();
                periodicBCMap.clear();
            }
        }
    }

    template<typename Real>
    void CGmsh<Real>::loadBoundary() {
        int dimension, order, nbNodes, numPrimaryNodes;
        std::string elementName;
        std::vector<double> localCoord;

        gmsh::vectorpair dimTags;
        std::vector<int> entityTags, elementTypes;
        std::vector<Integer> cellIds, nodeIds, nodeId;
        std::string physicalGroupNames;

        gmsh::model::getPhysicalGroups(dimTags, -1);

        for (const auto &dimTag : dimTags) {
            if (dimTag.first == 1) {
                gmsh::model::getEntitiesForPhysicalGroup(dimTag.first, dimTag.second, entityTags);

                for (const auto &entityTag : entityTags) {
                    gmsh::model::mesh::getElementTypes(elementTypes, dimTag.first, entityTag);

                    for (const auto &elementType : elementTypes) {
                        gmsh::model::getPhysicalName(dimTag.first, dimTag.second, physicalGroupNames);
                        _physicalGroupNames[dimTag.second] = physicalGroupNames;

                        gmsh::model::mesh::preallocateElementsByType(elementType, true, true, cellIds, nodeIds,
                                                                     entityTag);

                        gmsh::model::mesh::getElementProperties(elementType, elementName, dimension, order, nbNodes,
                                                                localCoord, numPrimaryNodes);

                        gmsh::model::mesh::getElementsByType(elementType, cellIds, nodeIds, entityTag, 0, 1);

                        for (Integer i = 0; i < nodeIds.size(); i += nbNodes) {
                            nodeId = std::vector<Integer>(nodeIds.begin() + i, nodeIds.begin() + i + nbNodes);
                            _facesBoundaryType[nodeId] = dimTag.second;
                        }
                    }
                }
            }
        }
    }

    template<typename Real>
    void CGmsh<Real>::loadCellsData() {
        std::vector<int> types;
        std::vector<Real> jacobians, determinants;

        gmsh::model::mesh::getElementTypes(types, 1, -1);
        _faceType = types.front();

        gmsh::model::mesh::getElementTypes(types, 2, -1);

        if (types.size() > 1) {
            gmsh::logger::write("Hybrid meshes not handled in this code!", "error");
            gmsh::finalize();
            exit(-1);
        } else {
            _cellType = types.front();

            gmsh::model::mesh::getElementEdgeNodes(_cellType, _cellsFace, -1, true, 0, 1);
            gmsh::model::mesh::getElementsByType(_cellType, _cellIds, _cellsNodeIds, -1, 0, 1);
        }

        gmsh::model::mesh::getIntegrationPoints(_cellType, "Gauss4", _quadratureLocalCoord, _quadratureWeights);
        gmsh::model::mesh::getJacobians(_cellType, _quadratureLocalCoord, jacobians,
                                        determinants, _quadratureCoord, -1, 0, 1);


        _nbQuadrature = _quadratureWeights.size();

        for (Integer i = 0; i < _nbQuadrature; ++i) {
            _quadratureWeights[i] = 2.0 * _quadratureWeights[i];
        }
    }

    template<typename Real>
    const int &CGmsh<Real>::dimension() const {
        return _dimension;
    }

    template<typename Real>
    const std::vector<Integer> &CGmsh<Real>::getNodeIds() const {
        return _nodeIds;
    }

    template<typename Real>
    const std::vector<Real> &CGmsh<Real>::getNodesCoord() const {
        return _nodeCoords;
    }

    template<typename Real>
    [[maybe_unused]] const std::vector<std::unordered_map<Integer, Integer>> &CGmsh<Real>::getPeriodicBCMap() const {
        return _periodicBCMap;
    }

    template<typename Real>
    const Integer &CGmsh<Real>::getNbPeriodicBC() const {
        return _nbPeriodicBC;
    }

    template<typename Real>
    const std::unordered_map<int, std::string> &CGmsh<Real>::getPhysicalGroupNames() const {
        return _physicalGroupNames;
    }

    template<typename Real>
    const std::map<std::vector<Integer>, int> &CGmsh<Real>::getFacesBoundaryType() const {
        return _facesBoundaryType;
    }

    template<typename Real>
    const int &CGmsh<Real>::getFaceType() const {
        return _faceType;
    }

    template<typename Real>
    const std::vector<Integer> &CGmsh<Real>::getCellIds() const {
        return _cellIds;
    }

    template<typename Real>
    const std::vector<Integer> &CGmsh<Real>::getCellsFace() const {
        return _cellsFace;
    }

    template<typename Real>
    const std::vector<Integer> &CGmsh<Real>::getCellNodeIds() const {
        return _cellsNodeIds;
    }

    template<typename Real>
    const int &CGmsh<Real>::getCellType() const {
        return _cellType;
    }

    template<typename Real>
    const Integer &CGmsh<Real>::getNbQuadrature() const {
        return _nbQuadrature;
    }

    template<typename Real>
    const std::vector<Real> &CGmsh<Real>::getQuadratureWeights() const {
        return _quadratureWeights;
    }

    template<typename Real>
    const std::vector<Real> &CGmsh<Real>::getQuadrature() const {
        return _quadratureCoord;
    }

    template<typename Real>
    void CGmsh<Real>::reset() {
        _dimension = 0;

        _nodeIds.clear();
        _nodeCoords.clear();

        _periodicBCMap.clear();
        _nbPeriodicBC = 0;

        _physicalGroupNames.clear();
        _facesBoundaryType.clear();

        _cellIds.clear();
        _cellsFace.clear();
        _cellsNodeIds.clear();

        _faceType = 0;
        _cellType = 0;
    }

    template<typename Real>
    CGmsh<Real>::~CGmsh() = default;

}

#endif //HERTZ_CGMSH_H
