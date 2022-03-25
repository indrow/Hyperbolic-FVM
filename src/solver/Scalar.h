//
// Created by indro on 10/03/21.
//

#ifndef HERTZ_SOLVER_SCALAR_H
#define HERTZ_SOLVER_SCALAR_H

#include "fvm/CBasis.h"
#include "fvm/CGridCreate.h"
#include "fvm/DCell.h"
#include "fvm/DFace.h"
#include "msh/CGmsh.h"
#include "physics/Scalar.h"
#include "schemes/RK.h"

#include <limits>
#include <memory>
#include <string>

namespace hertz::solver {
    using Eigen::Matrix;
    using Eigen::Dynamic;

    using hertz::msh::CGmsh;

    using hertz::geo::CCoordinate;
    using hertz::geo::CVector;
    using hertz::geo::CLine;

    using hertz::fvm::CBasis;
    using hertz::fvm::CGridCreate;
    using hertz::fvm::CGrid;
    using hertz::fvm::CFace;
    using hertz::fvm::CCell;
    using hertz::fvm::DCell;
    using hertz::fvm::DFace;
    using hertz::fvm::BoundaryType;

    using hertz::schemes::RK;

    template<typename Real>
    class Scalar {
    public:
        typedef Matrix<Real, Dynamic, 1> VectorXr;

        explicit Scalar(std::string  filename);

        void setFluxFunctions(Real (*f)(const Real &), Real (*df)(const Real &));

        void initVar();

        void initialize(Real (*u0)(const CCoordinate<Real> &), Real CFL, Real finalTime);

        void boundary();

        void construction();

        void timeStep();

        void timeStepping();

        void write();

        ~Scalar();

    private:
        CGrid<Real> _grid;

        Real (*_u0)(const CCoordinate<Real> &), (*_f)(const Real &), (*_df)(const Real &);

        Real _cfl, _finalTime, _t, _dt;

        Integer _cellId;

        CFace<Real> *_facePtr;
        CCell<Real> *_cellPtr, *_anotherCellPtr;

        DFace<Real, Real> *_dFacePtr, *_dFacePairPtr;
        DCell<Real, Real> *_dCellPtr, *_dCellNeighborPtr;

        std::unordered_map<Integer, DFace<Real, Real>> _dFaces;
        std::unordered_map<Integer, DCell<Real, Real>> _dCells;

        CBasis<Real> _basis;
        VectorXr _basisP3, _basisP3_0;

        const Real _inf;

    };

    template<typename Real>
    Scalar<Real>::Scalar(std::string filename) : _cellId(0), _inf(std::numeric_limits<Real>::infinity()) {
        auto cgmsh = std::make_unique<hertz::msh::CGmsh<Real>>();
        cgmsh->load(filename);

        hertz::fvm::CGridCreate<Real> gridCreate;
        gridCreate.process(*cgmsh);

        cgmsh.reset();

        _grid = gridCreate.grid();

        _basisP3.resize(10);
        _basisP3_0.resize(10);
    }

    template<typename Real>
    void Scalar<Real>::setFluxFunctions(Real (*f)(const Real &), Real (*df)(const Real &)) {
        _f = f;
        _df = df;
    }

    template<typename Real>
    void Scalar<Real>::initVar() {
        /* Construct the basis function */
        for (Integer i = 0; i < _grid.nbCells(); ++i) {
            _cellId = _grid.cellId(i);

            _cellPtr = &_grid.cell(_cellId);
            _dCells[_cellId] = DCell<Real, Real>();

            _dCellPtr = &_dCells[_cellId];

            if (_cellPtr->nbStencil() > 0) {
//                _basisP3_0.setZero();
//
//                for (Integer nq = 0; nq < _grid.nbCellQuadrature(); ++nq) {
//                    _basisP3_0 += _grid.cellQuadratureWeights()[nq] * _basis.p3(_cellPtr->quadrature(nq));
//                }

                for (Integer iStn = 0; iStn < _cellPtr->nbStencil(); ++iStn) {
                    if (iStn == _cellPtr->nbStencil() - 1) {
//                        std::cout << _cellId << " : ";

                        for (Integer iCell = 0; iCell < _cellPtr->nbCellInStencil(iStn); ++iCell) {
                            _basisP3.setZero();
                            _anotherCellPtr = &_grid.cell(_cellPtr->stencil(iStn)[iCell]);
//                            std::cout << _cellPtr->stencil(iStn)[iCell] << " ";

                            for (Integer nq = 0; nq < _grid.nbCellQuadrature(); ++nq) {
                                _basisP3 += _grid.cellQuadratureWeights()[nq] *
                                            _basis.p3(_cellPtr->scale(_anotherCellPtr->quadrature(nq)));
                            }

                            _dCellPtr->AP3.row(iCell) = _basisP3;
                        }
//                        std::cout << "\n";

                        _dCellPtr->AP3 = _basis.LeastSquare(_dCellPtr->AP3);
                    }

                    _dCellPtr->calculate = true;
                }

                for (Integer f = 0; f < _cellPtr->nbFace(); ++f) {
                    _facePtr = &_grid.face(_cellPtr->faceId(f));
                    _dFacePtr = &_dFaces[_cellPtr->faceId(f)];

                    for (Integer qn = 0; qn < _grid.nbFaceQuadrature(); ++qn) {
                        _basisP3 = _basis.p3(_cellPtr->scale(_facePtr->quadrature(qn)));
                        _dFacePtr->quadrature.push_back(_basisP3);
                    }
                }
            }
        }
    }

    template<typename Real>
    void Scalar<Real>::initialize(Real (*u0)(const CCoordinate<Real> &), Real CFL, Real finalTime) {
        _u0 = u0;
        _cfl = CFL;
        _finalTime = finalTime;

        for (Integer i = 0; i < _grid.nbInternalCell(); ++i) {
            _cellId = _grid.cellId(i);
            _cellPtr = &_grid.cell(_cellId);
            _dCellPtr = &_dCells[_cellId];

            _dCellPtr->u = _u0(_cellPtr->barycenter());
        }

        boundary();
    }

    template<typename Real>
    void Scalar<Real>::boundary() {
        for (Integer i = _grid.nbInternalCell(); i < _grid.nbCells(); ++i) {
            _cellId = _grid.cellId(i);

            if (_grid.cell(_cellId).boundaryType() == BoundaryType::BT_PERIODIC) {
                _dCells[_cellId].u = _dCells[_grid.cell(_cellId).linkedCellId()].u;
            } else {
                std::cout << _cellId << " : " << "Unsupported boundary type\n";
            }
        }
    }

    template<typename Real>
    void Scalar<Real>::construction() {
        VectorXr UP3(10);
//        Real test;
        boundary();

        for (Integer i = 0; i < _grid.nbCells(); ++i) {
            _cellId = _grid.cellId(i);

            _cellPtr = &_grid.cell(_cellId);
            _dCellPtr = &_dCells[_cellId];

            if (_dCellPtr->calculate) {
                for (Integer iStn = 0; iStn < _cellPtr->nbStencil(); ++iStn) {
                    if (iStn == _cellPtr->nbStencil() - 1) {

                        for (Integer iCell = 0; iCell < _cellPtr->nbCellInStencil(iStn); ++iCell) {
                            _dCellNeighborPtr = &_dCells[_cellPtr->stencil(iStn)[iCell]];
                            UP3(iCell) = _dCellNeighborPtr->u;
                        }

                        _dCellPtr->CP3 = _dCellPtr->AP3 * UP3;

//                        test = 0.0;
//                        for(Integer z = 0; z < _cellPtr->nbQuadrature(); ++z) {
//                            test += _grid.cellQuadratureWeights()[z] * (_dCellPtr->CP3.template dot(_basis.p3(_cellPtr->scale(_cellPtr->quadrature(z)))));
//                        }
//
//                        std::cout << test - _dCellPtr->u << "\n";

                    }
                }
            }
        }

        Real ul, ur, maxWaveSpeed;
        CVector<Real> maxWaveSpeedVector;

        for (Integer i = 0; i < _grid.nbInternalCell(); ++i) {
            _cellId = _grid.cellId(i);
            _cellPtr = &_grid.cell(_cellId);
            _dCellPtr = &_dCells[_cellId];

            for (Integer f = 0; f < _cellPtr->nbFace(); ++f) {
                _facePtr = &_grid.face(_cellPtr->faceId(f));

                _dFacePtr = &_dFaces[_cellPtr->faceId(f)];
                _dFacePairPtr = &_dFaces[_facePtr->pairFaceId()];

                if (!_dFacePtr->calculated) {
                    _dFacePtr->flux = 0.0;
                    _dFacePtr->maxWaveSpeed = -_inf;
                    _dCellNeighborPtr = &_dCells.at(_facePtr->neighborId());

                    for (Integer nq = 0; nq < _grid.nbFaceQuadrature(); ++nq) {
                        ul = _dCellPtr->CP3.template dot(_dFacePtr->quadrature[nq]);
                        ur = _dCellNeighborPtr->CP3.template dot(_dFacePtr->quadrature[nq]);

//                        ul = _dCellPtr->u;
//                        ur = _dCellNeighborPtr->u;

                        maxWaveSpeed = fmax(fabs(_df(ul)), fabs(_df(ur)));
                        maxWaveSpeedVector = CVector<Real>(maxWaveSpeed, maxWaveSpeed);

                        maxWaveSpeed = fabs(maxWaveSpeedVector.derived().template dot(_facePtr->normal()));

                        _dFacePtr->maxWaveSpeed = fmax(_dFacePtr->maxWaveSpeed, maxWaveSpeed);

                        _dFacePtr->flux += _grid.faceQuadratureWeights()[nq] *
                                hertz::physics::Scalar<Real>::Rusanov(ul, ur, maxWaveSpeed, _f, _f, _facePtr->normal());

                    }

                    _dFacePtr->flux *= _facePtr->area();
                    _dFacePtr->calculated = true;

                    _dFacePairPtr->maxWaveSpeed = _dFacePtr->maxWaveSpeed;
                    _dFacePairPtr->flux = -_dFacePtr->flux;
                    _dFacePairPtr->calculated = true;
                }

            }
        }

        for (Integer i = 0; i < _grid.nbInternalCell(); ++i) {
            _cellId = _grid.cellId(i);
            _cellPtr = &_grid.cell(_cellId);
            _dCellPtr = &_dCells[_cellId];

            _dCellPtr->maxWaveSpeed = -_inf;
            _dCellPtr->residual = 0.0;

            for (Integer f = 0; f < _cellPtr->nbFace(); ++f) {
                _dFacePtr = &_dFaces[_cellPtr->faceId(f)];

                _dCellPtr->maxWaveSpeed = fmax(_dCellPtr->maxWaveSpeed, _dFacePtr->maxWaveSpeed);
                _dCellPtr->residual -= _dFacePtr->flux;

                _dFacePtr->calculated = false;
            }

            _dCellPtr->residual /= _cellPtr->volume();
        }

    }

    template<typename Real>
    void Scalar<Real>::timeStep() {
        Real minTimeStep = _inf;

        for (Integer i = 0; i < _grid.nbInternalCell(); ++i) {
            _cellId = _grid.cellId(i);
            _cellPtr = &_grid.cell(_cellId);
            _dCellPtr = &_dCells[_cellId];

            minTimeStep = fmin(minTimeStep, _cellPtr->radii() / _dCellPtr->maxWaveSpeed);

        }

        _dt = _cfl * minTimeStep;
        _dt = _t + _dt > _finalTime ? (_finalTime - _t) : _dt;
    }

    template<typename Real>
    void Scalar<Real>::timeStepping() {

        while (_t < _finalTime) {

            construction();
            timeStep();

            for (Integer i = 0; i < _grid.nbInternalCell(); ++i) {
                _cellId = _grid.cellId(i);
                _dCellPtr = &_dCells[_cellId];

                _dCellPtr->urk[0] = _dCellPtr->u;
                _dCellPtr->urk[1] = RK<Real, Real>::RK_3_3(_dCellPtr->urk[0], _dCellPtr->urk[0],
                                                           _dCellPtr->residual, _dt, 1);

                _dCellPtr->u = _dCellPtr->urk[1];
            }

            construction();

            for (Integer i = 0; i < _grid.nbInternalCell(); ++i) {
                _cellId = _grid.cellId(i);
                _dCellPtr = &_dCells[_cellId];

                _dCellPtr->urk[2] = RK<Real, Real>::RK_3_3(_dCellPtr->urk[0], _dCellPtr->urk[1],
                                                           _dCellPtr->residual, _dt, 2);
                _dCellPtr->u = _dCellPtr->urk[2];
            }

            construction();

            for (Integer i = 0; i < _grid.nbInternalCell(); ++i) {
                _cellId = _grid.cellId(i);
                _dCellPtr = &_dCells[_cellId];

                _dCellPtr->u = RK<Real, Real>::RK_3_3(_dCellPtr->urk[0], _dCellPtr->urk[2],
                                                      _dCellPtr->residual, _dt, 3);
            }

            _t += _dt;
            std::cout << _t << " << time\n";
        }

    }

    template<typename Real>
    void Scalar<Real>::write() {
        const char *filename = R"(../user/out.dat)";
        std::ofstream out;

        out.open(filename);

        Real L1 = 0.0;
        CCoordinate<Real> coordinate;

        if (out.good() && out.is_open()) {
            goto lanjutmang;
        } else {
            out.close();
            exit(-1);
        }

        lanjutmang:
        for (Integer i = 0; i < _grid.nbCells(); ++i) {
            _cellId = _grid.cellId(i);
            _cellPtr = &_grid.cell(_cellId);
            _dCellPtr = &_dCells[_cellId];

            coordinate = (_cellPtr->barycenter().array() - _finalTime).matrix();

            L1 += fabs(_dCellPtr->u - _u0(coordinate));
            out << _grid.cell(_cellId).barycenter().transpose() << "\t" << _dCellPtr->u << "\t" << _cellId << "\n";
        }

        std::cout << "L1 Error = " << L1 / _grid.nbInternalCell() << " - N = " << _grid.nbInternalCell() << "\n";

        out.close();

        std::string fn = "../src/py/plot.py";
        std::string command = "python3 ";
        command += fn;

        FILE* in = popen(command.c_str(), "r");
        pclose(in);
    }

    template<typename Real>
    Scalar<Real>::~Scalar() = default;

}

#endif //HERTZ_SOLVER_SCALAR_H
