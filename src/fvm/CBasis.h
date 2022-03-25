//
// Created by indro on 10/03/21.
//

#ifndef HERTZ_CBASIS_H
#define HERTZ_CBASIS_H

namespace hertz::fvm {
    using Eigen::Matrix;
    using Eigen::Dynamic;
    using hertz::geo::CCoordinate;

    template<typename Real>
    class CBasis {
    private:
        Matrix<Real, 10, 1> _p3;
        Matrix<Real, 5, 1> _p2;
        Matrix<Real, 3, 1> _p1;

    public:
        typedef Matrix<Real, Dynamic, Dynamic> MatrixXr;

        CBasis();

        MatrixXr LeastSquare(const MatrixXr &sysMat);

        const Eigen::Matrix<Real, 10, 1> &p3(const CCoordinate<Real> &coordinate);

        const Eigen::Matrix<Real, 5, 1> &p2(const CCoordinate<Real> &coordinate);

        const Eigen::Matrix<Real, 3, 1> &p1(const CCoordinate<Real> &coordinate);

    };

    /*============================ Implementation =======================*/

    template<typename Real>
    CBasis<Real>::CBasis() {
        _p1.setZero();
        _p3.setZero();
    }

    template<typename Real>
    typename CBasis<Real>::MatrixXr CBasis<Real>::LeastSquare(const CBasis::MatrixXr &sysMat) {
        return (sysMat.transpose() * sysMat).inverse() * sysMat.transpose();
    }

    template<typename Real>
    const Eigen::Matrix<Real, 10, 1> &CBasis<Real>::p3(const CCoordinate<Real> &coordinate) {
        _p3(9) = 1.0;
        _p3(8) = coordinate.y();
        _p3(7) = coordinate.x();
        _p3(6) = coordinate.y() * coordinate.y();
        _p3(5) = coordinate.x() * coordinate.y();
        _p3(4) = coordinate.x() * coordinate.x();
        _p3(3) = _p3(5) * coordinate.y();
        _p3(2) = _p3(5) * coordinate.x();
        _p3(1) = _p3(6) * coordinate.y();
        _p3(0) = _p3(4) * coordinate.x();

        return _p3;
    }

    template<typename Real>
    const Eigen::Matrix<Real, 3, 1> &CBasis<Real>::p1(const CCoordinate<Real> &coordinate) {
        _p1(2) = 1.0;
        _p1(1) = coordinate.x();
        _p1(0) = coordinate.y();

        return _p1;
    }

    template<typename Real>
    const Matrix<Real, 5, 1> &CBasis<Real>::p2(const CCoordinate<Real> &coordinate) {
        _p2(4) = coordinate.y();
        _p2(3) = coordinate.x();
        _p2(2) = coordinate.y() * coordinate.y();
        _p2(1) = coordinate.x() * coordinate.y();
        _p2(0) = coordinate.x() * coordinate.x();

        return _p2;
    }
}

#endif //HERTZ_CBASIS_H
