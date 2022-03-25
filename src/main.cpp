#include <iostream>
#include "physics/Scalar.h"
#include "solver/Scalar.h"

int main() {
    std::string filename = R"(../user/test10.msh)";

    double cfl = 0.4, finalTime = 0.2;

    hertz::solver::Scalar<double> solver(filename);
    solver.setFluxFunctions(hertz::physics::Scalar<double>::linearAdvection, hertz::physics::Scalar<double>::dLinearAdvection);
    solver.initVar();
    solver.initialize(hertz::physics::Scalar<double>::sine, cfl, finalTime);
    solver.timeStepping();
    solver.write();

    return 0;
}
