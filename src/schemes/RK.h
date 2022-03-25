//
// Created by indro on 10/03/21.
//

#ifndef HERTZ_RK_H
#define HERTZ_RK_H

namespace hertz::schemes {

    template<typename Real, typename Derived>
    class RK {
    public:
        static Derived RK_3_3(const Derived &q0, const Derived &qn, const Derived &residual, const Real &dt,
                              const int &stage) {

            if (stage == 1) {
                return q0 + dt * residual;
            } else if (stage == 2) {
                return 0.75 * q0 + 0.25 * (qn + residual * dt);
            } else {
                return (q0 + 2.0 * (qn + dt * residual)) / 3.0;
            }
        }
    };

}

#endif //HERTZ_RK_H
