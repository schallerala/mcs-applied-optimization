#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>
#include <Functions/InteriorPointProblem.hh>
#include <Algorithms/NewtonMethods.hh>
#include <Utils/OptimizationStatistic.hh>

//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================
    class InteriorPoint {
    public:
        // LA typedefs
        using Vec = FunctionBaseSparse::Vec;

        static Vec
        solve(FunctionBaseSparse *_obj, const Vec &_initial_x, const std::vector<FunctionBaseSparse *> &_constraints,
              const double _eps = 1e-4, const double _mu = 10.0, const int _max_iters = 1000) {
            std::cerr << "******** Interior Point ********" << std::endl;

            // construct log-barrier problem
            InteriorPointProblem problem(_obj, _constraints);
            const auto opt_st = std::make_unique<AOPT::OptimizationStatistic>(&problem);

            // number of constraints
            const double m = _constraints.size();

            // points
            Vec x = _initial_x;
            Vec x_previous = x;

            //------------------------------------------------------//
            // implement the interior point method
            // Hint: Use projected newton method to solve for an approximated x.

            // Reference:
            //      https://slides.cgg.unibe.ch/aopt21/12-InequalityConstrainedOptimization-II-deck.html#/17/0/1

            // Assume feasible x_0

            // t := 1
            double t = 1.0;

            for (size_t i = 0; i < _max_iters; ++i) {
                problem.setT(t);

                bool converged;
                // hint:
                // Minimize f(x) - 1 / t * \sum_{i=1}^m log(-g_i(x))
                x = NewtonMethods::solve_with_projected_hessian(opt_st.get(), converged, x, 10, _eps, 1000);

                if (!converged) {
                    x = x_previous;
                    break;
                }

                // And update t := 10*t
                t *= _mu; // the given mu is actually by default 10

                x_previous = x;

                // Stop when m/t < ε
                if (m / t < _eps)
                    break;
            }

            //------------------------------------------------------//

            opt_st->print_statistics();

            return x;
        }
    };
    //=============================================================================

}





