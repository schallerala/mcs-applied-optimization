#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>
#include <Functions/AugmentedLagrangianProblem.hh>
#include <Utils/OptimizationStatistic.hh>
#include <Algorithms/NewtonMethods.hh>
#include "LBFGS.hh"


//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================
    class AugmentedLagrangian {
    public:
        // LA typedefs
        typedef FunctionBaseSparse::Vec Vec;

        static Vec
        solve(FunctionBaseSparse *_obj, const Vec &_initial_x, const std::vector<FunctionBaseSparse *> &_constraints,
              const std::vector<FunctionBaseSparse *> &_squared_constraints,
                // tolerances
              const double _eta = 1e-4, const double _tau = 1e-4,
              const int _max_iters = 20) {
            std::cout << "******** Augmented Lagrangian ********" << std::endl;

            // Initialize
            double mu = 10,
                    tau = 1. / mu,
                    eta = std::pow(mu, -0.1),
                    h_sqr_norm = 0.,
                    h_sqr_norm_previous = std::numeric_limits<double>::max();

            const auto tau2 = _tau * _tau;

            //vector of nu and vector of constraint value
            Vec nu(_constraints.size()), h(_constraints.size());
            nu.setZero(); // said in the video that it is fine to start with vector 0
            h.setZero();

            //initialize the augmented lagrangian problem for the unconstrained solver
            AugmentedLagrangianProblem problem(_obj, _constraints, _squared_constraints, nu, mu);
            auto opt_st = std::make_unique<AOPT::OptimizationStatistic>(&problem);

            auto stats_problem = opt_st.get();

            //get starting point
            Vec x = _initial_x;
            //store previous point
            Vec x_p = x;

            //allocate gradient storage
            Vec g(problem.n_unknowns());

            //------------------------------------------------------//
            // implement the augmented lagrangian method.
            // Hints: 1. Use projected newton method to solve for an approximated x.
            //           If the maximum iteration is reached or if the norm of the constraints
            //           gets larger, one can say it diverges for simplicity.
            //        2. Use set_mu(...) and set_nu(...) functions in AugmentedLagrangianProblem
            //           class to apply the change of nu and mu

            // References:
            // Algorithm: Slide 21 of Lecture week 10
            //      https://slides.cgg.unibe.ch/aopt21/10-EqualityConstrainedOptimization-II-deck.html#/augmented-lagrangian-method-1
            // Formulas: Slide 18 of Lecture week 10
            //      https://slides.cgg.unibe.ch/aopt21/10-EqualityConstrainedOptimization-II-deck.html#/17/0/7

            for (size_t i = 0; i < _max_iters; ++i) {
                // (update the constants to the underlying problem: hint 2)
                problem.set_nu(nu);
                problem.set_mu(mu);

                // Minimize sub-problem with fixed ν^k and μ_k
                //      Find approximate solution x_{k+1} of L_A(x, ν^k, μ_k) starting at
                //      x_k and stopping when norm(∇_x L_A(x_{k+1}, ν_k, μ_k)) <= τ_k
                //      (x_{k+1} = x_k if diverging)


                // ? Should iterate?
                //      With the text above "find approximate solution [...] and stopping when [...]"

                bool newton_converged;
                // x_{k+1} with projected newton method (as hint 1)
                x = NewtonMethods::solve_with_projected_hessian(stats_problem, newton_converged, x); // optimize like
                                                                                                        // for none-
                                                                                                        // constraints
                                                                                                        // problem like
                                                                                                        // said during
                                                                                                        // lecture 10 at
                                                                                                        // video record-
                                                                                                        // ing.
                                                                                                        // timestamp:
                                                                                                        // 1h03m40s

                // Remains of hint 1:
                // if
                //      A. maximum iteration is reached
                //      B. or if the norm of the constraints get larger:
                // then
                //      one can say it diverges for simplicity

                // A. (hint 1)
                if ( ! newton_converged)
                    break;

                // constraint violation:
                //  h_i(x_k) ~= (ν_i^* - ν_i^k) / μ_k
                problem.eval_constraints(x, h);
                h_sqr_norm = h.squaredNorm();

                // B. (hint 1)
                if (h_sqr_norm > h_sqr_norm_previous)
                    break;

                // update previous squared norm of constraints violations
                h_sqr_norm_previous = h_sqr_norm;

                // Augmented Lagrangian Function:
                //  L_A(x, ν, μ) = f(x) + sum_{i=1}^p ν_i * h_i(x) + μ/2 sum_{i=1}^p h_i^2(x)

                stats_problem->eval_gradient(x, g);
                const auto g_lagrangian_sqr_norm = g.squaredNorm();

                // if norm(h(x_{k+1})) <= η_k
                // then
                if (h_sqr_norm <= eta) {
                    // test for convergence

                    // if norm(h(x_{k+1})) <= η^* and norm(∇_x L_A(x_{k+1}, ν_k, μ_k)) <= τ^*
                    // then
                    if (h_sqr_norm <= _eta && g_lagrangian_sqr_norm <= tau2) {
                        // stop with approximate solution x_k
                        break;
                    }

                    // Update Multipliers, Tighten Tolerances
                    nu += mu * h; // ν^{k+1} = ν^k + μ_k * h(x_{k+1})
//                    mu = mu; // μ_{k+1} = μ_k
                    eta /= std::pow(mu, 0.9); // η_{k+1} = η_k / μ_{k+1}^0.9
                    tau /= mu; // τ_{k+1} = τ_k / μ_{k+1}
                }
                // norm(h(x_{k+1})) > η_k
                else {
                    // increase penalty parameter, reset tolerances
//                    nu = nu; // ν^{k+1} = ν^k
                    mu *= 100; // μ_{k+1} = 100*μ_k
                    eta = 1 / std::pow(mu, 0.1); // η_{k+1} = 1 / μ_{k+1}^0.1
                    tau = 1 / mu; // τ_{k+1} = 1 / μ_{k+1}
                }
            }

            //------------------------------------------------------//

            opt_st->print_statistics();

            return x;
        }
    };
    //=============================================================================

}



