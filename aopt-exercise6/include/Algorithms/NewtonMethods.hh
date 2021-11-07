#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>
#include "LineSearch.hh"

//== NAMESPACES ===============================================================

namespace AOPT {
    /**
    * @brief NewtonMethods is just a list of functions implementing several variations of the
    * newton's method */
    class NewtonMethods {
    public:
        typedef FunctionBaseSparse::Vec Vec;   // dense vector arbitrary size
        typedef FunctionBaseSparse::Mat Mat;   // dense matrix arbitrary size
        typedef FunctionBaseSparse::T T;        //Triplets
        typedef FunctionBaseSparse::SMat SMat;  // sparse matrix arbitrary size

        /**
         * @brief solve
         * \param _problem pointer to any function/problem inheriting from FunctionBaseSparse
         *        on which the basic Newton Method will be applied
         * \param _initial_x starting point of the method
         * \param _eps epsilon under which the method stops
         * \param _max_iters maximum iteration of the method*/
        static Vec solve(FunctionBaseSparse *_problem, const Vec &_initial_x, const double _eps = 1e-4,
                         const int _max_iters = 1000000) {
            std::cout << "******** Newton Method ********" << std::endl;

            // squared epsilon for stopping criterion
            const double e2 = 2 * _eps;
            const double ee = _eps * _eps;

            const int n = _problem->n_unknowns();

            // get starting point
            Vec x = _initial_x;

            // allocate gradient storage
            Vec g(n);

            // allocate hessian storage
            SMat H(n, n);

            // allocate search direction vector storage
            Vec delta_x(n);

            Eigen::SimplicialLLT<SMat> solver;

            //------------------------------------------------------//
            // implement Newton method

            auto fp = std::numeric_limits<double>::max();

            // /!\ don't compute the inverse of the hessian matrix
            //      --> use cholesky factorization
            //      --> solve left hand side
            //      --> check LDLT module of the Eigen library

            // repeat:
            for (size_t i = 0; i < _max_iters; ++i) {
                // 1. compute newton step: delta_x = -H^-1 * g
                //      /!\ H is costly to compute. Use cholesky factorization
                //      --> H = L * L^T
                _problem->eval_hessian(x, H);
                solver.compute(H);

                // solver.solve(b)
                //      returns the solution x of A x = b using the current decomposition of A.
                //      --> A being the hessian matrix we did plug in the `compute` method
                //      --> x being delta_x
                _problem->eval_gradient(x, g);
                delta_x = solver.solve(-g);

                // 2. compute newton decrement lambda^2 = -g^T * delta_x
                const auto lambda_2 = -g.transpose() * delta_x;

                // 3. Stopping criterion: quit if lambda^2/2 <= epsilon
                if (lambda_2 <= e2)
                    break;

                // Additional checks like given in the instructions
                const auto f = _problem->eval_f(x);
                if (
                        // * objective function does not increase
                        f >= fp ||
                        // * norm of gradient is not greater than epsilon
                        g.squaredNorm() <= ee
                )
                    break;

                // 4. Line search: choose step size: t > 0 (e.g. backtracking t_0 = 1)
                // use the alpha and tau constant as given in the slides
                const auto t_k = LineSearch::backtracking_line_search(_problem, x, g, delta_x, 1., .5, .75);

                // 5. Update: x = x + t * delta_x
                x += t_k * delta_x;

                // updates for stopping criteria
                fp = f;
            }

            //------------------------------------------------------//

            return x;
        }

        /**
         * @brief solve with the Projected Hessian method
         * \param _problem pointer to any function/problem inheriting from FunctionBaseSparse.
         *        This problem MUST provide a working eval_hession() function for this method to work.
         *
         * \param _initial_x starting point of the method
         * \param _tau_factor the evolution factor of the tau coefficient
         * \param _eps epsilon under which the method stops
         * \param _max_iters maximum iteration of the method*/
        static Vec
        solve_with_projected_hessian(FunctionBaseSparse *_problem, const Vec &_initial_x, const double _gamma = 10.0,
                                     const double _eps = 1e-4, const int _max_iters = 1000000) {
            bool converged = false;
            return solve_with_projected_hessian(_problem, converged, _initial_x, _gamma, _eps, _max_iters);
        }

        static Vec solve_with_projected_hessian(FunctionBaseSparse *_problem, bool &_converged, const Vec &_initial_x,
                                                const double _gamma = 10.0,
                                                const double _eps = 1e-4, const int _max_iters = 1000000) {
            std::cout << "******** Newton Method with projected hessian ********" << std::endl;

            // squared epsilon for stopping criterion
            double e2 = 2 * _eps * _eps;

            int n = _problem->n_unknowns();

            // get starting point
            Vec x = _initial_x;

            // allocate gradient storage
            Vec g(n);

            // allocate hessian storage
            SMat H(n, n);

            // allocate search direction vector storage
            Vec delta_x(n);
            int iter(0);

            // identity and scalar to add positive values to the diagonal
            SMat I(n, n);
            I.setIdentity();

            _converged = false;

            Eigen::SimplicialLLT<SMat> solver;

            //------------------------------------------------------//
            //TODO: implement Newton with projected hessian method
            //Hint: if the factorization fails, then add delta * I to the hessian.
            //      repeat until factorization succeeds (make sure to update delta!)

            //------------------------------------------------------//


            return x;
        }


    };

} // namespace AOPT
