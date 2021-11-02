#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>
#include "LineSearch.hh"

//== NAMESPACES ===============================================================

namespace AOPT {

    /* Performs a gradient descent on a given problem.
     * This can work with any Problem with a FunctionBase-style interface since
     * the gradient descent method is rather generic mathematically */
    class GradientDescent {
    public:
        typedef FunctionBaseSparse::Vec Vec; ///< Eigen::VectorXd


        /**
         * \param _problem a pointer to a specific Problem, which can be any type that
         *        has the same interface as FunctionBase's (i.e. with eval_f, eval_gradient, etc.)
         * \param _initial_x  the x starting point
         * \param _eps the stopping criterion below which we consider the method
         *             to be done
         * \param _max_iters a capping number of iterations in case you would end-up with a
         *             bad configuration where the successive attempts of finding the
         *             minimum kind of oscillate around the actual minimum without
         *             finding it
         *
         * \return the minimum found by the method. */
        template<class Problem>
        static Vec solve(Problem *_problem, const Vec &_initial_x, const double _eps = 1e-4, const int _max_iters = 1000000) {
            std::cout << "******** Gradient Descent ********" << std::endl;

            // squared epsilon for stopping criterion
            double e2 = _eps * _eps;

            // get starting point
            Vec x = _initial_x;
            Vec x_next;

            // allocate gradient storage
            Vec g(_problem->n_unknowns());
            Vec delta_x_k;

            //------------------------------------------------------//
            // implement the gradient descent

            // repeat
            for (size_t k = 0; k < _max_iters; ++k, x = x_next) {
                // 1. Determine descent direction \Delta x^{(k)}
                //      \Delta x^{(k)} = - \nabla f(x)
                _problem->eval_gradient(x, g);
                delta_x_k = -g;

                // 2. Line Search: choose a step size t^{(k)} > 0
                // TODO find step size? Or step with delta x?

                // 3. Update: x^{(k+1)} = x^{(k)} + t^{(k)} * \Delta x^{(k)}
                x_next = x + delta_x_k;

                // 4. k = k + 1

                // until \nabla f(x^{(k)}) <= epsilon^2
                if (g.norm() <= e2)
                    break;
            }

            //------------------------------------------------------//

            return x;
        }
    };
}



