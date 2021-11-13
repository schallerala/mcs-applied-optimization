#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>

//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================

    class LineSearch {
    public:
        typedef FunctionBaseSparse::Vec Vec;
        typedef FunctionBaseSparse::SMat SMat;

        /** Back-tracking line search method
         *
         * \param _problem a pointer to a specific Problem, which can be any type that
         *        has the same interface as FunctionBase's (i.e. with eval_f, eval_gradient, etc.)
         * \param _x starting point of the method. Should be of the same dimension as the Problem's
         * \param _g gradient at the starting point.
         * \param _dx delta x
         * \param _t0 inital step of the method
         * \param _alpha and _tau variation constant, as stated by the method's definition
         * \return the final step t computed by the back-tracking line search */
        template<class Problem>
        static double backtracking_line_search(Problem *_problem,
                                               const Vec &_x,
                                               const Vec &_g,
                                               const Vec &_dx,
                                               const double _t0,
                                               const double _alpha = 0.5,
                                               const double _tau = 0.75) {

            double t(0);

            //------------------------------------------------------//
            t = _t0;

            // pre-compute objective
            double fx = _problem->eval_f(_x);

            // pre-compute dot product
            double gtdx = _g.transpose() * _dx;

            // make sure dx points to a descent direction
            if (gtdx > 0) {
                std::cerr << "dx is in the direction that increases the function value. gTdx = " << gtdx << std::endl;
                return t;
            }

            // backtracking (stable in case of NAN)
            int i = 0;
            while (!(_problem->eval_f(_x + t * _dx) <= fx + _alpha * t * gtdx) && i < 1000) {
                t *= _tau;
                i++;
            }

            //------------------------------------------------------//

            return t;
        }


        /** Back-tracking line search for infeasible start Newton's method
        *
        * \param _problem a pointer to a specific Problem, which can be any type that
        *        has the same interface as FunctionBase's (i.e. with eval_f, eval_gradient, etc.)
        * \param _A matrix of the lhs of linear equality constraints
        * \param _b vector of the rhs of linear equality constraints
        * \param _x primal starting point of the method
        * \param _nu dual starting point of the method
        * \param _dx delta x
        * \param _dnu delta nu
        * \param _initial_res initial residual
        * \param _t0 inital step of the method
        * \param _alpha and _beta variation constant, as stated by the method's definition
        * \return the final step t computed by the back-tracking line search */
        template<class Problem>
        static double backtracking_line_search_newton_with_infeasible_start(Problem *_problem,
                                                                            const SMat &_A,
                                                                            const Vec &_b,
                                                                            const Vec &_x,
                                                                            const Vec &_nu,
                                                                            const Vec &_dx,
                                                                            const Vec &_dnu,
                                                                            const double _initial_res,
                                                                            const double _t0,
                                                                            const double _alpha = 0.01,
                                                                            const double _beta = 0.9) {
            //------------------------------------------------------//
            double t = _t0;

            //TODO: implement the algorithm
            int n = _problem->n_unknowns();
            Vec ga(n), Anua(n), rpr_a(n), rdual_a(n);
            double resa(0);

            // backtracking
            while (1) {
                _problem->eval_gradient(_x + t * _dx, ga);
                Anua = _A.transpose() * (_nu + t * _dnu);

                rpr_a = _A * (_x + t * _dx) - _b;
                rdual_a = ga + Anua;

                resa = rpr_a.squaredNorm() + rdual_a.squaredNorm();
                resa = sqrt(resa);

                if (resa <= (1. - _alpha * t) * _initial_res)
                    break;

                if (t < 1e-7) {
                    return t;
                }
                t *= _beta;
            }

            //------------------------------------------------------//

            return t;
        }


        template<class Problem>
        static double wolfe_line_search(Problem *_problem,
                                        const Vec &_0,
                                        const Vec &_g,
                                        const Vec &_d0,
                                        const double _alpha_0,
                                        const double _alpha_max = 100,
                                        const double _c_1 = 1e-4,
                                        const double _c_2 = 0.9) {
            //------------------------------------------------------//
            assert(0 < _c_1 && _c_1 < 1);
            assert(0 < _c_2 && _c_2 < 1);
            assert(_c_1 < _c_2);
            // choose alpha_{max} > 0
            assert(_alpha_max > 0);

            double alpha_i = _alpha_0;

            // choose alpha_1 \in (0, alpha_{max})
            assert(0 <= alpha_i && alpha_i <= _alpha_max);

            // increase rate
            const double inc = 2.;

            // save the function value at the alpha_i = 0
            const auto phi_0 = _problem->eval_f(_0);

            // projection of gradient on the search direction
            const auto dphi_0 = _g.dot(_d0);

            // make sure dx points to a descent direction
            if (dphi_0 > 0) {
                std::cerr << "dx is in the direction that increases the function value." << std::endl;
                return alpha_i;
            }

            // save constant across iteration related to _c_1 and _c_2 parameters
            const auto dphi_test = _c_1 * dphi_0;
            // used as stopping criterion
            const auto dphi_wolfe = - _c_2 * dphi_0;

            // first stage:
            // begins with a trial estimate alpha_i, and keeps increasing it until it finds either
            // an acceptable step length or an interval that brackets the desired step lengths
            Vec g(_problem->n_unknowns());

            // alpha_0 = 0
            double alpha_p = 0;
            // previous evaluation of phi(alpha_{i - 1})
            double eval_alpha_p = phi_0;
            // previous evaluation of phi'(x)
            double dphi_p = dphi_0;

            // i = 1
            for (size_t i = 1; i <= 20; ++i) {
                // Evaluate phi(alpha_i)
                double eval_alpha_i = _problem->eval_f(_0 + alpha_i * _d0);

                _problem->eval_gradient(_0 + alpha_i * _d0, g);
                const double dphi = g.dot(_d0);

                // if phi(alpha_i) > phi(0) + _c_1 * alpha_i * phi'(0) or (phi(alpha_i) >= phi(alpha_{i - 1}) and i > 1)
                if (eval_alpha_i > phi_0 + alpha_i * dphi_test || (i > 1 && eval_alpha_i >= eval_alpha_p)) {
                    // alpha_* = zoom(alpha_{i - 1}, alpha_i)
                    alpha_i = zoom(_problem, _0, _d0, phi_0, dphi_0, eval_alpha_p, eval_alpha_i, dphi_p, dphi,
                                   alpha_p, alpha_i);
                    // and stop
                    return alpha_i;
                }

                // Evaluate dphi = phi'(alpha_i) <--    already done as also used in the zoom in case we stop in the
                //                                      previous condition

                // if abs(phi'(alpha_i)) <= - _c_2 * phi'(0)
                if (std::abs(dphi) <= dphi_wolfe)
                    // alpha_* = alpha_i
                    // and stop
                    return alpha_i;

                // if phi'(alpha_i) >= 0
                if (dphi >= 0) {
                    // alpha_* = zoom(alpha_i, alpha_{i-1})
                    alpha_i = zoom(_problem, _0, _d0, phi_0, dphi_0, eval_alpha_i, eval_alpha_p, dphi, dphi_p,
                                   alpha_i, alpha_p);
                    // and stop
                    return alpha_i;
                }

                alpha_p = alpha_i;
                eval_alpha_p = eval_alpha_i;
                dphi_p = dphi;

                // Choose alpha_{i+1} \in (alpha_i, alpha_{max})

                // increase alpha_i by 2 in (alpha_i, t_max)
                if (alpha_i * inc > _alpha_max) {
                    alpha_i += _alpha_max;
                    alpha_i /= 2.;
                } else
                    alpha_i *= inc;
            }

            //------------------------------------------------------//

            return alpha_i;
        }


    private:
        template<class Problem>
        static double zoom(Problem *_problem,
                           const Vec &_0,
                           const Vec &_d0,
                           const double _phi_0,
                           double _dphi_0,
                           double _eval_alpha_lo,
                           double _eval_alpha_hi,
                           double _dphi_lo,
                           double _dphi_hi, // never used, but to avoid changing the signature, keep it
                           double _alpha_lo,
                           double _alpha_hi,
                           const double _c_1 = 1e-4,
                           const double _c_2 = 0.9) {
            // second stage:
            // successively decreases the size of the interval until
            // an acceptable step length is identified.
            // save constant across iteration related to _c_1 and _c_2 parameters
            const auto dphi_test = _c_1 * _dphi_0;
            // used as stopping criterion
            const auto dphi_wolfe = - _c_2 * _dphi_0;

            Vec g(_problem->n_unknowns());
            double alpha_j(1.);

            double alpha_hi;
            double alpha_lo;
//            double dphi_hi; <-- never used
            double dphi_lo;

            for (size_t i = 0; i <= 20; ++i) {
                if (i == 0) {
                    alpha_hi = _eval_alpha_hi;
//                    dphi_hi = _dphi_hi;
                    alpha_lo = _eval_alpha_lo;
                    dphi_lo = _dphi_lo;
                } else {
                    alpha_hi = _problem->eval_f(_0 + _alpha_hi * _d0);
                    alpha_lo = _problem->eval_f(_0 + _alpha_lo * _d0);
                    _problem->eval_gradient(_0 + _alpha_lo * _d0, g);
                    dphi_lo = g.dot(_d0);
                }


                // use {alpha_lo, alpha_hi, dphi_lo} to make a quadric interpolation of
                // the function said interpolation is used to estimate the minimum
                //
                // polynomial:    p(x) = a*(x - _t)^2 + b
                // conditions:  p(thi) = alpha_hi
                //              p(tlo) = alpha_lo
                //             p'(tlo) = dphi_lo
                alpha_j = (alpha_hi - alpha_lo) * _alpha_lo - (_alpha_hi * _alpha_hi - _alpha_lo * _alpha_lo) * dphi_lo / 2;
                double div = (alpha_hi - alpha_lo) - (_alpha_hi - _alpha_lo) * dphi_lo;
                if (div == 0)
                    alpha_j = _alpha_lo;
                else
                    alpha_j /= div;

                // if interpolation fails, bisection is used
                if (alpha_j <= std::min(_alpha_lo, _alpha_hi) || alpha_j >= std::max(_alpha_lo, _alpha_hi))
                    alpha_j = (_alpha_lo + _alpha_hi) / 2;

                // Evaluate phi(alpha_j)
                double eval_alpha_j = _problem->eval_f(_0 + alpha_j * _d0);
                _problem->eval_gradient(_0 + alpha_j * _d0, g);

                const double dphi = g.dot(_d0);

                // if phi(alpha_j) > phi(0) + c_1 * alpha_j * phi'(0) or phi(alpha_j) >= phi(alpha_{lo})
                if (eval_alpha_j > _phi_0 + alpha_j * dphi_test || eval_alpha_j >= alpha_lo) {

                    if (alpha_j == _alpha_hi) {
                        std::cerr << "alpha_j equals to thi, possibly due to insufficient numeric precision." << std::endl;
                        return alpha_j;
                    }

                    // alpha_{hi} = alpha_j
                    _alpha_hi = alpha_j;
                    alpha_hi = eval_alpha_j;
//                    dphi_hi = dphi;
                }
                // else
                else {
                    // Evaluate phi'(alpha_j)
                    // if abs(phi'(alpha_j)) <= - c_2 * phi'(0)
                    if (std::abs(dphi) <= dphi_wolfe)
                        // alpha_* = alpha_j
                        // and stop
                        return alpha_j;

                    // if phi'(alpha_j) * (alpha_{hi} - alpha_{lo}) >= 0
                    if (dphi * (_alpha_hi - _alpha_lo) >= 0) {
                        _alpha_hi = _alpha_lo;
                        // alpha_{hi} = alpha_{lo}
                        alpha_hi = alpha_lo;
//                        dphi_hi = dphi_lo;
                    }

                    if (alpha_j == _alpha_lo) {
                        std::cerr << "alpha_j equals to tlo, possibly due to insufficient numeric precision." << std::endl;
                        return alpha_j;
                    }


                    // alpha_{lo} = alpha_j
                    _alpha_lo = alpha_j;
                    alpha_lo = eval_alpha_j;
//                    dphi_lo = dphi;
                }
            }

            return alpha_j;
        }
    };
//=============================================================================
}



