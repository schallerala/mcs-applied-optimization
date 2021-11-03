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
         * \param _t0 initial step of the method
         * \param _alpha and _tau variation constant, as stated by the method's definition
         * \return the final step t computed by the back-tracking line search */
        template<class Problem>
        static double backtracking_line_search(Problem *_problem,
                                               const Vec &_x,
                                               const Vec &_g,
                                               const Vec &_dx,
                                               const double _t0,
                                               const double _alpha = 0.2,
                                               const double _tau = 0.6) {

            double t = _t0;

            //------------------------------------------------------//
            // implement the backtracking line search algorithm

            const auto f_k = _problem->eval_f(_x);

            while (true) { // bit scary, but don't be :)
                const auto lhs = _problem->eval_f(_x + t * _dx);
                const auto rhs = f_k + t * _alpha * _g.transpose() * _dx;

                // if lhs or rhs is NaN, it will in any case produce false.
                // also in this case, break
                if ( ! (lhs > rhs))
                    break;

                t *= _tau;
            }


            //------------------------------------------------------//

            return t;
        }


    private:
        template<class Problem>
        static double zoom(Problem *_problem,
                           const Vec &_x,
                           const Vec &_dx,
                           double _fx_init,
                           double _dg_init,
                           double _fx_hi,
                           double _fx_lo,
                           double _dg_hi,
                           double _dg_lo,
                           double _thi, double _tlo) {
            // second stage:
            // successively decreases the size of the interval until
            // an acceptable step length is identified.
            const double dg_test = 1e-4 * _dg_init,
                    dg_wolfe = -0.9 * _dg_init;

            int iter(0);
            Vec g(_problem->n_unknowns());
            double t(1.);

            double fx_hi, fx_lo, dg_hi, dg_lo;

            do {
                if (iter == 0) {
                    fx_hi = _fx_hi;
                    dg_hi = _dg_hi;
                    fx_lo = _fx_lo;
                    dg_lo = _dg_lo;
                } else {
                    fx_hi = _problem->eval_f(_x + _thi * _dx);
                    fx_lo = _problem->eval_f(_x + _tlo * _dx);
                    _problem->eval_gradient(_x + _tlo * _dx, g);
                    dg_lo = g.dot(_dx);
                }


                // use {fx_lo, fx_hi, dg_lo} to make a quadric interpolation of
                // the function said interpolation is used to estimate the minimum
                //
                // polynomial: p (x) = a*(x - _t)² + b
                // conditions: p (thi) = fx_hi
                //             p (tlo) = fx_lo
                //             p'(tlo) = dg_lo
                t = (fx_hi - fx_lo) * _tlo - (_thi * _thi - _tlo * _tlo) * dg_lo / 2;
                double div = (fx_hi - fx_lo) - (_thi - _tlo) * dg_lo;
                if (div == 0)
                    t = _tlo;
                else
                    t /= div;

                // if interpolation fails, bisection is used
                if (t <= std::min(_tlo, _thi) || t >= std::max(_tlo, _thi))
                    t = (_tlo + _thi) / 2;

                double fx = _problem->eval_f(_x + t * _dx);
                _problem->eval_gradient(_x + t * _dx, g);

                const double dg = g.dot(_dx);

                if (fx - _fx_init > t * dg_test || fx >= fx_lo) {
                    if (t == _thi) {
                        std::cerr << "t equals to thi, possibly due to insufficient numeric precision." << std::endl;
                        return t;
                    }

                    _thi = t;
                    fx_hi = fx;
                    dg_hi = dg;
                } else {
                    if (std::abs(dg) <= dg_wolfe)
                        return t;

                    if (dg * (_thi - _tlo) >= 0) {
                        _thi = _tlo;
                        fx_hi = fx_lo;
                        dg_hi = dg_lo;
                    }

                    if (t == _tlo) {
                        std::cerr << "t equals to tlo, possibly due to insufficient numeric precision." << std::endl;
                        return t;
                    }

                    _tlo = t;
                    fx_lo = fx;
                    dg_lo = dg;
                }

                iter++;
            } while (iter <= 20);

            return t;
        }
    };
//=============================================================================
}



