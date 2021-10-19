#pragma once

#include <FunctionBase/ParametricFunctionBase.hh>

//== NAMESPACES ===============================================================

namespace AOPT {

//== CLASS DEFINITION =========================================================

    /* This function evaluates the energy of an ideal spring with no length
     * going from x_a to x_b.
     * It is a Parametric Function because it requires the elastic constant
     * parameter k_ab for the energy computation. */
    class SpringElement2D : public ParametricFunctionBase {
    public:
        // E_ab(x) = 1/2 * k * ((x[0] - x[2])^2 + (x[1] - x[3])^2)
        // constructor
        SpringElement2D() : ParametricFunctionBase() {}

        // number of unknowns
        inline virtual int n_unknowns() override { return 4; }

        /** evaluates the spring element's energy
         * \param _x contains x_a and x_b contiguously,
         *           i.e. _x = [x_a, x_b], i.e. _x is of dimension 4
         * \param _coeffs stores the constant k,
         *                i.e. _coeffs[0] = k */
        inline virtual double eval_f(const Vec &_x, const Vec &_coeffs) override {
            //------------------------------------------------------//
            // implement the function f(x) = 1/2 * k * ((x[0] - x[2])^2 + (x[1] - x[3])^2)
            const double &k = _coeffs[0];
            return 0.5 * k * (std::pow(_x[0] - _x[2], 2) + std::pow(_x[1] - _x[3], 2));
            //------------------------------------------------------//
        }

        /** evaluates the spring element's energy gradient
         * \param _x contains x_a and x_b contiguously,
         *           i.e. _x = [x_a, x_b], i.e. _x is of dimension 4
         * \param _coeffs stores the constant k,
         *                i.e. _coeffs[0] = k
         * \param _g the output gradient, which should also be of dimension 4 */
        inline virtual void eval_gradient(const Vec &_x, const Vec &_coeffs, Vec &_g) override {
            //------------------------------------------------------//
            // implement the gradient and store in _g
            const double &k = _coeffs[0];

            const double ax = _x[0];
            const double ay = _x[1];

            const double bx = _x[2];
            const double by = _x[3];

            const double deriv_ax = k * (ax - bx);
            const double deriv_ay = k * (ay - by);

            const double deriv_bx = k * (bx - ax);
            const double deriv_by = k * (by - ay);

            _g << deriv_ax, deriv_ay, deriv_bx, deriv_by;

            //------------------------------------------------------//
        }


        /** evaluates the spring element's energy Hessian
         * \param _x contains x_a and x_b contiguously,
         *           i.e. _x = [x_a, x_b], i.e. _x is of dimension 4
         * \param _coeffs stores the constant k,
         *                i.e. _coeffs[0] = k
         * \param _H the output Hessian, which should be a 4x4 Matrix */
        inline virtual void eval_hessian(const Vec &_x, const Vec &_coeffs, Mat &_H) override {
            //------------------------------------------------------//
            // implement the hessian matrix and store in _H
            const double &k = _coeffs[0];

            const double ax = _x[0];
            const double ay = _x[1];

            const double bx = _x[2];
            const double by = _x[3];

            const double deriv_ax_ax = k;
            const double deriv_ax_ay = 0;
            const double deriv_ax_bx = -k;
            const double deriv_ax_by = 0;

            const double deriv_ay_ax = 0;
            const double deriv_ay_ay = k;
            const double deriv_ay_bx = 0;
            const double deriv_ay_by = -k;

            const double deriv_bx_ax = -k;
            const double deriv_bx_ay = 0;
            const double deriv_bx_bx = k;
            const double deriv_bx_by = 0;

            const double deriv_by_ax = 0;
            const double deriv_by_ay = -k;
            const double deriv_by_bx = 0;
            const double deriv_by_by = k;

            _H << deriv_ax_ax, deriv_ax_ay, deriv_ax_bx, deriv_ax_by,
                    deriv_ay_ax, deriv_ay_ay, deriv_ay_bx, deriv_ay_by,
                    deriv_bx_ax, deriv_bx_ay, deriv_bx_bx, deriv_bx_by,
                    deriv_by_ax, deriv_by_ay, deriv_by_bx, deriv_by_by;
            //------------------------------------------------------//
        }
    };

//=============================================================================
}


