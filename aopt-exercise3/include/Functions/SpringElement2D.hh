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

            const double deriv_ax = 0.5 * k * (2 * ax - 2 * bx);
            const double deriv_ay = 0.5 * k * (2 * ay - 2 * by);

            const double deriv_bx = 0.5 * k * (2 * bx - 2 * ax);
            const double deriv_by = 0.5 * k * (2 * by - 2 * ay);

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
            //Todo: implement the hessian matrix and store in _H


            //------------------------------------------------------//
        }
    };

//=============================================================================
}


