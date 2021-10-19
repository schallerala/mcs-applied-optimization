#pragma once

#include <FunctionBase/ParametricFunctionBase.hh>

//== NAMESPACES ===============================================================

namespace AOPT {


//== CLASS DEFINITION =========================================================

/* This class evaluates the energy with a spring with length.
 * It is very similar to the SpringElement2D except it has an addition parameter
 * l_ab which represents the length of the spring at rest. */
    class SpringElement2DWithLength : public ParametricFunctionBase {
    public:
        // E'_ab(x) = 1/2 * k * (((x[0] - x[2])^2 + (x[1] - x[3])^2) - l^2)^2
        // constructor
        SpringElement2DWithLength() : ParametricFunctionBase() {}

        // number of unknowns
        inline virtual int n_unknowns() override { return 4; }

        /** evaluates the spring element's energy
         * \param _x contains x_a and x_b contiguously,
         *           i.e. _x = [x_a, x_b], i.e. _x is of dimension 4
         * \param _coeffs stores the constants k and l,
         *                i.e. _coeffs[0] = k, _coeffs[1] = l */
        inline virtual double eval_f(const Vec &_x, const Vec &_coeffs) override {
            //------------------------------------------------------//
            // implement the function f(x) = 1/2 * k * (((x[0] - x[2])^2 + (x[1] - x[3])^2) - l^2)^2

            const double &k = _coeffs[0];
            const double &l = _coeffs[1];

            return 0.5 * k * std::pow(std::pow(_x[0] - _x[2], 2) + std::pow(_x[1] - _x[3], 2) - std::pow(l, 2), 2);
            //------------------------------------------------------//
        }


        /** evaluates the spring element's energy gradient
         * \param _x contains x_a and x_b contiguously,
         *           i.e. _x = [x_a, x_b], i.e. _x is of dimension 4
         * \param _coeffs stores the constants k and l,
         *                i.e. _coeffs[0] = k, _coeffs[1] = l
         * \param _g the output gradient, which should also be of dimension 4 */
        inline virtual void eval_gradient(const Vec &_x, const Vec &_coeffs, Vec &_g) override {
            //------------------------------------------------------//
            // implement the gradient and store in _g

            const double &k = _coeffs[0];
            const double &l = _coeffs[1];

            const double ax = _x[0];
            const double ay = _x[1];

            const double bx = _x[2];
            const double by = _x[3];

            const double rhs = std::pow(ax - bx, 2) + std::pow(ay - by, 2) - std::pow(l, 2);

            const double deriv_ax = 2 * k * (ax - bx) * rhs;
            const double deriv_ay = 2 * k * (ay - by) * rhs;

            const double deriv_bx = - deriv_ax;
            const double deriv_by = - deriv_ay;

            _g << deriv_ax, deriv_ay, deriv_bx, deriv_by;

            //------------------------------------------------------//
        }

        /** evaluates the spring element's energy Hessian
         * \param _x contains x_a and x_b contiguously,
         *           i.e. _x = [x_a, x_b], i.e. _x is of dimension 4
         * \param _coeffs stores the constants k and l,
         *                i.e. _coeffs[0] = k, _coeffs[1] = l
         * \param _H the output Hessian, which should be a 4x4 Matrix */
        inline virtual void eval_hessian(const Vec &_x, const Vec &_coeffs, Mat &_H) override {
            //------------------------------------------------------//
            // implement the hessian matrix and store in _H
            const double &k = _coeffs[0];
            const double &l = _coeffs[1];

            const double ax = _x[0];
            const double ay = _x[1];

            const double bx = _x[2];
            const double by = _x[3];

            // 2 * k * ((ax - bx) ^ 2 + (ay - by) ^ 2 - l ^ 2) + 4 * (ax - bx) ^ 2 * k
            const double deriv_ax_ax = 2 * k * (std::pow(ax - bx, 2) - std::pow(l, 2) + std::pow(ay - by, 2)) + 4 * std::pow(ax - bx, 2) * k;
            // 4 * k * (ax - bx) * (ay - by)
            const double deriv_ax_ay = 4 * k * (ax - bx) * (ay - by);
            // -(2 * k * ((ax - bx) ^ 2 + (ay - by) ^ 2 - l ^ 2)) - 4 * (ax - bx) ^ 2 * k
            const double deriv_ax_bx = -(2 * k * (std::pow(ax - bx, 2) - std::pow(l, 2) + std::pow(ay - by, 2))) - 4 * std::pow(ax - bx, 2) * k;
            // -(4 * k * (ax - bx) * (ay - by))
            const double deriv_ax_by = - deriv_ax_ay;

            // 4 * k * (ax - bx) * (ay - by)
            const double deriv_ay_ax = deriv_ax_ay;
            // 2 * k * ((ax - bx) ^ 2 + (ay - by) ^ 2 - l ^ 2) + 4 * (ay - by) ^ 2 * k
            const double deriv_ay_ay = 2 * k * (std::pow(ax - bx, 2) - std::pow(l, 2) + std::pow(ay - by, 2)) + 4 * std::pow(ay - by, 2) * k;
            // -(4 * k * (ax - bx) * (ay - by))
            const double deriv_ay_bx = - deriv_ax_ay;
            // -(2 * k * ((ax - bx) ^ 2 + (ay - by) ^ 2 - l ^ 2)) - 4 * (ay - by) ^ 2 * k
            const double deriv_ay_by = -(2 * k * (std::pow(ax - bx, 2) - std::pow(l, 2) + std::pow(ay - by, 2))) - 4 * std::pow(ay - by, 2) * k;

            // -(2 * k * ((ax - bx) ^ 2 + (ay - by) ^ 2 - l ^ 2) + 4 * (ax - bx) ^ 2 * k)
            const double deriv_bx_ax = -deriv_ax_ax;
            // -(4 * k * (ax - bx) * (ay - by))
            const double deriv_bx_ay = - deriv_ax_ay;
            // 4 * (ax - bx) ^ 2 * k + 2 * k * ((ax - bx) ^ 2 + (ay - by) ^ 2 - l ^ 2)
            const double deriv_bx_bx = 4 * std::pow(ax - bx, 2) * k + 2 * k * (std::pow(ax - bx, 2) - std::pow(l, 2) + std::pow(ay - by, 2));
            // 4 * k * (ax - bx) * (ay - by)
            const double deriv_bx_by = deriv_ax_ay;

            // -(4 * k * (ax - bx) * (ay - by))
            const double deriv_by_ax = - deriv_ax_ay;
            // -(2 * k * ((ax - bx) ^ 2 + (ay - by) ^ 2 - l ^ 2) + 4 * (ay - by) ^ 2 * k)
            const double deriv_by_ay = - deriv_ay_ay;
            // 4 * k * (ax - bx) * (ay - by)
            const double deriv_by_bx = deriv_ax_ay;
            // 4 * (ay - by) ^ 2 * k + 2 * k * ((ax - bx) ^ 2 - l ^ 2 + (ay - by) ^ 2)
            const double deriv_by_by = 4 * std::pow(ay - by, 2) * k + 2 * k * (std::pow(ax - bx, 2) - std::pow(l, 2) + std::pow(ay - by, 2));

            _H << deriv_ax_ax, deriv_ax_ay, deriv_ax_bx, deriv_ax_by,
                    deriv_ay_ax, deriv_ay_ay, deriv_ay_bx, deriv_ay_by,
                    deriv_bx_ax, deriv_bx_ay, deriv_bx_bx, deriv_bx_by,
                    deriv_by_ax, deriv_by_ay, deriv_by_bx, deriv_by_by;
            //------------------------------------------------------//
        }
    };

//=============================================================================
}


