#pragma once

#include <iostream>
#include <FunctionBase/FunctionBase.hh>

//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================
    /* This implements a generic non-convex function.
     * Its sole purpose is to test the ConvexitTest class */
    class Exercise1Function1D final : public FunctionBase {
    public:
        // f(x,y) = x^2

        // constructor
        Exercise1Function1D() {}

        // number of unknowns
        inline virtual int n_unknowns() { return 1; }

        /** funcion evaluation
         * \param _x the value at which to evaluate the function.
         *           It should be a 2D vector*/
        inline virtual double eval_f(const Vec &_x) {
            const double x = _x[0];

            //------------------------------------------------------//
            return std::pow(x, 2);
            //------------------------------------------------------//
        }

        // gradient evaluation. Not necessary for this function
        inline virtual void eval_gradient(const Vec &_x, Vec &_g) {}

        // hessian matrix evaluation. Not necessary for this function
        inline virtual void eval_hessian(const Vec &_x, Mat &_H) {}
    };

//=============================================================================
}


