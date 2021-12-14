#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>

//== NAMESPACES ===================================================================

#define x0 _x[2 * idx0_]
#define y0 _x[2 * idx0_ + 1]

#define x1 _x[2 * idx1_]
#define y1 _x[2 * idx1_ + 1]

#define x2 _x[2 * idx2_]
#define y2 _x[2 * idx2_ + 1]

namespace AOPT {

    //== CLASS DEFINITION =========================================================
    class AreaConstraint2D : public FunctionBaseSparse {
    public:
        // Area constraint: 1/2*det(v_01 | V02) >= eps
        // f = 1/2*((x1 - x0)(y2 - y0) - (x2 - x0)(y1 - y0)) + eps <= 0
        // constructor
        AreaConstraint2D(const int _n, const int _idx0, const int _idx1, const int _idx2,
                         const double _eps = 1e-10)
                : FunctionBaseSparse(), n_(_n), idx0_(_idx0), idx1_(_idx1), idx2_(_idx2), eps_(_eps) {}

        // number of unknowns
        inline virtual int n_unknowns() override { return n_; }

        // function evaluation
        // _x stores the coordinates of all nodes
        inline virtual double eval_f(const Vec &_x) override {
            //------------------------------------------------------//
            // implement the constraint function value

            return 0.5 * ((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)) + eps_;
            //------------------------------------------------------//
        }

        // gradient evaluation
        inline virtual void eval_gradient(const Vec &_x, Vec &_g) override {
            _g.setZero();
            //------------------------------------------------------//
            // implement the gradient and store in _g

            // d/dx0 = 0.5*y1 - 0.5*y2
            _g[2 * idx0_] = 0.5 * (y1 - y2);

            // d/dy0 = -0.5*x1 + 0.5*x2
            _g[2 * idx0_ + 1] = 0.5 * (x2 - x1);

            // d/dx1 = -0.5*y0 + 0.5*y2
            _g[2 * idx1_] = 0.5 * (y2 - y0);

            // d/dy1 = 0.5*x0 - 0.5*x2
            _g[2 * idx1_ + 1] = 0.5 * (x0 - x2);

            // d/dx2 = 0.5*y0 - 0.5*y1
            _g[2 * idx2_] = 0.5 * (y0 - y1);

            // d/dy2 = -0.5*x0 + 0.5*x1
            _g[2 * idx2_ + 1] = 0.5 * (x1 - x0);

            //------------------------------------------------------//
        }

        // hessian matrix evaluation
        inline virtual void eval_hessian(const Vec &_x, SMat &_h) override {
            _h.setZero();
            //------------------------------------------------------//
            //implement the hessian matrix and store in _h
            // d/dx0 dx0 = 0
//            _h.insert(2 * idx0_, 2 * idx0_) = 0;
            // d/dx0 dy0 = 0
//            _h.insert(2 * idx0_, 2 * idx0_ + 1) = 0;
            // d/dx0 dx1 = 0
//            _h.insert(2 * idx0_, 2 * idx1_) = 0;
            // d/dx0 dy1 = 0.5
            _h.insert(2 * idx0_, 2 * idx1_ + 1) = 0.5;
            // d/dx0 dx2 = 0
//            _h.insert(2 * idx0_, 2 * idx2_) = 0;
            // d/dx0 dy2 = -0.5
            _h.insert(2 * idx0_, 2 * idx2_ + 1) = -0.5;

            // d/dy0 dx0 = 0
//            _h.insert(2 * idx0_ + 1, 2 * idx0_) = 0;
            // d/dy0 dy0 = 0
//            _h.insert(2 * idx0_ + 1, 2 * idx0_ + 1) = 0;
            // d/dy0 dx1 = -0.5
            _h.insert(2 * idx0_ + 1, 2 * idx1_) = -0.5;
            // d/dy0 dy1 = 0
//            _h.insert(2 * idx0_ + 1, 2 * idx1_ + 1) = 0;
            // d/dy0 dx2 = 0.5
            _h.insert(2 * idx0_ + 1, 2 * idx2_) = 0.5;
            // d/dy0 dy2 = 0
//            _h.insert(2 * idx0_ + 1, 2 * idx2_ + 1) = 0;


            // d/dx1 dx0 = 0
//            _h.insert(2 * idx1_, 2 * idx0_) = 0;
            // d/dx1 dy0 = -0.5
            _h.insert(2 * idx1_, 2 * idx0_ + 1) = -0.5;
            // d/dx1 dx1 = 0
//            _h.insert(2 * idx1_, 2 * idx1_) = 0;
            // d/dx1 dy1 = 0
//            _h.insert(2 * idx1_, 2 * idx1_ + 1) = 0;
            // d/dx1 dx2 = 0
//            _h.insert(2 * idx1_, 2 * idx2_) = 0;
            // d/dx1 dy2 = 0.5
            _h.insert(2 * idx1_, 2 * idx2_ + 1) = 0.5;

            // d/dy1 dx0 = 0.5
            _h.insert(2 * idx1_ + 1, 2 * idx0_) = 0.5;
            // d/dy1 dy0 = 0
//            _h.insert(2 * idx1_ + 1, 2 * idx0_ + 1) = 0;
            // d/dy1 dx1 = 0
//            _h.insert(2 * idx1_ + 1, 2 * idx1_) = 0;
            // d/dy1 dy1 = 0
//            _h.insert(2 * idx1_ + 1, 2 * idx1_ + 1) = 0;
            // d/dy1 dx2 = -0.5
            _h.insert(2 * idx1_ + 1, 2 * idx2_) = -0.5;
            // d/dy1 dy2 = 0
//            _h.insert(2 * idx1_ + 1, 2 * idx2_ + 1) = 0;


            // d/dx2 dx0 = 0
//            _h.insert(2 * idx2_, 2 * idx0_) = 0;
            // d/dx2 dy0 = 0.5
            _h.insert(2 * idx2_, 2 * idx0_ + 1) = 0.5;
            // d/dx2 dx1 = 0
//            _h.insert(2 * idx2_, 2 * idx1_) = 0;
            // d/dx2 dy1 = -0.5
            _h.insert(2 * idx2_, 2 * idx1_ + 1) = -0.5;
            // d/dx2 dx2 = 0
//            _h.insert(2 * idx2_, 2 * idx2_) = 0;
            // d/dx2 dy2 = 0
//            _h.insert(2 * idx2_, 2 * idx2_ + 1) = 0;

            // d/dy2 dx0 = -0.5
            _h.insert(2 * idx2_ + 1, 2 * idx0_) = -0.5;
            // d/dy2 dy0 = 0
//            _h.insert(2 * idx2_ + 1, 2 * idx0_ + 1) = 0;
            // d/dy2 dx1 = 0.5
            _h.insert(2 * idx2_ + 1, 2 * idx1_) = 0.5;
            // d/dy2 dy1 = 0
//            _h.insert(2 * idx2_ + 1, 2 * idx1_ + 1) = 0;
            // d/dy2 dx2 = 0
//            _h.insert(2 * idx2_ + 1, 2 * idx2_) = 0;
            // d/dy2 dy2 = 0
//            _h.insert(2 * idx2_ + 1, 2 * idx2_ + 1) = 0;

            //------------------------------------------------------//
        }

    private:
        int n_;
        // index of the nodes
        int idx0_;
        int idx1_;
        int idx2_;

        double eps_;
    };

//=============================================================================
}






