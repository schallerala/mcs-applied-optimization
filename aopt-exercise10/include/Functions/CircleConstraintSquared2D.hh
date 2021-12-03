#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>

//== NAMESPACES ===================================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================
    class CircleConstraintSquared2D : public FunctionBaseSparse {
    public:
        using Vec = FunctionBaseSparse::Vec;
        using SMat = FunctionBaseSparse::SMat;

        // f(x,y) = ((x[2*idx]- center_x)^2 + (x[2*idx+1] - center_y)^2 - radius^2)^2
        // constructor
        CircleConstraintSquared2D(const int _n, const int _idx, const double _center_x, const double _center_y,
                                  const double _radius)
                : FunctionBaseSparse(), n_(_n), idx_(_idx), center_x_(_center_x), center_y_(_center_y),
                  radius_(_radius) {}

        // number of unknowns
        inline virtual int n_unknowns() override { return n_; }

        // funcion evaluation
        // _x stores the coordinates of all nodes
        inline virtual double eval_f(const Vec &_x) override {
            //------------------------------------------------------//
            // implement the constraint function value

            // f(x,y) = ((x[2*idx]- center_x)^2 + (x[2*idx+1] - center_y)^2 - radius^2)^2

            return std::pow(std::pow(_x[2 * idx_] - center_x_, 2) + std::pow(_x[2 * idx_ + 1] - center_y_, 2) - std::pow(radius_, 2), 2);

            //------------------------------------------------------//
        }

        // gradient evaluation
        // _g stores the gradient of all nodes
        inline virtual void eval_gradient(const Vec &_x, Vec &_g) override {
            _g.setZero();
            //------------------------------------------------------//
            // implement the gradient and store in _g

            const auto x = _x[2 * idx_];
            const auto y = _x[2 * idx_ + 1];

            // d / dx = (-4*cx + 4*x)*(-r**2 + (-cx + x)**2 + (-cy + y)**2)
            _g[2 * idx_]     = (-4*center_x_ + 4*x)*(-std::pow(radius_, 2) + std::pow(-center_x_ + x, 2) + std::pow(-center_y_ + y, 2));

            // d / dy = (-4*cy + 4*y)*(-r**2 + (-cx + x)**2 + (-cy + y)**2)
            _g[2 * idx_ + 1] = (-4*center_y_ + 4*y)*(-std::pow(radius_, 2) + std::pow(-center_x_ + x, 2) + std::pow(-center_y_ + y, 2));

            //------------------------------------------------------//
        }

        // hessian matrix evaluation
        // _h stores the hessian of all nodes
        inline virtual void eval_hessian(const Vec &_x, SMat &_h) override {
            _h.setZero();
            //------------------------------------------------------//
            // implement the hessian matrix and store in _h

            const auto x = _x[2 * idx_];
            const auto y = _x[2 * idx_ + 1];

            // d / dx dx = -4*r**2 + (-4*cx + 4*x)*(-2*cx + 2*x) + 4*(-cx + x)**2 + 4*(-cy + y)**2
            _h.insert(2 * idx_    , 2 * idx_)     = -4*std::pow(radius_, 2) + (-4*center_x_ + 4*x)*(-2*center_x_ + 2*x) + 4*std::pow(-center_x_ + x, 2) + 4*std::pow(-center_y_ + y, 2);
            // d / dx dy = (-4*cx + 4*x)*(-2*cy + 2*y)
            _h.insert(2 * idx_    , 2 * idx_ + 1) = (-4*center_x_ + 4*x)*(-2*center_y_ + 2*y);

            // d / dy dx = (-2*cx + 2*x)*(-4*cy + 4*y)
            _h.insert(2 * idx_ + 1, 2 * idx_)     = (-2*center_x_ + 2*x)*(-4*center_y_ + 4*y);
            // d / dy dy = -4*r**2 + 4*(-cx + x)**2 + (-4*cy + 4*y)*(-2*cy + 2*y) + 4*(-cy + y)**2
            _h.insert(2 * idx_ + 1, 2 * idx_ + 1) = -4*std::pow(radius_, 2) + 4*std::pow(-center_x_ + x, 2) + (-4*center_y_ + 4*y)*(-2*center_y_ + 2*y) + 4*std::pow(-center_y_ + y, 2);

            //------------------------------------------------------//
        }

    private:
        int n_;
        int idx_;
        double center_x_;
        double center_y_;
        double radius_;
    };

//=============================================================================
}





