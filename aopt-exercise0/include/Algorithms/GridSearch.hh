#pragma once

#include <Functions/FunctionQuadratic2D.hh>
#include <Functions/FunctionQuadraticND.hh>
#include <Functions/FunctionNonConvex2D.hh>
#include <vector>

//== NAMESPACES ===================================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================
    class GridSearch {
    private:
        int n_grid_;

    public:
        using Vec = FunctionBase::Vec;
        using Mat = FunctionBase::Mat;

        explicit GridSearch(const int _n_grid = 10) : n_grid_(_n_grid) {}

        ~GridSearch() = default;

    public:

        /** Evaluation of a 2D function over the whole grid to find its minimum
         *
         * \param _func a pointer to any 2D function inheriting from FunctionBase
         *              (see include/FunctionBase/FunctionBase.hh)
         * \param _x_l the coordinates of the lower corner of the grid
         * \param _x_u the coordinates of the upper corner of the grid.
         *             _x_l and _x_u together define a square in which the grid lies
         * \return 0 if all went well, -1 if not.*/
        int grid_search_2d(FunctionBase *_func, const Vec &_x_l, const Vec &_x_u, double &_f_min) const {
            std::cout << "Grid searching the minimum of a 2-D function..." << std::endl;
            double f = 0., fmin = std::numeric_limits<double>::max();

            Vec x_min(2);

            Vec xy(2);

            //------------------------------------------------------//
            //Todo: implement the 2d version of the grid search
            // algorithm to find minimum value of _func between _x_l and _x_u
            //------------------------------------------------------//
            Vec incr = (_x_u - _x_l) / n_grid_;
            Vec grid_xy(2);
            for (int x_i = 0; x_i <= n_grid_; x_i++) {
                for (int y_i = 0; y_i <= n_grid_; y_i++) {
                    grid_xy << x_i, y_i;

                    // coefficient wise multiplication
                    xy = _x_l + (incr.cwiseProduct(grid_xy));

                    f = _func->eval_f(xy);
                    if (f < fmin) {
                        fmin = f;
                        x_min = xy;
                    }
                }
            }

            //------------------------------------------------------//
            _f_min = fmin;
            std::cout << "Minimum value of the function is: " << fmin << " at x:\n" << x_min << std::endl;
            return 0;
        }


        /** Evaluation of an ND function over the whole grid to find its minimum
         *  using an iterative approach
         *
         * \param _func a pointer to any ND function inheriting from FunctionBase
         *              (see include/FunctionBase/FunctionBase.hh)
         * \param _x_l the coordinates of the lower corner of the grid
         * \param _x_u the coordinates of the upper corner of the grid.
         *             _x_l and _x_u together define an ND cuboid in which the grid lies
         * \return 0 if all went well, -1 if not.*/
        int grid_search_nd(FunctionBase *_func, const Vec &_x_l, const Vec &_x_u, double &_f_min) const {
            int n = _func->n_unknowns();
            if (_x_l.size() != n || _x_u.size() != n) {
                std::cout << "Error: input limits are not of correct dimension!" << std::endl;
                return -1;
            }
            std::cout << "Grid searching the minimum of a " << n << "-D function..." << std::endl;

            double f_min = std::numeric_limits<double>::max();
            Vec x_min(n);
            //------------------------------------------------------//
            //Todo: implement the nd version of the grid search
            // algorithm to find minimum value of a nd quadratic function
            // set f_min with the minimum, which is then stored in the referenced argument _f_min


            //------------------------------------------------------//
            // only 1 point to compute, which will be the only solution
            if (_x_l == _x_u) {
                f_min = _func->eval_f(_x_l);
                x_min = _x_l;
            }
            else {
                long iteration_limit = 1;
                for (size_t i = 0; i < n; i++) {
                    iteration_limit *= n_grid_;
                }

                // represent the size of an n dimension increment to move on the grid on all dimensions
                Vec incr = (_x_u - _x_l) / n_grid_;

                // represent the points on the grid
                Vec grid_x(n);
                grid_x.fill(0);

                // the corresponding x to feed into the function
                Vec x(n);
                x << _x_l;

                // value of 1 evaluation
                double f;

                const int grid_n_index = n_grid_ - 1;

                for (long i = 0; i < iteration_limit; i++) {

                    // coefficient wise multiplication
                    x = _x_l + (incr.cwiseProduct(grid_x));

                    f = _func->eval_f(x);
                    if (f < f_min) {
                        f_min = f;
                        x_min = x;
                    }

                    // -- progress on the grid position
                    // increment first and following if "reset" to 0
                    // for example with grid_n_index = 4
                    // 000 -> 001
                    // 004 -> 010
                    for (Eigen::Index j = 0; j < n; ++j) {
                        if (grid_x[j] == grid_n_index) {
                            // restart and increment next
                            grid_x[j] = 0;
                        } else {
                            grid_x[j]++;
                            break;
                        }
                    }
                }
            }

            _f_min = f_min;
            std::cout << "Minimum value of the function is: " << f_min << " at x:\n" << x_min << std::endl;

            return 0;
        }

    };

    //=============================================================================
}





