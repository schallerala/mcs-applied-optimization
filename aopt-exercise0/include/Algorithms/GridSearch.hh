#pragma once

#include <Functions/FunctionQuadratic2D.hh>
#include <Functions/FunctionQuadraticND.hh>
#include <Functions/FunctionNonConvex2D.hh>
#include <vector>

//== NAMESPACES ===================================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================
    class GridSearch {
    public:
        using Vec = FunctionBase::Vec;
        using Mat = FunctionBase::Mat;

        GridSearch(const int _n_cells = 10) : n_cells_(_n_cells){}
        ~GridSearch() {}

    public:

        /** Evaluation of a 2D function over the whole grid to find its minimum
         *
         * \param _func a pointer to any 2D function inheriting from FunctionBase
         *              (see include/FunctionBase/FunctionBase.hh)
         * \param _x_l the coordinates of the lower corner of the grid
         * \param _x_u the coordinates of the upper corner of the grid.
         *             _x_l and _x_u together define a square in which the grid lies
         * \return 0 if all went well, -1 if not.*/
        int grid_search_2d(FunctionBase* _func, const Vec& _x_l, const Vec& _x_u, double& _f_min) const {
            std::cout<<"Grid searching the minimum of a 2-D function..."<<std::endl;
            // std::cout<<"Vector size:"<<_x_l.rows()<<std::endl;
            double f = 0., fmin = std::numeric_limits<double>::max();
            
            Vec x_min(2);
            
            //------------------------------------------------------//
            //Todo: implement the 2d version of the grid search
            
            // Get grid steps
            double step_x0 = (_x_u(0) - _x_l(0)) / n_cells_;
            double step_x1 = (_x_u(1) - _x_l(1)) / n_cells_;
            
            // Loop over the two dimensions
            for (double i = _x_l(0); i <= _x_u(0); i += step_x0) {
            	for (double j = _x_l(1); j <= _x_u(1); j += step_x1) {
            	
            		// Construct the vector
            		FunctionQuadratic2D::Vec pt(2);
            		pt << i, j;
            		
            		// Check function
            		double new_min = (*_func).eval_f(pt);
            		
            		// Keep the minimum
            		if (new_min < fmin) {
            			fmin = new_min;
            			x_min = pt;
            		}
            		//fmin = std::min(fmin, new_min);
            	}
            }
            // algorithm to find minimum value of _func between _x_l and _x_u
            //------------------------------------------------------//
           
            
            //------------------------------------------------------//
            _f_min = fmin;
            std::cout<<"Minimum value of the function is: "<<fmin<<" at x:\n"<<x_min<<std::endl;
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
        int grid_search_nd(FunctionBase* _func, const Vec& _x_l, const Vec& _x_u, double& _f_min) const {
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

			Vec v;
			f_min = recursive_grid_search_nd(_func, _x_l, _x_u, f_min, &x_min, 0, v);
            
            //------------------------------------------------------//
            _f_min = f_min;
            std::cout << "Minimum value of the function is: " << f_min << " at x:\n" << x_min << std::endl;
            std::cout << "Number of tests: "<<n_tests_<<std::endl;

            return 0;
        }



    private:
        int n_cells_;
        mutable int n_tests_ = 0;
        
        // Recursive function
		double recursive_grid_search_nd(FunctionBase* _func, const Vec& _x_l, const Vec& _x_u, double current_min, Vec* x_min, int index, Vec vals) const {
			
			// Eval function
			if (index == _x_l.rows()) {
				n_tests_++;
				return (*_func).eval_f(vals);
			}
			
			// Should not be necessary, but prevent index errors
			if (index > _x_l.rows()) { return current_min; }
			
			
			double new_min = current_min;
			double step = (_x_u(index) - _x_l(index)) / n_cells_;
			Vec v;
			
			for (double i = _x_l(index); i <= _x_u(index); i += step) {
					
					
				// Add a new point to the vector
				v = vals;
				v.resize(v.rows() + 1);
				v(v.rows() - 1) = i;
				
        		// Recursive call to inspect lower dimensions
        		double m = recursive_grid_search_nd(_func, _x_l, _x_u, new_min, x_min, index + 1, v);
        		
        		// Check if a better result has been found
        		if (m < new_min) {
        		
        			// Save the new min
        			new_min = m;
        			
        			// Save the coordinates
        			// Bad, but... for some reason, *x_min = v only copies the first value
        			for (int j = 0; j < v.rows(); j++) {
        				(*x_min)(j) = v(j);
        			}
        		}
			}
		return new_min;
		}
    };

    //=============================================================================
}





