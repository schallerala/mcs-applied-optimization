#pragma once

#include <Utils/RandomNumberGenerator.hh>
#include <FunctionBase/FunctionBase.hh>

//== NAMESPACES ===================================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================

    class ConvexityTest {
    public:
        using Vec = FunctionBase::Vec; ///< Eigen::VectorXd
        using Mat = FunctionBase::Mat; ///< Eigen::MatrixXd

        ConvexityTest() {}

        ~ConvexityTest() {}

    public:

        /** Checks whether the function given as argument is convex or not.
         * If it is not, it should output a point not satisfying the convexity property
         * before returning false.
         * \param _function a function pointer that should be any class inheriting
         * from FunctionBase, e.g. FunctionQuadraticND
         * \param min the minimum value of all tested points' coordinate
         * \param max the maximum value of all tested points' coordinate
         * \param n_evals the number of evaluations/samples tested on the
         *        line between the two points on the function */
        static bool isConvex(FunctionBase* _function, const double min = -1000., const double max = 1000., const int n_evals = 10) {
            const int n = _function->n_unknowns();
            //randomly generate number from min to max
            RandomNumberGenerator rng(min, max);
            
            const int max_sampling_points(1000000);

            //------------------------------------------------------//
            //Todo: Add your code here
                        
            double lambda;
            Vec x,y;
            double step = 1.0 / n_evals;

            for(int j = 0; j < max_sampling_points; j++){
                
                //Initialize next sampling
                x = rng.get_random_nd_vector(n);
                y = rng.get_random_nd_vector(n);

                lambda = 0.0;
                
                for(int i = 0; i < n_evals; i++){

                    if(_function->eval_f(lambda * x + (1 - lambda) * y) > lambda * _function->eval_f(x) + (1 - lambda) * _function->eval_f(y)){

                        std::cout << "The function is not convex notably at point: " << lambda * x + (1 - lambda)* y << std::endl;
                        goto stop; 
                    }

                    lambda += step;
                    
                }
            }

            std::cout << "The function is supposed convex after trying 1000000 samples of 2 points with 10 values in between them" << std::endl;

            return true;
            
            stop:
            //------------------------------------------------------//
            return false;
        }


    private:
        static void printPathInfo(FunctionBase::Vec p1, FunctionBase::Vec p2, FunctionBase::Vec p, double t) {
            std::cout << "path: p(t) = (1 - t) * p1 + t * p2; \nwith:\n"
                      << "  p1 = (" << p1.transpose() << ")\n"
                      << "  p2 = (" << p2.transpose() << ")\n"
                      << "  p (t = " << t << ") = (" << p.transpose() << ")" << std::endl;
        }

    };




}
