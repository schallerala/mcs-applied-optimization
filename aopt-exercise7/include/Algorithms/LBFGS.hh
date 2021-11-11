#pragma once

#include <Eigen/Core>
#include <Algorithms/LineSearch.hh>
#include <FunctionBase/FunctionBaseSparse.hh>

//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================

    class LBFGS {
    public:
        using Vec = FunctionBaseSparse::Vec;
        using Mat = FunctionBaseSparse::Mat;
        using MapVec = Eigen::Map<Vec>;

        LBFGS(const int _m) : m_(_m) {}


        /**
        * @brief solve
        * \param _problem pointer to any function/problem inheriting from FunctionBaseSparse
        *        on which the basic Newton Method will be applied
        * \param _initial_x starting point of the method
        * \param _eps epsilon under which the method stops
        * \param _max_iters maximum iteration of the method*/
        template<class Problem>
        Vec solve(Problem *_problem, const Vec &_initial_x, const double _eps = 1e-4, const int _max_iters = 1000000) {
            std::cout << "******** LBFGS ********" << std::endl;

            int n = _problem->n_unknowns();

            //allocate storage
            init_storage(n);

            //get starting point
            Vec x = _initial_x;

            if (m_ < 1) {
                std::cout << "\nError: m should be larger than 0!" << std::endl;
                return x;
            }

            //squared epsilon for stopping criterion
            double e2 = _eps * _eps;

            //allocate gradient storage
            Vec g(n), sk(n), yk(n);

            //initialize k
            int k(0);

            //initialize
            double f = _problem->eval_f(x);
            _problem->eval_gradient(x, g);

            xp_ = x;
            gp_ = g;


            do {
                double g2 = g.squaredNorm();

                // print status
                std::cout << "iter: " << k <<
                          "   obj = " << f <<
                          "   ||g||^2 = " << g2 << std::endl;

                if (g2 < e2) {
                    std::cout << "Gradient norm converges!" << std::endl;
                    return x;
                }

                //compute r_
                //------------------------------------------------------//
                two_loop_recursion(g, sk, yk, k);
                //------------------------------------------------------//

                //compute the step size
                double t = LineSearch::backtracking_line_search(_problem, x, g, -r_, 1.);
//                double t = LineSearch::wolfe_line_search(_problem, x, g, -r_, 1.);

                if (t < 1e-16) {
                    std::cout << "The step length is too small!" << std::endl;
                    return x;
                }

                //update previous function value, x and gradient
                fp_ = f;
                xp_ = x;
                gp_ = g;

                //update x
                x -= t * r_;

                //evaluate current f
                f = _problem->eval_f(x);

                if (k > 0 && fp_ <= f) {
                    std::cout << "Function value converges!" << std::endl;
                    return x;
                }

                //current gradient
                _problem->eval_gradient(x, g);

                //update storage
                sk = x - xp_;
                yk = g - gp_;

                //------------------------------------------------------//
                update_storage(g, sk, yk, k);
                //------------------------------------------------------//

                k++;

            } while (k < _max_iters);


            return x;
        }


    private:
        void two_loop_recursion(const Vec &_g, const Vec &_sk, const Vec &_yk, const int _k) {
            //------------------------------------------------------//
            // implement the two-loop recursion as described in the lecture slides

            // juste so close to what is done on different project:
            // https://github.com/hjmshi/PyTorch-LBFGS/blob/master/functions/LBFGS.py#L323-L333
            // FIXME however, doesn't work, both test fail, with even once the direction going
            //      towards a greater objective function (error message from backtracking line search)

            // q = gradient f_k
            Vec q = _g;

            // for i=k - 1 ... k - m do
            // --> from `last` to `oldest`
            // --> navigate from right to left in the history
            for (int i = std::min(m_ - 1, _k - 1); i >= 0; --i) {
                // alpha_i = rho_i * s_i^T * q <-- computed in LBFGS::update_storage
                // q = q - alpha_i * y_i
                q -= alpha_[i] * mat_y_.col(i);
            }

            // H_k^0 = gamma_k * I
            // gamma_k = (s_{k-1}^T * y_{k - 1}) / (y_{k - 1}^T * y_{k - 1}) <-- computed in outer loop
            // gamma_k used to choose H_k^0
            const double gamma_k = _k > 0
                                 // wait for after the 1st iteration, as it the history will be empty at first
                                 ? (double)(_sk.transpose() * _yk) / (_yk.transpose() * _yk)
                                 : 1.;

            const Mat H_k_0 = gamma_k * Mat::Identity(_g.size(), _g.size());

            std::cout << "H_k_0\n" << H_k_0 << std::endl;

            // r = H_k^0 * q
            r_ = H_k_0 * q;

            // for i = k - m ... k - 1 do
            // --> from `oldest` to `last`
            // --> navigate from left to right in the history
            if (_k > 0) {
                for (size_t i = 0; i < std::min(_k, m_); ++i) {
                    // rho_k = 1 / (y_k^T * s_k) <-- computed in LBFGS::update_storage
                    // beta = rho_i * y_i^T * r
                    const double beta = rho_[i] * mat_y_.col(i).transpose() * r_;
                    // r = r + s_i * (alpha_i - beta);
                    r_ += mat_s_.col(i) * (alpha_[i] - beta);
                }
            }

            std::cout << "r_:\n" << r_ << std::endl;

            //------------------------------------------------------//
        }

        void update_storage(const Vec &_g, const Vec &_sk, const Vec &_yk, const int _k) {
            assert(m_ > 0);

            //------------------------------------------------------//
            // update the si and yi stored in the mat_s_ and mat_y_ respectively
            const auto shift_left_matrix = [_k](Mat &matrix, const int m) {
                const int last = std::min(m - 1, _k);
                if (_k < m)
                    return last;

                for (size_t i = 0; i < last; ++i) {
                    matrix.col(i) = matrix.col(i + 1);
                }
                return last;
            };

            const auto shift_up_vector = [_k](Vec &vector, const int m) {
                const int last = std::min(m - 1, _k);
                if (_k < m)
                    return last;

                for (size_t i = 0; i < last; ++i) {
                    vector[i] = vector[i + 1];
                }
                return last;
            };

            // 1. shift sk in history
            const auto next_col = shift_left_matrix(mat_s_, m_);
            // 2. shift yk in history
            shift_left_matrix(mat_y_, m_);

            // 3. append sk in history
            mat_s_.col(next_col) = _sk;
            // 4. append yk in history
            mat_y_.col(next_col) = _yk;

            // 5. compute new rho: rho_k = 1 / (y_k^T * s_k)
            const auto next_rho = 1 / (_yk.transpose() * _sk);
            // 7. shift rho in history
            const auto next_row = shift_up_vector(rho_, m_); // next_row will be same as next_col, but still store
                                                                // variable for clarity
            // 8. add new rho in history
            rho_[next_row] = next_rho;

            // 9. compute new alpha: alpha_i = rho_i * s_i^T * q, with q = g
            const auto next_alpha = next_rho * _sk.transpose() * _g;
            // 10. shift alpha in history
            shift_up_vector(alpha_, m_);
            // 10. add new alpha in history
            alpha_[next_row] = next_alpha;

            std::cout << "history sk\n" << mat_s_ << std::endl;
            std::cout << "history yk\n" << mat_y_ << std::endl;
            std::cout << "history rho\n" << rho_ << std::endl;
            std::cout << "history alpha\n" << alpha_ << std::endl;
            //------------------------------------------------------//
        }


        void init_storage(const int _n) {
            mat_y_.resize(_n, m_);
            mat_s_.resize(_n, m_);
            rho_.resize(m_);
            alpha_.resize(m_);

            // for readability
            rho_.setZero();
            mat_y_.setZero();
            mat_s_.setZero();
            alpha_.setZero();
        }

    private:
        //variables' name convention follow the lecture slides
        //number of most recent pairs of s and y to store
        int m_;
        //store y
        Mat mat_y_;
        //store s
        Mat mat_s_;
        //previous function value
        double fp_;
        //previous x
        Vec xp_;
        //previous gradient
        Vec gp_;
        //move direction
        Vec r_;

        Vec alpha_;
        Vec rho_;

    };

//=============================================================================
}





