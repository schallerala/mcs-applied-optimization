#pragma once

#include <iostream>
#include <iomanip>
#include <limits>
#include <cfloat>
#include <cmath>
#include "RandomNumberGenerator.hh"


//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================

    class DerivativeChecker {
    public:
        using Vec = Eigen::VectorXd;
        using SMat = Eigen::SparseMatrix<double>;

        DerivativeChecker() : rng_(-1., 1.) {}

        ~DerivativeChecker() {}

        struct Config {
            Config() : x_min(-1.0), x_max(1.0), n_iters(1), dx(1e-5), eps(1e-3),
                       relativeEps(std::numeric_limits<double>::quiet_NaN()) {}

            double x_min;
            double x_max;
            int n_iters;
            double dx;
            double eps;
            double relativeEps;
        };

        template<class ProblemInterface>
        bool check_all(ProblemInterface &_np, double _dx, double _eps) {
            conf_.dx = _dx;
            conf_.eps = _eps;
            return check_all(_np);
        }

        template<class ProblemInterface>
        bool check_all(ProblemInterface &_np) {
            bool d1_ok = check_d1(_np);
            bool d2_ok = check_d2(_np);

            return (d1_ok && d2_ok);
        }

        template<class ProblemInterface>
        bool check_d1(ProblemInterface &_np) {
            int n_ok = 0;
            int n_errors = 0;


            int n = _np.n_unknowns();
            Vec x(n), g(n);

            for (int i = 0; i < conf_.n_iters; ++i) {
                // get random x
                get_random_x(x, conf_.x_min, conf_.x_max);
//                std::cout << "\nx: " << x.transpose() << std::endl;
                // gradient
                _np.eval_gradient(x, g);
//                std::cout << "\ng: " << g.transpose() << std::endl;

                for (int j = 0; j < n; ++j) {
                    x[j] += conf_.dx;
                    double f1 = _np.eval_f(x);
                    x[j] -= 2.0 * conf_.dx;
                    double f0 = _np.eval_f(x);
                    x[j] += conf_.dx;
                    double fd = (f1 - f0) / (2.0 * conf_.dx);

                    if ((!std::isnan(conf_.relativeEps) &&
                         std::abs(fd - g[j]) > std::max(std::abs(g[j]), 1.0) * conf_.relativeEps) ||
                        std::abs(fd - g[j]) > conf_.eps) {
                        ++n_errors;
                        std::cerr << "Gradient error in component " << j << ": " << g[j]
                                  << " should be " << fd << " (" << std::abs(fd - g[j]) << ")" << f0 << " " << f1 << " "
                                  << std::endl;
                    } else ++n_ok;
                }
            }

            std::cerr << "############## Gradient Checker Summary #############\n";
            std::cerr << "#ok   : " << n_ok << std::endl;
            std::cerr << "#error: " << n_errors << std::endl;

            return (n_errors == 0);
        }

        template<class MatrixType>
        inline double getCoeff(const MatrixType &m, int r, int c) {
            return m.coeff(r, c);
        }

        template<class ProblemInterface>
        bool check_d2(ProblemInterface &_np) {
            int n_ok = 0;
            int n_errors = 0;

            const int n = _np.n_unknowns();
            Vec x(n);
            SMat H(n, n);

            for (int i = 0; i < conf_.n_iters; ++i) {
                // get random x
                get_random_x(x, conf_.x_min, conf_.x_max);
                // gradient
                _np.eval_hessian(x, H);

                for (int j = 0; j < n; ++j)
                    for (int k = 0; k < n; ++k) {
                        x[j] += conf_.dx;
                        x[k] += conf_.dx;
                        double f0 = _np.eval_f(x);
                        x[j] -= 2.0 * conf_.dx;
                        double f1 = _np.eval_f(x);
                        x[j] += 2.0 * conf_.dx;
                        x[k] -= 2.0 * conf_.dx;
                        double f2 = _np.eval_f(x);
                        x[j] -= 2.0 * conf_.dx;
                        double f3 = _np.eval_f(x);

                        double fd = (f0 - f1 - f2 + f3) / (4.0 * conf_.dx * conf_.dx);


                        if ((!std::isnan(conf_.relativeEps) && std::abs(fd - H.coeff(j, k)) >
                                                               std::max(std::abs(getCoeff(H, j, k)), 1.0) *
                                                               conf_.relativeEps) ||
                            std::abs(fd - getCoeff(H, j, k)) > conf_.eps) {
                            ++n_errors;
                            std::cerr << "Hessian error in component " << j << "," << k << ": " << getCoeff(H, j, k)
                                      << " should be (following FiniteDifferences) " << fd
                                      << " (absolute delta: " << std::abs(fd - getCoeff(H, j, k))
                                      << ", relative delta:"
                                      << (std::abs(fd - getCoeff(H, j, k)) / std::max(std::abs(getCoeff(H, j, k)), 1.0))
                                      << ")" << std::endl;
                        } else ++n_ok;
                    }
            }

            std::cerr << "############## Hessian Checker Summary #############\n";
            std::cerr << "#ok   : " << n_ok << std::endl;
            std::cerr << "#error: " << n_errors << std::endl;

            return (n_errors == 0);
        }

        Config &config() { return conf_; }

    protected:
        void get_random_x(Vec &_x, double _xmin, double _xmax) {
            // get random values in [min, max]
            _x = rng_.get_random_nd_vector(_x.size());

            double range = _xmax - _xmin;
            for (unsigned int i = 0; i < _x.size(); ++i)
                _x[i] = (((_x[i] + 1.0) / 2.0) * range + _xmin);
        }

    private:
        Config conf_;
        RandomNumberGenerator rng_;
    };
//=============================================================================
}


