#pragma once

#include "SpringElement2DWithLength.hh"

#include <Eigen/Eigenvalues>
//== NAMESPACES ===============================================================

namespace AOPT {

//== CLASS DEFINITION =========================================================

/**
*   Class that overrides the hessian of the non-convex energy of the spring element
 * by fixing the negative eigen values of the hessian matrix
*/

    class SpringElement2DWithLengthPSDHess : public SpringElement2DWithLength {
    public:

        SpringElement2DWithLengthPSDHess() : SpringElement2DWithLength() {}

        inline virtual void eval_hessian(const Vec &_x, const Vec &_coeffs, Mat &_H) override {
            SpringElement2DWithLength::eval_hessian(_x, _coeffs, _H);

            // slide 18 and http://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html
            //------------------------------------------------------//
            // compute the hessian matrix and project it to a positve definite matrix
            //Hint: 1. to compute the eigen decomposition, use
            //          Eigen::SelfAdjointEigenSolver<Mat> solver(A);
            Eigen::SelfAdjointEigenSolver<Mat> solver(_H);

            //          Mat evecs = solver.eigenvectors();  //this matrix contains the eigenvectors in its columns
            //          Vec evals = solver.eigenvalues();
            //      2. to convert a vector d to a (dense) diagonal matrix D, use
            //          D = d.asDiagonal()
            const auto& D = solver.eigenvalues();
            const auto& Q = solver.eigenvectors();

            Vec new_m_i(D.size());

            // M = Q^T * diag(m_i) * Q
            //        +- 0,                     lambda_i >= epsilon
            // m_i = <
            //        +- epsilon - lambda_i,    lambda_i < epsilon
            for (size_t i = 0; i < D.size(); ++i) {
                const auto lambda_i = D[i];
                const double m_i = lambda_i >= m_eps
                             ? 0
                             : m_eps - lambda_i;
                new_m_i[i] = m_i;
            }

            _H = Q.transpose() * new_m_i.asDiagonal() * Q;

            //------------------------------------------------------//
        }

        static constexpr double m_eps = 1e-7;
    };

//=============================================================================

}
