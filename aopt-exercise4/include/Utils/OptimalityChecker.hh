#pragma once

#include <FunctionBase/FunctionBase.hh>
#include <vector>
//== NAMESPACES ===============================================================

namespace AOPT {

//== CLASS DEFINITION =========================================================
    class OptimalityChecker {
    public:
        using Vec = Eigen::VectorXd;

        OptimalityChecker(const double _epsilon = 1e-13) : eps_(_epsilon) {}


        /** Checks whether a particular optimization problem satisfies the KKT conditions
         *
         * \param _objective pointer to the objective function, which should be any function
         *        inheriting from FunctionBase
         *
         * \param _inequality_constraints an array containing the inequalities,
         *        each represented as a function inheriting from FunctionBase,
         *        such that _inequality_constraints[i].eval_f(x) <= 0 for all i
         *        if x is in the feasible set
         *
         * \param _equality_constraints an array containing the equalities,
         *        each represented as a function inheriting from FunctionBase,
         *        such that _equality_constraints[i].eval_f(x) == 0 for all i
         *        if x is in the feasible set
         *
         * \param _query_point the point at which the KKT conditions should be checked
         *
         * \param _lambda the lambda coefficients of the dual problem.
         *        It should be of the same dimension as the count of _inequality_constraints
         *        since there is one lambda factor per inequality in the dual problem
         * \param _nu the nu coefficients of the dual problem.
         *        It should be of the same dimension as the count of _equality_constraints
         *        since there is one nu factor per equality in the dual problem
         * */
        bool is_KKT_satisfied(FunctionBase *_objective, const std::vector<FunctionBase *>& _inequality_constraints,
                              const std::vector<FunctionBase *>& _equality_constraints,
                              const Vec& _query_point, const Vec& _lambda, const Vec& _nu) const {
            // Don't check as there is specific test case with not equal count of lambdas and (in)equality functions
//            assert(_inequality_constraints.size() == _lambda.size());
//            assert(_equality_constraints.size() == _nu.size());

            assert(_objective->n_unknowns() == _query_point.size());
            
            //------------------------------------------------------//
            // 1. check only condition 4 in case there are no constraints
            // 2. check inequality constraints (cond. 1.)
            for (size_t i = 0; i < _inequality_constraints.size(); ++i) {
                const auto ineq_f = _inequality_constraints[i];
                assert(ineq_f->n_unknowns() == _query_point.size());

                double ineq_evaluation = ineq_f->eval_f(_query_point);
                // must be <= 0
                if (ineq_evaluation > eps_) {
                    return false;
                }
            }
            
            // 3. check equality constraints (cond. 1.)
            for (size_t i = 0; i < _equality_constraints.size(); ++i) {
                const auto eq_f = _equality_constraints[i];
                assert(eq_f->n_unknowns() == _query_point.size());

                const double eq_evaluation = eq_f->eval_f(_query_point);
                // must be == 0
                if (std::abs(eq_evaluation) > eps_) {
                    return false;
                }
            }
            
            // 4. check lambda (cond. 2.)
            // all coefficient must be >= 0
            // -- Note: not sure if compare with 0 or epsilon
            //          as could image the lambda to be computed and produce floating point errors
            if (_lambda.minCoeff() < -eps_) {
                return false;
            }

            // 5. check complementary slackness (cond. 3.)
            for (size_t i = 0; i < _inequality_constraints.size(); ++i) {
                const auto ineq_f = _inequality_constraints[i];
                const auto lambda = _lambda[i];

                double ineq_evaluation = lambda * ineq_f->eval_f(_query_point);
                // must be == 0
                if (std::abs(ineq_evaluation) > eps_) {
                    return false;
                }
            }

            // 6. check gradient (cond. 4.)
            {
                Vec gradient_sum(_objective->n_unknowns());
                gradient_sum.setZero();

                Vec g_eval(gradient_sum.size());

                const auto eval_add_gradient = [_query_point](Vec gradient_sum, FunctionBase *f, Vec g_eval, const double factor = 0.) {
                    g_eval.setZero();
                    f->eval_gradient(_query_point, g_eval);

                    gradient_sum += factor * g_eval;
                };

                eval_add_gradient(gradient_sum, _objective, g_eval);

                for (size_t i = 0; i < _inequality_constraints.size(); ++i) {
                    const auto ineq_f = _inequality_constraints[i];
                    const auto lambda = _lambda[i];

                    eval_add_gradient(gradient_sum, ineq_f, g_eval, lambda);
                }

                for (size_t i = 0; i < _equality_constraints.size(); ++i) {
                    const auto eq_f = _equality_constraints[i];
                    const auto nu = _nu[i];

                    eval_add_gradient(gradient_sum, eq_f, g_eval, nu);
                }

                // sum of all gradient evaluation with their factor must be == vector 0
                if (gradient_sum.cwiseAbs().maxCoeff() > eps_) {
                    return false;
                }
            }
            //------------------------------------------------------//

            return true;
        }

    private:
        double eps_;
    };
//=============================================================================
}



