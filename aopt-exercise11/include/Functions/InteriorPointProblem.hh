#pragma once


#include <FunctionBase/FunctionBaseSparse.hh>


//== NAMESPACES ===============================================================

namespace AOPT {

    //== CLASS DEFINITION =========================================================
    class InteriorPointProblem : public FunctionBaseSparse {
    public:
        // default constructor
        InteriorPointProblem(FunctionBaseSparse *_obj, const std::vector<FunctionBaseSparse *> &_constraints)
                : FunctionBaseSparse(), obj_(_obj), constraints_(_constraints), t_(1.0) {
            v_ = Vec(obj_->n_unknowns());
            M_ = SMat(obj_->n_unknowns(), obj_->n_unknowns());
            N_ = SMat(obj_->n_unknowns(), obj_->n_unknowns());
        }

        // defualt destructor
        virtual ~InteriorPointProblem() {};

        // access log-barrier parameter
        double &t() { return t_; }

        // number of unknowns
        virtual int n_unknowns() override { return obj_->n_unknowns(); }

        // function evaluation
        virtual double eval_f(const Vec &_x) override {
            //------------------------------------------------------//
            // add function value (objective function + barrier function)

            // minimize f(x) - 1 / t \sum_{i=1}^m log(-g_i(x))

            const auto eval_f = obj_->eval_f(_x);

            double constraints_eval_f = 0;

            for (const auto constraint: constraints_) {
                constraints_eval_f += log_of_minus_function(constraint, _x);
            }

            return eval_f - 1 / t_ * constraints_eval_f;

            //------------------------------------------------------//
        }

        // gradient evaluation
        virtual void eval_gradient(const Vec &_x, Vec &_g) override {
            //------------------------------------------------------//
            // add gradients (objective function + barrier function)

            // d/dx central path = d/dx f(x) - 1/t * \sum_{i=1}^m d/dx (log(-g_i(x)))

            _g.setZero();
            // 1. gradient of objective function
            obj_->eval_gradient(_x, _g);

            // to overcome the necessity to multiple all gradient summed iteratively by -1 / t,
            // divide in advance the current result and multiply the final result.
            _g /= -1 / t_;

            // 2. get all constraints' gradient and sum it
            for (const auto constraint: constraints_) {
                add_grad_of_log_of_function(constraint, _x, _g);
            }

            // 3. multiply the sum of constraints
            _g *= -1 / t_;

            //------------------------------------------------------//
        }

        // hessian matrix evaluation
        virtual void eval_hessian(const Vec &_x, SMat &_H) override {
            //------------------------------------------------------//
            // add hessian matrices (objective function + barrier function)

            // almost same as for gradient

            // d/d²x central path = d/d²x f(x) - 1/t * \sum_{i=1}^m d/d²x (log(-g_i(x)))
            // =>           d/d²x = -1/t * ( (d/d²x f(x))/(-1/t) + \sum_{i=1}^m d/d²x (log(-g_i(x))) )

            _H.setZero();
            // 1. hessian of objective function
            obj_->eval_hessian(_x, _H);

            // to overcome the necessity to multiple all hessian summed iteratively by -1 / t,
            // divide in advance the current result and multiply the final result.
            _H /= -1 / t_;

            // 2. get all constraints' hessian and sum it
            for (const auto constraint: constraints_) {
                add_hess_of_log_of_function(constraint, _x, _H);
            }

            // 3. multiply the sum of constraints
            _H *= -1 / t_;

            //------------------------------------------------------//
        }

        void setT(double t) {
            t_ = t;
        }


    private:
        static double log_of_minus_function(FunctionBaseSparse *_o, const Vec &_x) {
            auto f = _o->eval_f(_x);
            return log(-f);
        }

        void add_grad_of_log_of_function(FunctionBaseSparse *_o, const Vec &_x, Vec &_g) {
            v_.setZero();
            _o->eval_gradient(_x, v_);

            _g += 1.0 / _o->eval_f(_x) * v_;
        }

        void add_hess_of_log_of_function(FunctionBaseSparse *_o, const Vec &_x, SMat &_H) {
            // get f, grad, hess
            double d = _o->eval_f(_x);

            _o->eval_gradient(_x, v_);
            _o->eval_hessian(_x, M_);

            triplets_.clear();
            triplets_.reserve(36);

            //N = v*v.transpose()
            N_.setZero();
            for (auto i = 0u; i < v_.size(); ++i) {
                if (v_[i] != 0) {
                    for (auto j = 0u; j < v_.size(); ++j) {
                        if (v_[j] != 0) {
                            triplets_.emplace_back(i, j, v_[i] * v_[j]);
                        }
                    }
                }
            }
            N_.setFromTriplets(triplets_.begin(), triplets_.end());
            _H += (1.0 / d) * M_ - (1.0 / (d * d)) * N_;
        }

    private:

        // objective function
        FunctionBaseSparse *obj_;

        // constraint functions
        std::vector<FunctionBaseSparse *> constraints_;

        // log barrier parameter
        double t_;


        // temp variables
        Vec v_;
        SMat M_;
        SMat N_;
        std::vector<T> triplets_;
    };
    //=============================================================================

}







