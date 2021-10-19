#pragma once

#include <FunctionBase/FunctionBaseSparse.hh>
#include <FunctionBase/ParametricFunctionBase.hh>
//#include "ConstrainedSpringElement2D.hh"

//== NAMESPACES ===============================================================

namespace AOPT {


//== CLASS DEFINITION =========================================================

/* This problem is very similar to the MassSpringProblem2DDense except it uses
 * a sparse-matrix-based representation of its data.
 * For more information on sparse matrices, please refer to Eigen's documentation:
 * https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
 *
 * The eval_f() and eval_gradient() should be the same as MSP2DDense since sparse-ness
 * only concerns matrices. */
    class MassSpringProblem2DSparse : public FunctionBaseSparse {
    public:
        using Vec = FunctionBaseSparse::Vec;
        using Mat = Eigen::MatrixXd;

        using Edge = std::pair<int, int>;
        // sparse matrix type
        using SMat = FunctionBaseSparse::SMat;
        // triplet
        using T = FunctionBaseSparse::T;

        MassSpringProblem2DSparse(ParametricFunctionBase &_spring, const int _n_unknowns) :
                FunctionBaseSparse(),
                n_(_n_unknowns),
                func_(_spring) {
            xe_.resize(func_.n_unknowns());
            ge_.resize(func_.n_unknowns());
            he_.resize(func_.n_unknowns(), func_.n_unknowns());

        }

        ~MassSpringProblem2DSparse() {}

        virtual int n_unknowns() override {
            return n_;
        }

        /** evaluates the spring element's energy, which is the sum of the energy
         * of all its springs.
         *
         * \param _x the problem's springs positions.
         *           It should contain the positions of all nodes of the system.
         *           i.e. (_x[2*i], _x[2*i+1]) is the position of the i-th node
         * \return the sum of the energy of all the springs */
        virtual double eval_f(const Vec &_x) override {
            double energy(0);

            //used to store the value of k and l, i.e. coeff[0] = ks_[i];
            Vec coeff(2);

            //------------------------------------------------------//
            // assemble function values of all spring elements
            //use vector xe_ to store the local coordinates of two nodes of every spring
            //then pass it to func_.eval_f(...)

            for (size_t i = 0; i < springs_.size(); ++i) {
                const auto& spring = springs_[i];

                const int from_i = spring.first;
                const int to_i = spring.second;

                xe_ << _x[2 * from_i], _x[2 * from_i + 1], _x[2 * to_i], _x[2 * to_i + 1];
                coeff << ks_[i], ls_[i];
                
                energy += func_.eval_f(xe_, coeff);
            }

            //------------------------------------------------------//


            return energy;
        }


        /** The problem's energy gradient is a composition of the individual gradient
         * of each of its springs.
         *
         * \param _x the problem's springs positions.
         *           It should contain the positions of all nodes of the system.
         *           i.e. (_x[2*i], _x[2*i+1]) is the position of the i-th node
         *
         * \param _g should contain the successive per-node gradients.
         *           i.e. (_g[2*i], _g[2*i+1]) is the sum of gradients of all the
         *           springs connected to the i-th node (see handout) */
        virtual void eval_gradient(const Vec &_x, Vec &_g) override {
            _g.resize(n_unknowns());
            _g.setZero();
            Vec coeff(2);

            //------------------------------------------------------//
            // assemble local gradient vector to the global one
            //use ge_ to store the result of the local gradient
            // gonna be the same as the dense!

            for (size_t i = 0; i < springs_.size(); ++i) {
                const auto& spring = springs_[i];

                const int index_ax = 2 * spring.first;
                const int index_ay = index_ax + 1;

                const int index_bx = 2 * spring.second;
                const int index_by = index_bx + 1;

                xe_ << _x[index_ax], _x[index_ay], _x[index_bx], _x[index_by];
                coeff << ks_[i], ls_[i];
                ge_.setZero();

                func_.eval_gradient(xe_, coeff, ge_);

                _g[index_ax] += ge_[0];
                _g[index_ay] += ge_[1];
                _g[index_bx] += ge_[2];
                _g[index_by] += ge_[3];
            }

            //------------------------------------------------------//
        }


        /** The problem's energy Hessian is a composition of the individual Hessian
         * of each of its springs.
         * This should be the same as in MassSpringProblem2DDense except for
         * how the resulting matrix is filled. Please refer to Eigen's official
         * documentation if the usage of sparse matrices is still not clear to you
         *
         * \param _x the problem's springs positions.
         *           It should contain the positions of all nodes of the system.
         *           i.e. (_x[2*i], _x[2*i+1]) is the position of the i-th node **/
        virtual void eval_hessian(const Vec &_x, SMat &_h) override {
            std::vector<T> triplets;
            triplets.reserve(16 * springs_.size());

            //used to store the value of k and l, i.e. coeff[0] = ks_[i], coeff[0] = ls_[i];
            Vec coeff(2);
            //------------------------------------------------------//
            // assemble local hessian matrix to the global one
            //use he_ to store the local hessian matrix

            for (size_t i = 0; i < springs_.size(); ++i) {
                const auto& spring = springs_[i];

                const int index_ax = 2 * spring.first;
                const int index_ay = index_ax + 1;

                const int index_bx = 2 * spring.second;
                const int index_by = index_bx + 1;

                xe_ << _x[index_ax], _x[index_ay], _x[index_bx], _x[index_by];
                coeff << ks_[i], ls_[i];

                func_.eval_hessian(xe_, coeff, he_);

                auto push_block = [&triplets](const int i, const int j, const Eigen::Block<Mat, 2, 2, false>& block_view)
                {
                    triplets.push_back(T(i, j, block_view(0, 0)));
                    triplets.push_back(T(i, j + 1, block_view(0, 1)));
                    triplets.push_back(T(i + 1, j, block_view(1, 0)));
                    triplets.push_back(T(i + 1, j + 1, block_view(1, 1)));
                };

                {
                    const auto ii = he_.block<2, 2>(0, 0);
                    push_block(index_ax, index_ax, ii);
                }
                {
                    const auto ij = he_.block<2, 2>(0, 2);
                    push_block(index_ax, index_bx, ij);
                }
                {
                    const auto ji = he_.block<2, 2>(2, 0);
                    push_block(index_bx, index_ax, ji);
                }
                {
                    const auto jj = he_.block<2, 2>(2, 2);
                    push_block(index_bx, index_bx, jj);
                }
            }


            //------------------------------------------------------//



            _h.resize(n_unknowns(), n_unknowns());
            _h.setZero();
            _h.setFromTriplets(triplets.begin(), triplets.end());
        }


        void add_spring_element(const int _v_idx0, const int _v_idx1, const double _k = 1., const double _l = 1.) {
            if (2 * _v_idx0 > (int) n_ || _v_idx0 < 0 || 2 * _v_idx1 >= (int) n_ || _v_idx1 < 0)
                std::cout << "Warning: invalid spring element was added... " << _v_idx0 << " " << _v_idx1 << std::endl;
            else {
                springs_.emplace_back(_v_idx0, _v_idx1);
                ks_.push_back(_k);
                ls_.push_back(_l);
            }
        }

//        void add_constrained_spring_element(const int _v_idx, const double _w = 1., const double _px = 0.,
//                                            const double _py = 0.) {
//            if (2 * _v_idx > (int) n_ || _v_idx < 0)
//                std::cout << "Warning: invalid node constraint element was added... " << _v_idx << std::endl;
//            else {
//                attached_node_indices_.push_back(_v_idx);
//                weights_.push_back(_w);
//                desired_points_.push_back(_px);
//                desired_points_.push_back(_py);
//            }
//        }


    private:
        int n_;
        std::vector<Edge> springs_;

        ParametricFunctionBase &func_;

        //vector of constants
        std::vector<double> ks_;
        std::vector<double> ls_;

        // coordinates of spring nodes (local)
        Vec xe_;
        // gradient of each spring element
        Vec ge_;
        // hessian of each spring element
        Mat he_;


//        std::vector<int> attached_node_indices_;

//        //vector of constants
//        std::vector<double> weights_;
//        std::vector<double> desired_points_;

        // coordinates of each attached node (local)
//        Vec cs_xe_;
//        // gradient of each attached node
//        Vec cs_ge_;
//        // hessian of each node constraint element
//        Mat cs_he_;
    };

//=============================================================================
}


