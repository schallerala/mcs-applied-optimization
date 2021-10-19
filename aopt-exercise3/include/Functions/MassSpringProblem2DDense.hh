#pragma once

#include <FunctionBase/FunctionBase.hh>
#include <FunctionBase/ParametricFunctionBase.hh>

//#include "FixingNodeElement.hh"

//== NAMESPACES ===============================================================

namespace AOPT {


//== CLASS DEFINITION =========================================================

    /* Defines a 2D Mass Spring Problem based on a dense matrix representation of its data.
     * It is represented as an energy function where its value is the sum of the energy
     * of all of its springs.
     * Its gradient and Hessian, however, are a bit more complicated as they are assembled
     * from the gradient and Hessian of each spring individually.
     * See the two functions below for more details.
     *
     * The MSP uses a Parametric Function to represent its inner springs, which is given
     * to the constructor. This function can, for instance, be a SpringElement2D or a
     * SpringElement2DWithLength, which themselves extend the ParametricFunctionBase interface.
     *
     * Additionally, the set of springs is represented by an array of Edges, which are simply
     * pairs or indices of the vertices connected by each edge. */
    class MassSpringProblem2DDense : public FunctionBase {
    public:
        using Vec = FunctionBase::Vec;
        using Mat = FunctionBase::Mat;
        using Edge = std::pair<int, int>;

        /** Base constructor defining the Parametric Function that will be used to
         * evaluate the springs' energy*/
        MassSpringProblem2DDense(ParametricFunctionBase &_spring, const int _n_unknowns) :
                FunctionBase(),
                n_(_n_unknowns),
                func_(_spring) {
            xe_.resize(func_.n_unknowns());
            ge_.resize(func_.n_unknowns());
            he_.resize(func_.n_unknowns(), func_.n_unknowns());
        }

        ~MassSpringProblem2DDense() {}

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

            //used to store the value of k and l, i.e. coeff[0] = ks_[i], coeff[0] = ls_[i];
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

            //used to store the value of k and l, i.e. coeff[0] = ks_[i];
            Vec coeff(2);

            //------------------------------------------------------//
            // assemble local gradient vector to the global one
            //use ge_ to store the result of the local gradient

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
         * \param _x the problem's springs positions.
         *           It should contain the positions of all nodes of the system.
         *           i.e. (_x[2*i], _x[2*i+1]) is the position of the i-th node **/
        virtual void eval_hessian(const Vec &_x, Mat &_h) override {
            _h.resize(n_unknowns(), n_unknowns());
            _h.setZero();

            //used to store the value of k and l, i.e. coeff[0] = ks_[i];
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

                {
                    const auto ii = he_.block<2, 2>(0, 0);
                    _h(index_ax, index_ax) += ii(0, 0);
                    _h(index_ax, index_ax + 1) += ii(0, 1);
                    _h(index_ax + 1, index_ax) += ii(1, 0);
                    _h(index_ax + 1, index_ax + 1) += ii(1, 1);
                }
                {
                    const auto ij = he_.block<2, 2>(0, 2);
                    _h(index_ax, index_bx) += ij(0, 0);
                    _h(index_ax, index_bx + 1) += ij(0, 1);
                    _h(index_ax + 1, index_bx) += ij(1, 0);
                    _h(index_ax + 1, index_bx + 1) += ij(1, 1);
                }
                {
                    const auto ji = he_.block<2, 2>(2, 0);
                    _h(index_bx, index_ax) += ji(0, 0);
                    _h(index_bx, index_ax + 1) += ji(0, 1);
                    _h(index_bx + 1, index_ax) += ji(1, 0);
                    _h(index_bx + 1, index_ax + 1) += ji(1, 1);
                }
                {
                    const auto jj = he_.block<2, 2>(2, 2);
                    _h(index_bx, index_bx) += jj(0, 0);
                    _h(index_bx, index_bx + 1) += jj(0, 1);
                    _h(index_bx + 1, index_bx) += jj(1, 0);
                    _h(index_bx + 1, index_bx + 1) += jj(1, 1);
                }
            }

            //------------------------------------------------------//
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


    private:
        int n_;
        std::vector<Edge> springs_;

        ParametricFunctionBase &func_;

        //vector of constants
        std::vector<double> ks_;
        std::vector<double> ls_;

        // coordinates of each node
        Vec xe_;
        // gradient of each spring element
        Vec ge_;
        // hessian of each spring element
        Mat he_;
    };

//=============================================================================
}


