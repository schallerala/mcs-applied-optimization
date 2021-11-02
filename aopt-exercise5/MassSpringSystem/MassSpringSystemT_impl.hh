#define MASSSPRINGSYSTEM_C

#include "MassSpringSystemT.hh"

namespace AOPT {


    template<class MassSpringProblem>
    double MassSpringSystemT<MassSpringProblem>::initial_system_energy() const {
        if (msp_ != nullptr) {
            Vec points = get_spring_graph_points();
            return msp_.get()->eval_f(points);
        }

        return -1;
    }

    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::set_spring_graph_points(const Vec &_points) {
        int n_vertices = sg_.n_vertices();

        for (size_t i = 0; i < n_vertices; ++i)
            sg_.set_vertex(i, Point(_points[2 * i], _points[2 * i + 1]));
    }

    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::save_spring_system(const char *_filename) const {
        sg_.save_to_files(_filename);
    }

    template<class MassSpringProblem>
    std::shared_ptr<MassSpringProblem> MassSpringSystemT<MassSpringProblem>::get_problem() const {
        return msp_;
    }

    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::setup_problem(const int _spring_element_type, const bool _least_square) {
        //set unknown variable number
        n_unknowns_ = 2 * sg_.n_vertices();

        //initialize the problem pointer
        //for least square problem (Gauss-Newton)
        if (_least_square) {

        } else { //for normal problem
            if (_spring_element_type == WITH_LENGTH) {
                msp_ = std::make_shared<MassSpringProblem>(sewl_, n_unknowns_);
            } else if (_spring_element_type == WITHOUT_LENGTH) {
                msp_ = std::make_shared<MassSpringProblem>(se_, n_unknowns_);
            } else {
                std::cout << "Error: spring function index should be 0, 1 or 2!" << std::endl;
                return;
            }
        }



        //add spring elements
        for (size_t i = 0; i < sg_.n_edges(); ++i)
            msp_.get()->add_spring_element(sg_.from_vertex(i), sg_.to_vertex(i), sg_.coefficient(i), sg_.length(i));
    }


    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::add_constrained_spring_elements(const enum ConstraintsScenario _scenario) {
        //------------------------------------------------------//
        // add constrained spring elements to the problem
        //implement both scenarios here.

        const auto WEIGHT = 1E5;

        // Didn't get why I didn't get to do a switch case of the enums value
        if (_scenario == FOUR_CORNERS) {
            // Attach the four corner nodes in a m x n spring graph to target coordinates
            // (0,0), (m, 0), (0, n) and (m, n) respectively.
            const auto tl = get_grid_index(0, 0);
            const auto tr = get_grid_index(0, n_grid_x_);
            const auto bl = get_grid_index(n_grid_y_, 0);
            const auto br = get_grid_index(n_grid_y_, n_grid_x_);

            msp_->add_constrained_spring_element(tl, WEIGHT, 0, 0);
            msp_->add_constrained_spring_element(tr, WEIGHT, n_grid_y_, 0);
            msp_->add_constrained_spring_element(bl, WEIGHT, 0, n_grid_x_);
            msp_->add_constrained_spring_element(br, WEIGHT, n_grid_y_, n_grid_x_);

        } else if (_scenario == SIDES) {
            for (size_t i = 0; i <= n_grid_y_; ++i) {
                {
                    const auto n_i_0_index = get_grid_index(i, 0);
                    const auto &point_i_0 = sg_.point(n_i_0_index);
                    msp_->add_constrained_spring_element(n_i_0_index, WEIGHT, point_i_0.x(), point_i_0.y());
                }
                {
                    msp_->add_constrained_spring_element(get_grid_index(i, n_grid_x_), WEIGHT, i, n_grid_x_);
                }
            }

        } else {
            throw std::invalid_argument("Expecting valid scenario enum as argument!");
        }


        //------------------------------------------------------//


    }


    template<class MassSpringProblem>
    void MassSpringSystemT<MassSpringProblem>::setup_spring_graph() {
        //------------------------------------------------------//
        //add vertices
        for (int j = 0; j <= n_grid_y_; ++j)
            for (int i = 0; i <= n_grid_x_; ++i)
                sg_.add_vertex(Point(i, j));

        //add edges
        for (int j = 0; j < n_grid_y_; ++j) {
            for (int i = 0; i < n_grid_x_; ++i) {
                //horizontal edge
                sg_.add_edge(get_grid_index(i, j), get_grid_index(i + 1, j), 1., 1.);
                //vertical edge
                sg_.add_edge(get_grid_index(i, j), get_grid_index(i, j + 1), 1., 1.);
                //diagonal edge
                sg_.add_edge(get_grid_index(i, j), get_grid_index(i + 1, j + 1), 1., sqrt(2.));
                //diagonal edge
                sg_.add_edge(get_grid_index(i + 1, j), get_grid_index(i, j + 1), 1., sqrt(2.));
            }
        }

        //add right most
        for (int j = 0; j < n_grid_y_; ++j)
            sg_.add_edge(get_grid_index(n_grid_x_, j), get_grid_index(n_grid_x_, j + 1), 1., 1.);

        //add top cap
        for (int i = 0; i < n_grid_x_; ++i)
            sg_.add_edge(get_grid_index(i, n_grid_y_), get_grid_index(i + 1, n_grid_y_), 1., 1.);

        //------------------------------------------------------//
    }

    template<class MassSpringProblem>
    typename MassSpringSystemT<MassSpringProblem>::Vec
    MassSpringSystemT<MassSpringProblem>::get_spring_graph_points() const {
        Vec points(n_unknowns_);
        int n_vertices = sg_.n_vertices();

        for (size_t i = 0; i < n_vertices; ++i) {
            points[2 * i] = sg_.point(i)[0];
            points[2 * i + 1] = sg_.point(i)[1];
        }

        return points;
    }

    template<class MassSpringProblem>
    int MassSpringSystemT<MassSpringProblem>::get_grid_index(const int _i, const int _j) const {
        assert(_i <= n_grid_x_ && _j <= n_grid_y_);
        return (n_grid_x_ + 1) * _j + _i;
    }

    template<class MassSpringProblem>
    size_t MassSpringSystemT<MassSpringProblem>::n_grid_points() const {
        return (n_grid_x_ + 1) * (n_grid_y_ + 1);
    }

    template<class MassSpringProblem>
    size_t MassSpringSystemT<MassSpringProblem>::n_edges() const {
        return sg_.n_edges();
    }
}
