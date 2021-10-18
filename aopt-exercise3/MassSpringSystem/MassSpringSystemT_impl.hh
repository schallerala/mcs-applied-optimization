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
    void MassSpringSystemT<MassSpringProblem>::setup_spring_graph() {
        //------------------------------------------------------//
        // set up the spring graph of n_grid_x by n_grid_y ()

        const double diag_l = std::sqrt(2);

        const int y_limit = n_grid_y_ + 1;
        const int x_limit = n_grid_x_ + 1;

        for (size_t y = 0; y < y_limit; ++y) {
            for (size_t x = 0; x < x_limit; ++x) {
                auto *pt = new Point(x, y);
                int last_edge_index = sg_.add_vertex(*pt);

                // connect horizontally to previous
                if (x > 0)
                    sg_.add_edge(last_edge_index - 1, last_edge_index);

                if (y > 0) {
                    // connect vertically to edge on top
                    sg_.add_edge(last_edge_index - x_limit, last_edge_index);

                    if (/* not first row && */ last_edge_index % x_limit > 0) {
                        // add top-left 2 bottom-right arrow (anchor tl)
                        sg_.add_edge(last_edge_index - x_limit - 1, last_edge_index, 1, diag_l);
                    }
                    if (/* not first row && */ (last_edge_index + 1) % x_limit > 0) {
                        // add bottom-left 2 top-right arrow (anchor bl)
                        sg_.add_edge(last_edge_index, last_edge_index - x_limit + 1, 1, diag_l);
                    }
                }
            }
        }

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
