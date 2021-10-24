#ifndef AOPT_EXERCISES_EXERCISE2_H
#define AOPT_EXERCISES_EXERCISE2_H

#include <vector>
#include <iostream>
#include <Utils/OptimalityChecker.hh>
#include <Functions/FunctionQuadraticND.hh>
#include <Functions/FunctionQuadratic2D.hh>

bool evaluateExercise2(int _argc = 0,
                       const char *const *_argv = {}) {
    //-------------------------------------------------------------------------------//
    // set up the optimization problem by populating the various
    //matrices and vectors coefficients

    // --------
    // 1. set objective function: x_1^2 - 2 * x_2^2
    const auto obj_func = new AOPT::FunctionQuadratic2D(-2);

    // --------
    // 2. inequality constraints

    // 2.1 inequality 1: (x_1 + 4)^2 - 2 <= x_2
    AOPT::FunctionBase::Mat A_ineq1(2, 2);
    A_ineq1 << 2, 0,
            0, 0;

    AOPT::FunctionBase::Vec b_ineq1(2);
    b_ineq1 << 8, -1;

    double c_ineq1 = 14;

    // 2.2 inequality 2: x_1 >= -10
    AOPT::FunctionBase::Mat A_ineq2(2, 2);
    A_ineq2.setZero();

    AOPT::FunctionBase::Vec b_ineq2(2);
    b_ineq2 << -1, 0;

    double c_ineq2 = -10;

    const std::vector<AOPT::FunctionBase *> ineq_cons = {
            new AOPT::FunctionQuadraticND(A_ineq1, b_ineq1, c_ineq1),
            new AOPT::FunctionQuadraticND(A_ineq2, b_ineq2, c_ineq2)
    };

    // --------
    // 3. equality constraints

    // 3.1 equality: x_1 - x_2 + 4 = 0
    AOPT::FunctionBase::Mat A_eq(2, 2);
    A_eq.setZero();

    AOPT::FunctionBase::Vec b_eq(2);
    b_eq << 1, -1;

    double c_eq = 4;
    const std::vector<AOPT::FunctionBase *> eq_cons = {
            new AOPT::FunctionQuadraticND(A_eq, b_eq, c_eq)
    };

    // --------
    // 4. set lambdas and vs
    AOPT::FunctionBase::Vec lambda(2);
    lambda << 4, 0;

    AOPT::FunctionBase::Vec nu(1);
    nu << -12;

    // 5. set query point
    AOPT::FunctionBase::Vec query_point(2);
    query_point << -2, 2;

    //-------------------------------------------------------------------------------//

    auto oc = _argc < 2
              // not given any argument, use default constructor
              ? AOPT::OptimalityChecker()
              : AOPT::OptimalityChecker(strtod(_argv[1], nullptr));

    if (oc.is_KKT_satisfied(obj_func, ineq_cons, eq_cons, query_point, lambda, nu)) {
        std::cout << "\nThe query point satisfies the KKT condition." << std::endl;
        return true;
    } else {
        std::cout << "\nThe query point does NOT satisfy the KKT condition." << std::endl;
        return false;
    }
}

#endif //AOPT_EXERCISES_EXERCISE2_H
