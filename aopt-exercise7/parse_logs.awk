function print_group_header(group, headers_t) {
    n = split(headers_t, headers, " ")
    for (i = 1; i <= n; i++) {
        printf "\t" group " " headers[i]
    }
}

#BEGIN {
#    printf "method\targ1\targ2\tgrid side\ttotal time[s]"
#    print_group_header("total time evaluation", "time[s] percentage[%%]")
#    print_group_header("eval_f", "time[s] evals avg[s]")
#    print_group_header("eval_grad", "time[s] evals avg[s] factor")
#    print_group_header("eval_hess", "time[s] evals avg[s] factor")
#    printf "\n"
#}

$1 == "ARGUMENTS:" {
    EXECUTABLE=$2
    ARG1=$3
    ARG2=$4
    GRID_SIDE=$5

    # set FN_INDEX based on argument position for specific executable
    if (EXECUTABLE == "NewtonMethods")
        FN_INDEX = ARG2
    else
        FN_INDEX = ARG1

#    if (EXECUTABLE == "NEWTON_METHOD") {
#        EXECUTABLE_ID=1
#    }
#    else if (EXECUTABLE == "GRADIENT_DESCENT") {
#        EXECUTABLE_ID=2
#    }
#    else if (EXECUTABLE == "LBGFS") {
#        EXECUTABLE_ID=3
#    }
#    else { #
#        EXECUTABLE_ID=4
#    }

    next
}

$1 == "ITERATION:" {
    ITERATION=$2

    # If NetwonMethods, has an additional argument:
    #   arg1:
    #       0: standard Newton
    #       1: projected hessian

    # If Gradient descent, select either spring constraint scenario:
    #   arg2:
    #       1: corners
    #       2: sides
    # else
    #   other method will default to scenario 1: corners
    NEWTON_METHOD = ""
    if (EXECUTABLE == "NewtonMethods") {
        if (ARG1 == "0") {
            NEWTON_METHOD = "Standard Newton"
        } else {
            NEWTON_METHOD = "Newton with projected Hessian"
        }
    } else if (EXECUTABLE == "GaussNewton") {
        # Gauss Newton use standard Newton
        NEWTON_METHOD = "Standard Newton"
    }

    # CONSTRAINT_SCENARIO
    CONSTRAINT_SCENARIO = "Corners" # default
    if (EXECUTABLE == "GradientDescent") {
        if (ARG1 == "1")
            CONSTRAINT_SCENARIO = "Corners"
        else
            CONSTRAINT_SCENARIO = "Sides"
    }

    # FN_INDEX Description
    if (FN_INDEX == "0") {
        FN_INDEX_DESC = "0: f without length"
    } else if (FN_INDEX == "1") {
        FN_INDEX_DESC = "1: f with length"
    } else {
        FN_INDEX_DESC = "2: f with length with positive local hessian"
    }

    # For LBFGS, ARG2 holds the history
    HISTORY="-"
    if (EXECUTABLE == "LBFGS")
        HISTORY=ARG2

    next
}

$0 ~ /total time +:/ {
    total_time = $NF
    gsub(/s/, "", total_time)
    next
}

$0 ~ /total time evaluation +:/ {
    eval_time = $5
    gsub(/s/, "", eval_time)
    percentage = $6
    gsub(/\(/, "", percentage)
    next
}

$0 ~ /eval_f time +:/ {
    f_time = $4
    gsub(/s/, "", f_time)
    f_evals = $7
    f_avg = $10
    gsub(/s/, "", f_avg)
    next
}

$0 ~ /eval_grad time:/ {
    g_time = $3
    gsub(/s/, "", g_time)
    g_evals = $6
    g_avg = $9
    gsub(/s,/, "", g_avg)
    g_factor = $11
    gsub(/)/, "", g_factor)
    next
}

$0 ~ /eval_hess time:/ {
    h_time = $3
    gsub(/s/, "", h_time)
    h_evals = $6
    h_avg = $9
    gsub(/s,/, "", h_avg)
    h_factor = $11
    gsub(/)/, "", h_factor)

    # Executable\tArg1\tArg2\tFunction index description\tNewton method\tSpring constraint scenario\tHistory (m)\tGrid side\tIterations\tTotal time [s]\tTotal time evaluation [s]\tTotal time evaluation percentage [%]\teval_f time [s]\teval_f evaluations\teval_f avg [s]\teval_g time [s]\teval_g evaluations\teval_g avg [s]\teval_g factor\teval_h time [s]\teval_h evaluations\teval_h avg [s]\teval_h factor

    printf EXECUTABLE "\t" ARG1 "\t" ARG2 "\t" FN_INDEX_DESC "\t" NEWTON_METHOD "\t" CONSTRAINT_SCENARIO "\t" HISTORY "\t" GRID_SIDE "\t" ITERATION
    printf "\t" total_time
    printf "\t" eval_time "\t" percentage
    printf "\t" f_time "\t" f_evals "\t" f_avg
    printf "\t" g_time "\t" g_evals "\t" g_avg "\t" g_factor
    printf "\t" h_time "\t" h_evals "\t" h_avg "\t" h_factor "\n"
    next
}

# Example of timing stats:
#   total time    : 0.027232s
#   total time evaluation : 0.003943s  (14.4793 %)
#   eval_f time   : 0.00058s  ( #evals: 285 -> avg 0.00000s )
#   eval_grad time: 0.00018s  ( #evals: 96 -> avg 0.00000s, factor: 0.90086)
#   eval_hess time: 0.00319s  ( #evals: 96 -> avg 0.00003s, factor: 16.31277)
