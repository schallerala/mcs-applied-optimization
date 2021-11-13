# https://stackoverflow.com/a/27158086/3771148
function ltrim(s) { sub(/^[ \t\r\n]+/, "", s); return s }
function rtrim(s) { sub(/[ \t\r\n]+$/, "", s); return s }
function trim(s) { return rtrim(ltrim(s)); }

function print_group_header(group, headers_t) {
    n = split(headers_t, headers, " ")
    for (i = 1; i <= n; i++) {
        printf "\t" group " " headers[i]
    }
}

BEGIN {
    printf "method\targ1\targ2\tgrid side\ttotal time[s]"
    print_group_header("total time evaluation", "time[s] percentage[%%]")
    print_group_header("eval_f", "time[s] evals avg[s]")
    print_group_header("eval_grad", "time[s] evals avg[s] factor")
    print_group_header("eval_hess", "time[s] evals avg[s] factor")
    printf "\n"
}

$1 == "Arguments:" {
    ARG1=$2
    ARG2=$3
    GRID_SIDE=$4
}

$0 ~ /\*\*\*\*\*\*\*\*/ {
    gsub(/\*/, "")
    method = trim($0)
    printf method "\t" ARG1 "\t" ARG2 "\t" GRID_SIDE
    next
}

$0 ~ /total time +:/ {
    time = $NF
    gsub(/s/, "", time)
    printf "\t" time
    next
}

$0 ~ /total time evaluation +:/ {
    time = $5
    gsub(/s/, "", time)
    percentage = $6
    gsub(/\(/, "", percentage)
    printf "\t" time "\t" percentage
    next
}

$0 ~ /eval_f time +:/ {
    f_time = $4
    gsub(/s/, "", f_time)
    f_evals = $7
    f_avg = $10
    gsub(/s/, "", f_avg)
    printf "\t" f_time "\t" f_evals "\t" f_avg
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
    printf "\t" g_time "\t" g_evals "\t" g_avg "\t" g_factor
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
    printf "\t" h_time "\t" h_evals "\t" h_avg "\t" h_factor "\n"
    next
}

# total time    : 0.027232s
# total time evaluation : 0.003943s  (14.4793 %)
# eval_f time   : 0.00058s  ( #evals: 285 -> avg 0.00000s )
# eval_grad time: 0.00018s  ( #evals: 96 -> avg 0.00000s, factor: 0.90086)
# eval_hess time: 0.00319s  ( #evals: 96 -> avg 0.00003s, factor: 16.31277)
