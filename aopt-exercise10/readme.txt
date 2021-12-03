Series S10 aopt: Inequality Constrained Optimization

Carrel Vincent, Fontana Jonas, Schaller Alain


Active set method for convex quadratic problem
===
- See document Exercise9_FontanaSchallerCarrel.pdf.




"Phase I" method for initial feasible point
===
- See document Exercise9_FontanaSchallerCarrel.pdf.





TODO! Done bonus?!

"Big M" method for initial feasible point
===
- See document Exercise9_FontanaSchallerCarrel.pdf.







Programming
===


# Augmented Lagrangian Method

Regret a bit the presentation of the exercise in the instruction as it doesn't
feel the natural order of implementing things if we follow strictly the bullet
points.

Also, no certain to totally get if we should iterate at the start of the
augmented lagrangian optimization loop, as I didn't get over what should
we iterate; therefore, didn't. The following wording in the lecture notes
confused me:

     find approximate solution [...] and stopping when [...]

Also, regret that "details" like the possibility to start with the nu vector
with zero, is valid. Did research for more than an hour to only re-watch the
theory and hear the teacher's brief comment about it.

You can find the output mass spring problem plot in the "results" folder.


## Executable logs

.../Build/bin/AugmentedLagrangian 2 12 12 10000 augmented-lagrangian-2-12-12
Saving initial spring graph to augmented-lagrangian-2-12-12_*.csv
******** Augmented Lagrangian ********
******** Newton Method with projected hessian ********
******** Newton Method with projected hessian ********
******** Newton Method with projected hessian ********
######## Timing statistics ########
total time    : 3.71651s
total time evaluation : 3.35294s  (90.2174 %)
eval_f time   : 0.04436s  ( #evals: 87 -> avg 0.00051s )
eval_grad time: 0.06660s  ( #evals: 34 -> avg 0.00196s, factor: 3.84198)
eval_hess time: 3.24198s  ( #evals: 31 -> avg 0.10458s, factor: 205.10486)
Saving optimized spring graph to augmented-lagrangian-2-12-12_opt_*.csv

Process finished with exit code 0
