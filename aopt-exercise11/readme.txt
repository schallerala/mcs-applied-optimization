Series S11 aopt: Inequality Constrained Optimization

Carrel Vincent, Fontana Jonas, Schaller Alain


Algorithms Overview
===
- Not done






Programming
===

# Interior Point Method

I did all my comments already on the forum. Really regret the confusion of
µ and the orientation.

You can find the output mass spring problem plot in the "results" folder.



## Executable logs

.../Build/bin/InteriorPoint 2 1 10 10 10000 interior-point-2-1-10-10
Saving initial spring graph to interior-point-2-1-10-10_*.csv
******** Interior Point ********
******** Newton Method with projected hessian ********
******** Newton Method with projected hessian ********
******** Newton Method with projected hessian ********
******** Newton Method with projected hessian ********
******** Newton Method with projected hessian ********
******** Newton Method with projected hessian ********
******** Newton Method with projected hessian ********
******** Newton Method with projected hessian ********
######## Timing statistics ########
total time    : 13.6988s
total time evaluation : 13.4756s  (98.3709 %)
eval_f time   : 0.15606s  ( #evals: 277 -> avg 0.00056s )
eval_grad time: 0.18925s  ( #evals: 53 -> avg 0.00357s, factor: 6.33785)
eval_hess time: 13.13033s  ( #evals: 53 -> avg 0.24774s, factor: 439.72345)
iter: 2   obj = 11886.8   ||lambda||^2 = 1.80789e-09   n_projection_steps = 0
Saving optimized spring graph to interior-point-2-1-10-10_opt_*.csv

Process finished with exit code 0
