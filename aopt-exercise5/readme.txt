Series S05 aopt: Line Search Methods

Schaller Alain, Fontana Jonas, Carrel Vincent


Exact line search for the convex quadratic function
===
- See the attached part1.pdf file.






Gradient descent with exact line search
===
- See the attached part1.pdf file.






Programming Exercise: Modified Mass Spring System
===
Was soon done with the task but the failing test and believing we could
trust the default parameter values for `alpha` and `tau` made us lose a
lot of time always making us reconsider what was implemented.

Also, regret that related data, like the points under constraints are presented
under different fields, while they should always be initialized and consumed
together. In the same sense, it is unclear wether or not we should use the
implemented `ConstrainedSpringElement2D` in the `MassSpringProblem2DSparse`
as no field is declared before hand. Ending up with some duplicated code
(this has been emphasised in the code too).

Finally, was hard to follow the instruction sheet, like the constraints scenario
number 2, ones targeting an other point and in the second phase a fix coordinate.
Hope the results are really what was expected.


# Starting points affecting the gradient descent method
Functions with multiple local minimum, depending on the starting point, we might end
in one minimum or the other which might not be the optimal value. Also, the
backtracking line search will output related distances, due to the evaluation of
the function with this starting point.


- See the part2.pdf for the results