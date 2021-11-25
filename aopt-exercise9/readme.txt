Series S09 aopt: Equality Constrained Optimization II

Carrel Vincent, Fontana Jonas, Schaller Alain


Problem 1: Constraints Elimination
===
- See document Exercise9_FontanaSchallerCarrel.pdf.




Problem 2: Solving the KKT system
===
- See document Exercise9_FontanaSchallerCarrel.pdf.







Programming
===


# Infeasible Start Newton Method

Some defined variables where maybe confusing, but after some reflection/work
got around it. End up passing all tests and all execution produce realistic
mass spring systems resolution.




# Hybrid Newton Method

After following the first part of the programming exercise, the bonus was more
than self-explanatory. Just too bad with the current structure of the code, we
can't reuse the first part of the exercise to progress towards feasible x; then
just call the previous implementation.
