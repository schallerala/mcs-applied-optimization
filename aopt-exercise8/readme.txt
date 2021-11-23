Series S08 aopt: Trust Region and Equality Constrained Optimization

Carrel Vincent, Fontana Jonas, Schaller Alain


Trust Region Subproblem: Exact Solution
===
- See document Exercise8_FontanaSchallerCarrel.pdf.







Programming
===

# Equality Constrained Newton Method
Too bad there is no dedicated test for solving without needing a projection.

For the projection, you can find details of the reflection for the solution in the
code as comments.
Additionally, emphasizing here again that the use of the `COLAMDOrdering< Index >`
as template parameter and pass `int` as underlying template parameter is a bit
confusing, but on the surface, makes sense, plus it passes the test and executing
the mass spring system with that solver produce meaningful results.