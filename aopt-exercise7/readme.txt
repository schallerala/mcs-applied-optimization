Series S07 aopt: L-BFGS and Gauss Newton Method

Schaller Alain, Fontana Jonas, Carrel Vincent


Secant Equation
===
- See the first pages of the attached Exercise7_CarrelFontanaSchaller.pdf file.





Curvature Condition
===
- See the first pages of the attached Exercise7_CarrelFontanaSchaller.pdf file.






Programming
===


# L-BFGS
Having this class level field to hold values is somewhat not
conventional. Could expect to use vectors and "shift_left"
from the C++20 standard.

Also, having this alpha class level field made us do a mistake
by computing it in the update storage method while it should
only be accessible in the two loop method as it should always
be recomputed in the first loop.

The hint in the instruction sheet about:

    "whenever the back-tracking line search fails to find a step
    length that meets the curvature condition, you need to skip
    the update"

Was unclear if we should change anything in the "main" solve method.
Doing nothing felt like the right thing to do as it passes the tests.
Certainly, there was already a condition for very small "t",
but it is stopping the algorithm and not "skipping the update",
like described in the hint.

Finally, having those class fields is not thread safe and is
unconventional especially compared to the way to the structure
of the wolfe zoom method with a considerable number of parameters.




# Gauss-Newton Method
This was really not easy but a lot a unknown. As announced, it was
tough to get it to work.

It was somewhat unclear to follow the way we should consider eval_f
to be eval_r for this problem. Could make it somewhat clearer to
have a intermediate abstract class which would declare an abstract
method eval_r which would be implemented by the child classes and this
method would be called by eval_f.

Would have been great to have also tests for the intermediate eval_r,
eval_gradient or eval_jacobian of the mass spring problem.

Also, the comment about rj^2 in mass spring problem eval_r for the
constraints was unclear if the value retrieved from the spring
element eval_f value should be scared or not.

Finally, it is really regrettable that their as issues with the tests,
like said last week already, it make us rethink all we did.



# Comparison of different methods
The Gauss Newton executable is based on the standard Netwon Method without
any Hessian modification for any LLT factorization. For this reason
its execution "without length" won't be able to optimize.

- See the attached Exercise7_CarrelFontanaSchaller.pdf file starting from page 3.


# A line search algorithm for the Wolfe conditions
- No comment