Series S07 aopt: L-BFGS and Gauss Newton Method

Schaller Alain, Fontana Jonas, Carrel Vincent


Affine Invariance
===
- See the attached part1.pdf file.






Programming
===
- See the attached part1.pdf file.






Programming
===
After changing a few arguments in the tests, everything looks good,
same for the plots.
-- Did retake the code from the web app, to produce the same plots
locally and more easily iterate through results.
-- Finally, also added a few targets to automatically produce all
output for algorithms comparison and process them with an AWK program.

# Standard Newton Method
Like said on the forum, feel a bit weird to give a start a
code, where variables are not necessarily used like is in our
implementation, same for the tests to fail. Those small things
make it misleading and make us rethink all the things.



# Newton method with modified Hessian
Again, the prepared identity variable isn't even necessary
when we read properly the documentation of the LLT class.


# Newton Methods vs Gradient descent
Find the different results in the "results" folder, especially
the "newton_vs_gradient_results.pdf" file which includes all plots
and finally the table of all values.

# Bonus
Implemented and passing tests.
- not additional comments
