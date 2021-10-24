Series S04 aopt: Duality

Schaller Alain, Fontana Jonas, Carrel Vincent


Lagrange Duality
===
- See the attached part1.pdf file.






KKT Condition
===
- See the attached part1.pdf file.








Programming
===
The exercise 2 checking the KKT conditions have been extracted
from the main method into its own method to be able to include
its execution in the gtest suite.

Concerning the `OptimalityChecker::is_KKT_satisfied` method,
to our opinion, it makes little sense to call the method with
no inequality functions and/or no equality functions
but with a none empty lambda and/or nu vector.
As the iteration of an empty vector will be "equal" to check
its size and "jump" to the last portion of the method.


# Remark:
Would it be possible to remove all previous "todos" comments
for the future exercises?
To avoid them to be index by the IDE and for us to jump to all of
them and check if there is anything to do there.
A bit in the same sens, could you take out the eigen tutorial part?
So we don't have their test waiting with a user input for example.
Thank you.