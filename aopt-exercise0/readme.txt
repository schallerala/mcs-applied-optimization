serie S01 aopt

Schaller Alain, Fontana Jonas, Carrel Vincent

1. Implementing objective functions

We just had to implement the evaluate functions, so we just put the correct return formula
The screenshot can be found in the /output folder

2. Grid Search

Searching over 2D is pretty simple, just iterate over both dimensions to test each "steps" at each grid corner, check for min and keep the lowest

For ND, we can't just nest N for loop together, so we had to use recursion to ensure testing each cuboïd's corner. 
We opted for a recursive approch because it is straight forward, do the job and is easy (relatively) to understand when code reading


##### WARNING ####

We modified the unit_test.cc file on line 104. The matrix A was supposed to be a 5x5 Matrix, but 30 coefficients were given (6x5). We first tried to remove the additionnal column at the right side, 
but we were getting another value when testing, and then realized that when removing the last 5 coefficients instead, we were getting the expected result