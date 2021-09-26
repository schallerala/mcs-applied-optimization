serie S01 aopt

Schaller Alain, Fontana Jonas, Carrel Vincent

1. Implementing objective functions

We just had to implement the evaluate functions, so we just put the correct return formula

2. Grid Search

Searching over 2D is pretty simple, just iterate over both dimensions to test each "steps" at each grid corner, check for min and keep the lowest

For ND, we can't just nest N for loop together, so we had to use recursion to ensure testing each cuboïd's corner. Struggled a bit with the recursion due to our C++ level being a bit rusty. 