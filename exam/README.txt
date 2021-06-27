The last two digits of my student number is 98, so my project is 10, yet another cubic spline

I have written a small report where I do the theoretical derivation of the expression for the coefficients.
This is contained in the file "projectderivation.pdf"

The function that builds the quadratic interpolating polynomial that is then used to estimate the first derivative of
each point is named "estimate_derivative" and is found in "subSpline.c". Functions for building and evaluating the cubic
spline are also found in "subSpline.c", and are named "initialize_sub_spline" and "evaluate_sub_spline", respectively.

I included the functions for building quadratic and cubic splines from the homework "Interpolation", since I needed
some of the functions I had made earlier, and I wanted to compare to my cubic spline.

I used the spline on tabulated data from the cosine function, and some random data with a jump that I made up.
The results can be seen in subSplineplot.png and subSplinePlotJump.png

The project could have been tested on more data, but I am happy with the result and evaluate myself at 6 or 7 out of 10

