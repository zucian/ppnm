#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>

double myErf(double);
double myGamma(double);

int main(){
    double xmin=-2.1, xmax=4.2;
    double x=xmin;
    while(x<=xmax){
        printf("%10g %10g %10g %10g %10g %10g %10g\n",x, erf(x), gsl_sf_erf(x), myErf(x),tgamma(x), gsl_sf_gamma(x), myGamma(x));
        x+=1.0/8;
    }

    return 0;
}