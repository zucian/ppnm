#include <math.h>
double myGamma(double x){
///single precision gamma function (Gergo Nemes, from Wikipedia)
    if(x<0)return M_PI/sin(M_PI*x)/myGamma(1-x);
    if(x<9)return myGamma(x+1)/x;
    double lnGamma=x*log(x+1/(12*x-1/x/10))-x+log(2*M_PI/x)/2;
    return exp(lnGamma);
}