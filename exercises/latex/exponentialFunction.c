#include <math.h>

double exponential(double x)
{
    if (x < 0)
    {
        double result = 1.0 / exponential(-x); //Negative value more stable this way
        return result;
    }
    if (x > 1. / 8)
    {
        double result = pow(exponential(x / 2), 2);
        return result;
    }
    double result = 1 + x * (1 + x / 2 * (1 + x / 3 * (1 + x / 4 * (1 + x / 5 * (1 + x / 6 * (1 + x / 7 * (1 + x / 8 *
                                                                                                               (1 +
                                                                                                                x / 9 *
                                                                                                                (1 + x /
                                                                                                                     10)))))))));
    return result;
}