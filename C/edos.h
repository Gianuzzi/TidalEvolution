# include <math.h>
# include <stdio.h>
# include "aux_funcs.h"

// Test derivate. y^power
extern double dydt (double t, double *y)
{
    double power = 2.;
    return power * pow (t, power - 1);
}

extern double dadt (double t, double *y) 
{   
    double k0aux = Ki(k0cte, y[0]);
    double k1aux = Ki(k0cte, y[0]);
    return
     2. / m1prime * pow (y[0], -7)  * 
     (
        (f2 (y[1]) / ni (y[0])) * (k0aux * cos (y[5]) * y[4] + k1aux * cos (y[3]) * y[2]) -
        f3 (y[1]) * (k0aux + k1aux)
     );
}

extern double dedt (double t, double *y) 
{
    return
     9. / m1prime * pow (y[0], -8) * y[1] *
     (
        (11. / 18. * f4 (y[1]) / ni (y[0])) * (k0 * cos (y[5]) * y[4] + k1 * cos (y[3]) * y[2]) -
        f5 (y[1]) * (k0 + k1)
     );
}

extern double dspin0dt (double t, double *y) 
{
    return
     - k0 * ni(y[0]) / (c0 * pow(y[0], 6)) *
     (
      f1 (y[1]) * 0.5 * (1. + pow (cos (y[5]), 2)) * y[4] / ni (y[0]) - \
      f2 (y[1]) * cos (y[5])
     );
}

extern double dspin1dt (double t, double *y) 
{
    return
     - k1 * ni (y[0]) / (c1 * pow (y[0], 6)) *
     (
      f1 (y[1]) * 0.5 * (1. + pow (cos (y[3]), 2)) * y[2] / ni (y[0]) - \
      f2 (y[1]) * cos (y[3])
     );
}

extern double depsilon0dt (double t, double *y) 
{
    return
     k0 * ni (y[0]) / (c0 * y[4] * pow (y[0], 6)) *
     (
      f1 (y[1]) * 0.5 * y[4] / ni (y[0]) * (cos (y[5]) - (c0 * y[4]) / (m1 * AngMom (y[0], y[1]))) - \
      f2 (y[1])
     ) * sin (y[5]);
}

extern double depsilon1dt (double t, double *y) 
{
    return
     k1 * ni (y[0]) / (c1 * y[2] * pow (y[0], 6)) *
     (
      f1 (y[1]) * 0.5 * y[2] / ni (y[0]) * (cos (y[3]) - (c1 * y[2]) / (m0 * AngMom (y[0], y[1]))) - \
      f2 (y[1])
     ) * sin (y[3]);
}