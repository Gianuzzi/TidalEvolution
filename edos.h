# include <math.h>
# include <stdio.h>
# include "aux_funcs.h"

// Test derivate. y^power
extern double dydt (double t, double *y)
{
    double power = 2.;
    return power * pow(t, power - 1);
}
derivator derydt = &dydt;

extern double dadt (double t, double *y) 
{   
    return
     2. / m1prime * pow(y[0], -7)  * 
     (
        (f2(y[1]) / ni(y[0])) * (k0 * cos(y[5]) * y[4] + k1 * cos(y[3]) * y[2]) -
        f3(y[1]) * (k0 + k1)
     );
}
derivator deriva1 = &dadt;

extern double dedt (double t, double *y) 
{
    return
     9. / m1prime * pow(y[0], -8) * y[1] *
     (
        (11. / 18. * f4(y[1]) / ni(y[0])) * (k0 * cos(y[5]) * y[4] + k1 * cos(y[3]) * y[2]) -
        f5(y[1]) * (k0 + k1)
     );
}
derivator derive1 = &dedt;

extern double dspin0dt (double t, double *y) 
{
    return
     - k0 * ni(y[0]) / (c0 * pow(y[0], 6)) *
     (
      f1(y[1]) * 0.5 * (1. + pow(cos(y[5]), 2)) * y[4] / ni(y[0]) - \
      f2(y[1]) * cos(y[5])
     );
}
derivator derivs0 = &dspin0dt;

extern double dspin1dt (double t, double *y) 
{
    return
     - k1 * ni(y[0]) / (c1 * pow(y[0], 6)) *
     (
      f1(y[1]) * 0.5 * (1. + pow(cos(y[3]), 2)) * y[2] / ni(y[0]) - \
      f2(y[1]) * cos(y[3])
     );
}
derivator derivs1 = &dspin1dt;

extern double depsilon0dt (double t, double *y) 
{
    return
     k0 * ni(y[0]) / (c0 * y[4] * pow(y[0], 6)) *
     (
      f1(y[1]) * 0.5 * y[4] / ni(y[0]) * (cos(y[5]) - (c0 * y[4]) / (m1 * AngMom(y[0], y[1]))) - \
      f2(y[1])
     ) * sin(y[5]);
}
derivator derivo0 = &depsilon0dt;

extern double depsilon1dt (double t, double *y) 
{
    return
     k1 * ni(y[0]) / (c1 * y[2] * pow(y[0], 6)) *
     (
      f1(y[1]) * 0.5 * y[2] / ni(y[0]) * (cos(y[3]) - (c1 * y[2]) / (m0 * AngMom(y[0], y[1]))) - \
      f2(y[1])
     ) * sin(y[3]);
}
derivator derivo1 = &depsilon1dt;




// extern double da2dt (double t, double *y) 
// {
//     return 2. / (m1prime * pow(B1_ptr->a, 7)) + auxdadt(B0_ptr) + auxdadt(B1_ptr);
// }

// extern double de2dt (double t, double *y) 
// {
//     return 9. * B1_ptr->a / (m1prime * pow(B1_ptr->a, 8)) + auxdedt(B0_ptr) + auxdedt(B1_ptr);
// }

// extern double dspindt (double t, double *y) 
// {
//     return - 9. * B1_ptr->a / (m1prime * pow(B1_ptr->a, 8)) + auxdedt(B0_ptr) + auxdedt(B1_ptr);
// }