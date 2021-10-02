# include <math.h>
# include <stdio.h>
# include "aux_funcs.h"
# include "allvars.h"

// Test derivate. y^power
extern double dydt (double t, double y)
{
    double power = 2.;
    return power * pow (t, power - 1);
}

extern double dadt (double t, double a) 
{   
    double k0aux = Ki(k0cte, a);
    double k1aux = Ki(k0cte, a);
    return
     2. / m1prime * pow (a, -7)  * 
     (
        (f2 (e1) / ni (a)) * (k0aux * cos (o0) * s0 + k1aux * cos (o1) * s1) -
        f3 (e1) * (k0aux + k1aux)
     );
}

extern double dedt (double t, double e) 
{
    return
     9. / m1prime * pow (a1, -8) * e *
     (
        (11. / 18. * f4 (e) / n1) * (k0 * cos (o0) * s0 + k1 * cos (o1) * s1) -
        f5 (e) * (k0 + k1)
     );
}

extern double dspin0dt (double t, double s) 
{
    return
     - k0 * n1 / (c0 * pow(a1, 6)) *
     (
      f1 (e1) * 0.5 * (1. + pow (cos (o0), 2)) * s / n1 - \
      f2 (e1) * cos (o0)
     );
}

extern double dspin1dt (double t, double s) 
{
    return
     - k1 * n1 / (c1 * pow (a1, 6)) *
     (
      f1 (e1) * 0.5 * (1. + pow (cos (o1), 2)) * s / n1 - \
      f2 (e1) * cos (o1)
     );
}

extern double depsilon0dt (double t, double o) 
{
    return
     k0 * n1 / (c0 * s0 * pow (a1, 6)) *
     (
      f1 (e1) * 0.5 * s0 / n1 * (cos (o) - (c0 * s0) / (m1 * AngMom (a1, e1))) - \
      f2 (e1)
     ) * sin (o);
}

extern double depsilon1dt (double t, double o) 
{
    return
     k1 * n1 / (c1 * s1 * pow (a1, 6)) *
     (
      f1 (e1) * 0.5 * s1 / n1 * (cos (o) - (c1 * s1) / (m0 * AngMom (a1, e1))) - \
      f2 (e1)
     ) * sin (o);
}