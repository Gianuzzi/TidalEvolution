# include <math.h>

extern double f1 (double e)
{
    return
    (1. +
    3. * pow (e, 2) +
    0.375 * pow (e, 4)
    ) *
    pow ((1. - pow (e, 2)), -4.5
    );
}

extern double f2 (double e)
{
    return 
    (1. +
    7.5 * pow (e, 2) + 
    5.625 * pow (e, 4) +
    0.3125 * pow (e, 6)
    ) *
    pow ((1. - pow (e, 2)), -6
    ); 
}

extern double f3 (double e)
{
    return 
    (1. +
    15.5 * pow (e, 2) +
    31.875 * pow (e, 4) +
    11.5625 * pow (e, 6) +
    0.390625 * pow (e, 8)) *
    pow ((1. - pow (e, 2)), -7.5
    );
}

extern double f4 (double e)
{
    return 
    (1. +
    1.5 * pow (e, 2) +
    0.125 * pow (e, 4)) *
    pow ((1. - pow (e, 2)), -5
    );
}

extern double f5 (double e)
{
    return 
    (1. +
    3.75 * pow (e, 2) +
    1.875 * pow (e, 4) +
    0.078125 * pow (e, 6)) *
    pow ((1. - pow (e, 2)), -6.5
    );  
}


extern double f1_short (double e)
{
    return
    (1. + 7.5 * SQUARE (e));
}

extern double f2_short (double e)
{
    return 
    (1. + 13.5 * SQUARE (e));
}

extern double f3_short (double e)
{
    return 
    (1. + 23. * SQUARE (e));
}

extern double f4_short (double e)
{
    return 
    (1. + 6.5 * SQUARE (e));
}

extern double f5_short (double e)
{
    return 
    (1. + 10.25 * SQUARE (e)); 
}


// Mean motion
double ni (double a)
{
    return sqrt(MU / pow (a, 3));
}

// k2 * Deltat
extern double k2dti (double q, double a)
{
    return 1.5 / (q * ni (a));
}

// Tidal force magnitude constant (ki * a⁻¹/²)
extern double Kicte (double mass, double radius, double q)
{
    return 4.5 * G * SQUARE (mass) * pow (radius, 5) / (q * sqrt (MU)) ;
}

// Tidal force magnitude
extern double Ki (double kicte, double a)
{
    return kicte * pow (a, 1.5);
}

// // Tidal force magnitude
// extern double Ki (double mass, double radius, double q, double a)
// {
//     return 3. * G * SQUARE (mass) * pow (radius, 5) * k2dti (q, a);
// }


// Inertia moment constant
extern double Ci (double mass, double radius, double alpha)
{
    return alpha * mass * SQUARE (radius);
}

extern double AngMom (double a, double e)
{
    return sqrt(MU * a * (1 - SQUARE (e)));
}

// t_[i+1] / t_[i]
extern double CteLogt (double total_t, int npoints)
{
    return exp (log (total_t) / (npoints - 1)); 
}

// Print parameters of array
extern void print_n (double *y, int i, double t)
{
    if (i < 0)
    {
        if (N == 6)
        {
                    printf (
                        " a1:%.3e e1:%.3e spin1:%.3e " \
                        "eps1:%.3e spin0:%.3e eps0:%.3e t:%.3e [yr]\n",
                        y[0], y[1], y[2], y[3], y[4], y[5], DAY2YR (t)
                   );
        }
        else
        {
            for (unsigned int j = 0; j < N; j++)
            {
                printf (" %.4e,  ", y[j]);
            }
            printf (" t = %.5f [yr]\n", DAY2YR (t));
        }
    }
    else
    {
        printf (" y[%i] = %.5f, t = %.5f [yr]\n", i, y[i], DAY2YR (t));
    }
}