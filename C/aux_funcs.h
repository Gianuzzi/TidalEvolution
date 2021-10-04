# include <math.h>

extern double f1_l (double e)
{
    return
    (1. +
    3. * pow (e, 2) +
    0.375 * pow (e, 4)
    ) *
    pow ((1. - pow (e, 2)), -4.5
    );
}

extern double f2_l (double e)
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

extern double f3_l (double e)
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

extern double f4_l (double e)
{
    return 
    (1. +
    1.5 * pow (e, 2) +
    0.125 * pow (e, 4)) *
    pow ((1. - pow (e, 2)), -5
    );
}

extern double f5_l (double e)
{
    return 
    (1. +
    3.75 * pow (e, 2) +
    1.875 * pow (e, 4) +
    0.078125 * pow (e, 6)) *
    pow ((1. - pow (e, 2)), -6.5
    );  
}


extern double f1_s (double e)
{
    return
    (1. + 7.5 * SQUARE (e));
}

extern double f2_s (double e)
{
    return 
    (1. + 13.5 * SQUARE (e));
}

extern double f3_s (double e)
{
    return 
    (1. + 23. * SQUARE (e));
}

extern double f4_s (double e)
{
    return 
    (1. + 6.5 * SQUARE (e));
}

extern double f5_s (double e)
{
    return 
    (1. + 10.25 * SQUARE (e)); 
}

extern void set_f_funcs (char t)
{
    if (t == 0)
    {
        f1 = f1_l;
        f2 = f2_l;
        f3 = f3_l;
        f4 = f4_l;
        f5 = f5_l;
        
    }
    else
    {
        f1 = &f1_s;
        f2 = &f2_s;
        f3 = &f3_s;
        f4 = &f4_s;
        f5 = &f5_s;
    }    
}

extern void eval_f (double e)
{
    fe1 = f1 (e);
    fe2 = f2 (e);
    fe3 = f3 (e);
    fe4 = f4 (e);
    fe5 = f5 (e);
}


// Mean motion
double ni (double a)
{
    return sqrt (MU / pow (a, 3));
}

// k2 * Deltat
extern double k2dti (double q, double a)
{
    return 1.5 / (q * ni (a));
}

// Tidal force magnitude
extern double Ki (double mass, double radius, double q, double a)
{
    return 3. * G * SQUARE (mass) * pow (radius, 5) * k2dti (q, a);
}

// Inertia moment constant
extern double Ci (double mass, double radius, double alpha)
{
    return alpha * mass * SQUARE (radius);
}

// Angular momentum
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
extern void print_n (double t)
{
    printf (" a1:%.3e e1:%.3e spin1:%.3e " \
            "eps1:%.3e spin0:%.3e eps0:%.3e t:%.3e [yr]\n",
            a1, e1, s1, o1, s0, o0, DAY2YR (t));
}