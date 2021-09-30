# include <math.h>

extern double f1 (double e)
{
    return 
    (
        1. +
        3. * pow(e, 2) +
        0.375 * pow(e, 4)
        ) *
        pow((1. - pow(e, 2)), -4.5);
}

extern double f2 (double e)
{
    return 
    (
        1. + 
        7.5 * pow(e, 2) + 
        5.625 * pow(e, 4) +
        0.3125 * pow(e, 6)
        ) *
        pow((1. - pow(e, 2)), -6);  
}

extern double f3 (double e)
{
    return 
    (
        1. +
        15.5 * pow(e, 2) +
        31.875 * pow(e, 4) +
        11.5625 * pow(e, 6) +
        0.390625 * pow(e, 8)) *
        pow((1. - pow(e, 2)), -7.5);
}

extern double f4 (double e)
{
    return 
    (
        1. + 
        1.5 * pow(e, 2) +
        0.125 * pow(e, 4)) *
        pow((1. - pow(e, 2)), -5);
}

extern double f5 (double e)
{
    return 
    (
        1. +
        3.75 * pow(e, 2) +
        1.875 * pow(e, 4) +
        0.078125 * pow(e, 6)) *
        pow((1. - pow(e, 2)), -6.5);  
}

// Mean motion
double ni (double a)
{
    return sqrt(MU / pow(a, 3));
}

// Tidal force magnitude
extern double Ki (double mass, double radius, double k2dt)
{
    return 3.f * G * SQUARE(mass) * pow(radius, 5) * k2dt;
}

// Inertia moment constant
extern double Ci (double mass, double radius, double alpha)
{
    return alpha * mass * SQUARE(radius);
}

extern double k2dti (double q, double n)
{
    return 1.5 / (q * n);
}

extern double AngMom (double a, double e)
{
    return sqrt(MU * a * (1 - SQUARE(e)));
}

// t_[i+1] / t_[i]
extern double CteLogt (double total_t, int npoints)
{
    return exp(log(total_t) / npoints); 
}

// t[1] -> 100
// --> t[2] = cte * t[1]
// --> t[3] = cte * t[2] = cte^2 * t[1] 

// STRUCT MODE
// # include "allvars.h"

// extern double auxdadt(double *y, int i)
// {
//     if (i == 0)
//     {
//         return k0 * (f2(B1_ptr->e) * cos(obj->epsilon) * obj->spin / B1_ptr->n - f3(B1_ptr->e));
//     }
// }

// extern double auxdedt(Object_ptr obj)
// {
//     (11. / 18.f * f4(B1_ptr->e) * cos(obj->epsilon) * obj->spin / B1_ptr->n - f5(B1_ptr->e));
// }


// // Tidal force magnitude
// extern double Ki (Object_ptr obj)
// {
//     return 3.f * G * pow(obj->mass, 2), pow(obj->radius, 5) * obj->k2 * obj->dt;
// }

// // Inertia moment constant
// extern double Ci (Object_ptr obj)
// {
//     return obj->alpha * obj->mass * pow(obj->radius, 2);
// }

// extern double auxdadt(Object_ptr obj)
// {
//     return Ki(obj) * (f2(B1_ptr->e) * cos(obj->epsilon) * obj->spin / B1_ptr->n - f3(B1_ptr->e));
// }

// extern double auxdedt(Object_ptr obj)
// {
//     return Ki(obj) * (11. / 18.f * f4(B1_ptr->e) * cos(obj->epsilon) * obj->spin / B1_ptr->n - f5(B1_ptr->e));
// }

// Print parameters of array
extern void print_n (double *arr, int i, double t)
{
    if (i < 0)
    {
        if (N == 6)
        {
                    printf 
                    (
                        "a1:%.3e e1:%.3e spin1:%.3e " \
                        "eps1:%.3e spin0:%.3e eps0:%.3e t:%.3e [yr]\n",
                        y[0], y[1], y[2], y[3], y[4], y[5], t0/365.25
                   );
        }
        else
        {
            for (unsigned int j = 0; j < N; j++)
            {
                printf("%.4e,  ", arr[j]);
            }
            printf("t = %.5f [yr]\n", t/365.25);
        }
    }
    else
    {
        printf("y[%i] = %.5f, t = %.5f [yr]\n", i, arr[i], t/365.25);
    }
}