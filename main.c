# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <omp.h>
# include <string.h>
# include "allvars.h"
# include "init.h"
# include "integrators.h"
# include "edos.h"

// Make step of polomial to test integrators.
void test_step (
    double t,
    double *y,
    double dt,
    integrator integ
    )
{
    derivator derydt = &dydt;
    for (unsigned int j = 0; j < N; j++)
    {
        y[j] = integ (t, y, dt, j, derydt);
    }
}


// Step and update parameters in tidal interaction
void tidal_step (
    double t,
    double *y,
    double dt,
    integrator integ
    )
{
    derivator deriva1 = &dadt;
    derivator derive1 = &dedt;
    derivator derivs1 = &dspin1dt;
    derivator derivo1 = &depsilon1dt;    
    derivator derivs0 = &dspin0dt;
    derivator derivo0 = &depsilon0dt;
    double a1n, e1n, s1n, o1n, s0n, o0n;

    a1n = integ (t, y, dt, 0, deriva1);
    e1n = integ (t, y, dt, 1, derive1);
    s1n = integ (t, y, dt, 2, derivs1);
    if (y[3] == 0.)
    {
        o1n = 0.;
    }
    else
    {
        o1n = integ (t, y, dt, 3, derivo1);
    }  
    s0n = integ (t, y, dt, 4, derivs0);
    if (y[5] == 0.)
    {
        o0n = 0.;
    }
    else
    {
        o0n = integ (t, y, dt, 5, derivo0);
    }  

    y[0] = a1n;
    y[1] = e1n;
    y[2] = fmod(s1n, TWO_PI);
    y[3] = o1n;
    y[4] = fmod(s0n, TWO_PI);
    y[5] = o0n;
}

int main () 
{
    start_time = omp_get_wtime (); // Timestamp

    // Objects
    /// Object 0
    m0      = 1.0f;
    spin0   = TWO_PI / 28.;
    oblic0  = 0.0;
    radius0 = KM2AU(695700.);
    alpha0  = 0.4;
    q0      = 1.e6;

    /// Object 1
    m1      = MJ2MS(1.);
    a1      = 0.05;
    e1      = 0.1;
    spin1   = TWO_PI / 1.;
    oblic1  = 0.0;
    radius1 = KM2AU(69911.);
    alpha1  = 0.4;
    q1      = 1.e5;

    /// Calculated
    n1    = ni(a1);
    k2dt0 = k2dti(q0, n1);
    k2dt1 = k2dti(q1, n1);
    k0    = Ki(m1, radius0, k2dt0);
    c0    = Ci(m0, radius0, alpha0);
    k1    = Ki(m0, radius1, k2dt1);
    c1    = Ci(m1, radius1, alpha1);

    // Run conditions
    t0       = DAY2YR(0.0);
    dt       = DAY2YR(1e5);
    tf       = DAY2YR(1e15);
    n_points = 2000;
    t0_old   = 0.;

    n_iter = (int)((tf - t0) / dt);
    Logt   = CteLogt(tf - t0, n_points);

    printf("ITER: %i\n", n_iter);

    // RUN
    y = malloc (y_size);
    init_params (y);
    print_n (y, -1, t0);

    integrator integ = &runge_kutta4;

    fp = fopen ("salida2.txt", "w");
    fprintf (fp, "%e, %e, %e, %e, %e, %e, %e\n",
             y[0], y[1], y[2], y[3], y[4], y[5], t0/365.25
            );
        
    for (unsigned int i = 0; i < n_iter; i++)
    {
        t0 = t0 + dt;
        tidal_step (t0, y, dt, integ);
        // if (fmod (i, 5000) == 0)
        if (t0 >= (t0_old * Logt))
        {
            if (isnan (y[0]))
            {
                printf ("End of RUN. [Encounter]\n");
                break;
            }
            t0_old = t0;
            print_n (y, -1, t0);
            fprintf (fp, "%e, %e, %e, %e, %e, %e, %e\n",
                        y[0], y[1], y[2], y[3], y[4], y[5], t0/365.25
                    );
        }        
    }

    fclose(fp);
    
    free (y);

    elapsed = omp_get_wtime () - start_time;
    printf("Tiempo total: %f\n.", elapsed);
return 0;
}


