# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <omp.h>
# include <string.h>
# include <unistd.h>
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
    derydt = &dydt;
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
    double a1n, e1n, s1n, o1n, s0n, o0n;

    k0 = Ki(k0cte, y[0]);
    k1 = Ki(k1cte, y[0]);

    a1n = integ (t, y, dt, 0, deriva1);
    e1n = integ (t, y, dt, 1, derive1);
    s1n = integ (t, y, dt, 2, derivs1);
    o1n = integ (t, y, dt, 3, derivo1);
    s0n = integ (t, y, dt, 4, derivs0);
    o0n = integ (t, y, dt, 5, derivo0);

    // #pragma omp parallel sections num_threads(6)
    // {
    //     #pragma omp section
    //         { a1n = integ (t, y, dt, 0, deriva1);}
    //     #pragma omp section
    //         { e1n = integ (t, y, dt, 1, derive1);}
    //     #pragma omp section
    //         { s1n = integ (t, y, dt, 2, derivs1);}
    //     #pragma omp section
    //         { o1n = integ (t, y, dt, 3, derivo1);}
    //     #pragma omp section
    //         { s0n = integ (t, y, dt, 3, derivs0);}
    //     #pragma omp section
    //         { o0n = integ (t, y, dt, 5, derivo0);}
    // }

    y[0] = a1n;
    y[1] = e1n;
    y[2] = fmod (s1n, TWO_PI);
    y[3] = o1n;
    y[4] = fmod (s0n, TWO_PI);
    y[5] = o0n;
}

int main () 
{
    start_time = omp_get_wtime (); // Timestamp

    // Objects
    /// Object 0
    m0      = 1.0;
    spin0   = TWO_PI / 28.;
    oblic0  = 0.0;
    radius0 = KM2AU (695700.);
    alpha0  = 0.4;
    q0      = 1.e6;

    /// Object 1
    m1      = MJ2MS (1.);
    a1      = 0.05;
    e1      = 0.1;
    spin1   = TWO_PI / 1.;
    oblic1  = 0.0;
    radius1 = KM2AU (69911.);
    alpha1  = 0.4;
    q1      = 1.e5;
    
    // Run conditions
    t0       = 0.;
    dt       = YR2DAY (25.);
    tf       = YR2DAY (2.1e9);
    n_points = 3000;

    /// Calculated
    n1    = ni (a1);
    k0cte = Kicte (m1, radius0, q0);
    k1cte = Kicte (m0, radius1, q1);
    c0    = Ci (m0, radius0, alpha0);
    c1    = Ci (m1, radius1, alpha1);
    k0    = Ki (k0cte, a1);
    k1    = Ki (k1cte, a1);

    // Define Equations
    integ   = &runge_kutta4;
    deriva1 = &dadt;
    derive1 = &dedt;
    derivs0 = &dspin0dt;
    derivs1 = &dspin1dt;
    derivo0 = &depsilon0dt;
    derivo1 = &depsilon1dt;

    
    // Init memory
    y = malloc (y_size);

    // Define output
    output  = "Salida.txt";
    strcpy (mode, "w");

    if (access (output, F_OK) == 0)
    {
        printf ("File %s already exists.\n", output);
        printf ("Enter an option:\n");
        printf ("\t   0: Rewrite file.\n");
        printf ("\t   1: Continue the run, using this physical parameters.\n");
        printf ("\t >=2: Exit.\n");
        scanf ("%i", &selection);
        switch (selection)
        {
        case 0:
            printf ("Rewriting...\n");
            break;
        
        case 1:
            printf ("Reading parameters from file...\n");
            t0 = set_from_file(output);
            strcpy (mode, "a");
            break;
        
        default:
            printf ("Exiting.\n");
            exit (0);
            break;
        }
    }

    /// Check times
    if (tf < t0)
    {
        printf ("Final time %f [days] already reached.\n", tf);
        printf ("Exiting.\n");
        exit (0);
    }

    // Open file
    fp = fopen (output, mode);
    
    /// Set and Print parameters array
    init_params (y, a1, e1, spin1, oblic1, spin0, oblic0);
    print_n (y, -1, t0);
    if (strcmp(mode, "w"))
    {
        fprintf (fp, "%e, %e, %e, %e, %e, %e, %e\n",
            y[0], y[1], y[2], y[3], y[4], y[5], t0
        );
    }  

    // LOOP parameters
    n_iter = (int)((tf - t0) / dt);    
    Logt   = CteLogt (tf - t0, n_points);
    t_add  = Logt;
    printf("Approximate Iterations: %i\n", n_iter);
    
    /// LOOP
    t      = t0;
    t_out  = t0 + t_add;
    for (unsigned int i = 0; i < n_iter; i++)
    {
        t = t + dt;
        tidal_step (t, y, dt, integ);
        if (t >= t_out)
        {
            if ((isnan (y[0])) || (y[0] < 0))
            {
                printf ("End of RUN. [Encounter]\n");
                printf ("Total iterations: %i\n", i);
                break;
            }
            t_add = t_add * Logt;
            t_out = t0 + t_add;
            printf ("Iteration: %i\n", i);
            print_n (y, -1, t);
            fprintf (fp, "%e, %e, %e, %e, %e, %e, %e\n",
                        y[0], y[1], y[2], y[3], y[4], y[5], t
            );
        }        
    }

    fclose (fp);    
    
    free (y);

    elapsed = omp_get_wtime () - start_time;
    printf ("Total running time: %f [sec]\n.", elapsed);
return 0;
}


