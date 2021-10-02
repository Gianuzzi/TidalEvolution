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
    double y,
    double dt,
    integrator integ
    )
{
    y = integ (t, y, dt, &dydt);
}


// Step and update parameters in tidal interaction
void tidal_step (
    double t,
    double dt,
    integrator integ
    )
{
    double a1n, e1n, s1n, o1n, s0n, o0n;

    n1 = ni (a1);
    k0 = Ki (k0cte, a1);
    k1 = Ki (k1cte, a1);

    a1n = integ (t, a1, dt, deriva1);
    e1n = integ (t, e1, dt, derive1);
    s1n = integ (t, s1, dt, derivs1);
    o1n = integ (t, o1, dt, derivo1);
    s0n = integ (t, s0, dt, derivs0);
    o0n = integ (t, o0, dt, derivo0);

    // #pragma omp parallel sections num_threads(6)
    // {
    //     #pragma omp section
    //         { a1n = integ (t, a1, dt, deriva1);}
    //     #pragma omp section
    //         { e1n = integ (t, e1, dt, derive1);}
    //     #pragma omp section
    //         { s1n = integ (t, s1, dt, derivs1);}
    //     #pragma omp section
    //         { o1n = integ (t, o1, dt, derivo1);}
    //     #pragma omp section
    //         { s0n = integ (t, s0, dt, derivs0);}
    //     #pragma omp section
    //         { o0n = integ (t, o0, dt, derivo0);}
    // }

    a1 = a1n;
    e1 = e1n;
    s1 = s1n;
    o1 = fmod (o1n, TWO_PI);
    s0 = s0n;
    o0 = fmod (o0n, TWO_PI);
}

int main () 
{
    start_time = omp_get_wtime (); // Timestamp

    // Objects
    /// Object 0
    m0      = 1.0;
    s0      = TWO_PI / 28.;
    o0      = 0.0;
    radius0 = KM2AU (695700.);
    alpha0  = 0.4;
    q0      = 1.e6;

    /// Object 1
    m1      = MJ2MS (1.);
    a1      = 0.05;
    e1      = 0.1;
    s1      = TWO_PI / 1.;
    o0      = 0.0;
    radius1 = KM2AU (69911.);
    alpha1  = 0.4;
    q1      = 1.e5;
    
    // Run conditions
    t0       = 0.;
    dt       = YR2DAY (25.);
    tf       = YR2DAY (2.5e9);
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
    print_n(t0);

    // Define output
    output = "Salida.txt";
    strcpy (mode, "w");;

    if (access (output, F_OK) == 0)
    {
        printf ("File %s already exists.\n", output);
        printf ("Enter an option:\n");
        printf ("\t   0: Overwrite file.\n");
        printf ("\t   1: Continue the run, using this physical parameters.\n");
        printf ("\t >=2: Exit.\n");
        scanf ("%i", &selection);
        switch (selection)
        {
        case 0:
            printf ("Overwriting...\n");
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
    print_n (t0);
    if (strcmp(mode, "w"))
    {
        fprintf (fp, "%e, %e, %e, %e, %e, %e, %e\n",
                 a1, e1, s1, o1, s0, o0, t0);
    }  

    // LOOP parameters
    n_iter = (int)((tf - t0) / dt);    
    Logt   = CteLogt (tf - t0, n_points);
    t_add  = Logt;
    printf ("Approximate Iterations: %i\n", n_iter);
    
    /// LOOP
    t      = t0;
    t_out  = t0 + t_add;
    for (unsigned int i = 0; i < n_iter; i++)
    {
        t = t + dt;
        tidal_step (t, dt, integ);
        if (t >= t_out)
        {
            if ((isnan (a1)) || (a1 < 0))
            {
                printf ("End of RUN. [Encounter]\n");
                printf ("Total iterations: %i\n", i);
                break;
            }
            t_add = t_add * Logt;
            t_out = t0 + t_add;
            printf ("Iteration: %i\n", i);
            print_n (t);
            fprintf (fp, "%e, %e, %e, %e, %e, %e, %e\n",
                        a1, e1, s1, o1, s0, o0, t);
        }        
    }

    fclose (fp);

    elapsed = omp_get_wtime () - start_time;
    printf ("Total running time: %f [sec]\n.", elapsed);
return 0;
}


