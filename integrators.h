# include <math.h>


extern double euler_forward(double t, double *y, double dt, int i, derivator dydt)
{
    return y[i] + dt * (dydt)(t, y);
}


extern double euler_backward (double t, double *y, double dt, int i, derivator dydt)
{
    double *y1 = (double*) malloc (y_size); // Aux array 
    memcpy (y1, y, y_size);

    for (unsigned int k = 0; k < MAX_ITER; k++)
    {
        y1[i] = y1[i] + (dydt)(t + dt, y1) * dt;
        if ((dydt)(t + dt, y1) * dt <= MIN_ERR)
        {
            break;
        }
    }

    double ans = y[i] + (dydt)(t + dt, y1) * dt;
    free (y1);

    return ans;
}


extern double euler_centred (double t, double *y, double dt, int i, derivator dydt)
{
    double *y1 = (double*) malloc (y_size); // Aux array 
    memcpy (y1, y, y_size);

    for (unsigned int k = 0; k < MAX_ITER; k++)
    {
        y1[i] = y1[i] + (dydt)(t + dt, y1) * dt;
        if ((dydt)(t + dt, y1) * dt <= MIN_ERR)
        {
            break;
        }
    }

    double ans = y[i] + 0.5 * ((dydt)(t + dt, y1) + (dydt)(t, y)) * dt;
    free (y1);

    return ans;
}


extern double runge_kutta2 (double t, double *y, double dt, int i, derivator dydt)
{
    double rk1, rk2;
    double *y1 = (double*) malloc (y_size); // Aux array 
    memcpy (y1, y, y_size);

    rk1   = (dydt)(t, y);
    y1[i] = y[i] * rk1 * dt;
    rk2   = (dydt)(t + dt, y1);

    free (y1);

    return y[i] + dt * 0.5 * (rk1 + rk2);
}


extern double midpoint (double t, double *y, double dt, int i, derivator dydt)
{
    double rk1, rk2;
    double *y1 = (double*) malloc (y_size); // Aux array 
    memcpy (y1, y, y_size);

    rk1   = (dydt)(t, y);
    y1[i] = y[i] * rk1 * dt * 0.5;
    rk2   = (dydt)(t + dt * 0.5, y1);

    free (y1);

    return y[i] + rk2 * dt;
}


extern double ralston (double t, double *y, double dt, int i, derivator dydt)
{
    double rk1, rk2;
    double *y1 = (double*) malloc (y_size); // Aux array 
    memcpy (y1, y, y_size);

    rk1   = (dydt)(t, y);
    y1[i] = y[i] * rk1 * dt * 0.75;
    rk2   = (dydt)(t + dt * 0.75, y1);

    free (y1);

    return y[i] + (rk1 + 2. * rk2) / 3. * dt;
}


extern double runge_kutta4 (double t, double *y, double dt, int i, derivator dydt)
{
    double rk1, rk2, rk3, rk4;
    double *y1 = (double*) malloc (y_size); // Aux array 
    memcpy (y1, y, y_size);

    rk1   = (dydt)(t, y);
    y1[i] = y[i] + 0.5 * dt * rk1;
    rk2   = (dydt)(t + 0.5 * dt, y1);
    y1[i] = y[i] + 0.5 * dt * rk2;
    rk3   = (dydt)(t + 0.5 * dt, y1);
    y1[i] = y[i] + dt * rk3;
    rk4   = (dydt)(t + dt, y1);

    free (y1);
    
    return y[i] + (rk1 + 2. * (rk2 + rk3) + rk4) / 6. * dt;
}

