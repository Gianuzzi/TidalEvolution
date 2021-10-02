# include <math.h>

extern double euler_forward(double t, double y, double dt, derivator dydt)
{
    return y + dt * (dydt)(t, y);
}


extern double euler_backward (double t, double y, double dt, derivator dydt)
{
    double y1 = y;
    double dy1;

    for (unsigned int k = 0; k < MAX_ITER; k++)
    {
        dy1 = (dydt)(t + dt, y1) * dt;
        y1 = y1 + dy1;
        if (fabs ((dydt)(t + dt, y1) * dt) <= MIN_ERR)
        {
            break;
        }
    }

    return y + dy1 * dt;
}


extern double euler_centred (double t, double y, double dt, derivator dydt)
{
    double y1 = y;
    double dy1;

    for (unsigned int k = 0; k < MAX_ITER; k++)
    {
        dy1 = (dydt)(t + dt, y1) * dt;
        y1 = y1 + dy1;
        if (fabs (dy1) <= MIN_ERR)
        {
            break;
        }
    }

    return y + 0.5 * (dy1 + (dydt)(t, y)) * dt;
}


extern double runge_kutta2 (double t, double y, double dt, derivator dydt)
{
    double rk1, rk2;
    double y1;

    rk1 = (dydt)(t, y);
    y1  = y + rk1 * dt;
    rk2 = (dydt)(t + dt, y1);

    return y + dt * 0.5 * (rk1 + rk2);
}


extern double midpoint (double t, double y, double dt, derivator dydt)
{
    double rk1, rk2;
    double y1;

    rk1 = (dydt)(t, y);
    y1  = y + rk1 * dt * 0.5;
    rk2 = (dydt)(t + dt * 0.5, y1);

    return y + rk2 * dt;
}


extern double ralston (double t, double y, double dt, derivator dydt)
{
    double rk1, rk2;
    double y1;

    rk1 = (dydt)(t, y);
    y1  = y + rk1 * dt * 0.75;
    rk2 = (dydt)(t + dt * 0.75, y1);

    return y + (rk1 + 2. * rk2) / 3. * dt;
}


extern double runge_kutta4 (double t, double y, double dt, derivator dydt)
{
    double rk1, rk2, rk3, rk4;
    double y1;

    rk1 = (dydt)(t, y);
    y1  = y + 0.5 * dt * rk1;
    rk2 = (dydt)(t + 0.5 * dt, y1);
    y1  = y + 0.5 * dt * rk2;
    rk3 = (dydt)(t + 0.5 * dt, y1);
    y1  = y + dt * rk3;
    rk4 = (dydt)(t + dt, y1);
    
    return y + (rk1 + 2. * (rk2 + rk3) + rk4) / 6. * dt;
}

