# include <math.h>


extern void init_params_zeros (double *arr)
{
    for (unsigned int j = 0; j < N; j++)
    {
        arr[j] = 0.;
    }
}

// 0 ,  1,   2  ,   3   ,   4  ,   5
// a1, e1, spin1, oblic1, spin0, oblic0
extern void init_params (double *arr)
{
    if (N == 6)
    {
        arr[0] = a1;
        arr[1] = e1;
        arr[2] = spin1;
        arr[3] = oblic1;
        arr[4] = spin0;
        arr[5] = oblic0;
    }
    else
    {
        init_params_zeros (arr);       
    }    
}


