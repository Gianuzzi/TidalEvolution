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
extern void init_params (double *arr,
                         double a1, double e1, double spin1, double oblic1,
                         double spin0, double oblic0)
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

// Return t0, and sets initial conditions from file 
extern double set_from_file (char *filename)
{
    fp = fopen (filename, "r");
    static const long line_length = 96 + 1; // Last is EOF
    char buff[line_length - 1];
    fseek (fp, -line_length, SEEK_END);
    fread (buff, line_length, 1, fp);
    buff[line_length - 1] = '\0'; // NULL terminator
    fclose (fp);

    printf ("%s\n", buff);
    sscanf (buff, "%lf, %lf, %lf, %lf, %lf, %lf, %lf",
                   &a1, &e1, &spin1, &oblic1, &spin0, &oblic0, &t0);
    return t0;
}

