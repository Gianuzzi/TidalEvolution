# include <math.h>
# include "allvars.h"

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
                   &a1, &e1, &s1, &o1, &s0, &o0, &t0);
    return t0;
}

