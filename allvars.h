//---------------------- Constants ---------------------------------
// INTEGRATION
static const unsigned int N = 6;                 // Amount of parameters
static const size_t y_size = N * sizeof(double); // Size (bytes) of parameters

// Functions
# define SQUARE(x) (x * x)

// Physics
// # define G 39.478 // AU³ yr⁻² Ms⁻¹
# define PI 3.1415927            // Pi
# define TWO_PI 6.2831853        // Two Pi
# define G SQUARE(0.01720209895) // AU³ days⁻² Ms⁻¹
# define KM2AU(x) (x / 1.496E8)
# define MJ2MS(x) (x * 9.54e-4)
# define MT2MS(x) (x * 3.04e-6)
# define DAY2YR(x) (x / 365.2563)
# define YR2DAY(x) (x * 365.2563)
# define MU (G * (m0 + m1))

// Implicit integrations methods
# define MAX_ITER 10000 // Max iter
# define MIN_ERR 1e-9   // Min err

// User defined constants
unsigned int n_iter; // Amount of iterations
static double dt;    // TimeStep [days]
static double t0;    // Initial itegration time [days]
static double tf;    // Final itegration time  [days]
static double Logt;  // Cte for ratio (t_[i+1] / t_[i])
static int n_points; // Amount of total output points
static char *output;  // Output file name

# define m1prime (m0 * m1) / (m0 + m1)

//---------------------- Dynamics ---------------------------------
// Time calculations
double start_time; // Initial cpu time
double elapsed;    // Total cpu time

FILE *fp;      // File with solution array
char mode[1];  // File with solution array

//---------------------- Input ---------------------------------
unsigned int selection; // Dummy int if output file already exists


//---------------------- Function types ----------------------------
typedef double (*derivator)(double, double *);
typedef double (*integrator)(double, double *, double, int, derivator);

derivator derydt;  // dydt (for testing)
derivator deriva1; // da1dt
derivator derive1; // de1dt
derivator derivs1; // ds1dt
derivator derivo1; // do1dt
derivator derivs0; // ds1dt
derivator derivo0; // do0dt
integrator integ;  // Integrator used

//---------------------- Objects parametes -------------------------

// Object 0
double m0;       // Mass [Ms]
double spin0;    // Spin [rad day⁻¹]
double oblic0;   // Obliquity [rad]
double radius0;  // Radius [AU]
double k20;      // Love coefficient
double dt0;      // Delta t (obliquity) [seg]
double k2dt0;    // Love coefficient * Delta t (obliquity) [seg]
double alpha0;   // Alpha (inertia)

double k0;       // Tidal force magnitude
double c0;       // Moment of inertia
double q0;       // Tidal dissipation factor

// Object 1
double m1;       // Mass [Ms]
double a1;       // Semimajor axis [AU]
double e1;       // Eccentricity
double n1;       // Mean movement [rad seg⁻¹]
double spin1;    // Spin [rad day⁻¹]
double oblic1;   // Obliquity [rad]
double radius1;  // Radius [AU]
double k21;      // Love coefficient
double dt1;      // Delta t (obliquity) [seg]
double k2dt1;    // Love coefficient * Delta t (obliquity) [seg]
double alpha1;   // Alpha (inertia)

double k1;       // Tidal force magnitude
double c1;       // Moment of inertia
double q1;       // Tidal dissipation factor

// ---------------------- Run -------------------------
static double *y; // Params pointer [a1, e1, spin1, oblic1, spin0, oblic0]
static double t0_old; // Time checkpoint for logarithmic scale output

//---------------------------------------------------------------

