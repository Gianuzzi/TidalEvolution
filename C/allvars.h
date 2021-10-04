//---------------------- Constants ---------------------------------

// Functions
# define SQUARE(x) (x * x)

// Physics
// # define G 39.478 // AU³ yr⁻² Ms⁻¹
# define PI 3.1415927            // Pi
# define TWO_PI 6.2831853        // Two Pi
# define G SQUARE(0.01720209895) // AU³ days⁻² Ms⁻¹
# define KM2AU(x) (x * 6.68459e-9)
# define MJ2MS(x) (x * 9.543e-4)
# define MT2MS(x) (x * 3.04e-6)
# define DAY2YR(x) (x / 365.2563)
# define YR2DAY(x) (x * 365.2563)
# define MU (double)(G * (m0 + m1))

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
double start_time;             // Initial cpu time
double elapsed;                // Total cpu time

FILE *fp;                      // Output file
FILE *ip;                      // Output initial params file
char mode[2];                  // Char for possible restart

static double t_out;            // Time checkpoint for logarithmic scale output
static double t_add;            // Time adition for logarithmic scale output

double AM;                      // Angular momentum
double fe1, fe2, fe3, fe4, fe5; // f_i(e1)
unsigned char shortf;           // Use long f_i [=0] or not

//---------------------- Input ---------------------------------
unsigned int selection; // Dummy int if output file already exists

//---------------------- Function types ----------------------------
typedef double (*derivator)(double, double);
typedef double (*integrator)(double, double, double, derivator);
typedef double (*f_func)(double);

derivator derydt;            // dydt (for testing)
derivator deriva1;           // da1dt
derivator derive1;           // de1dt
derivator derivs1;           // ds1dt
derivator derivo1;           // do1dt
derivator derivs0;           // ds1dt
derivator derivo0;           // do0dt
integrator integ;            // Integrator used
f_func f1, f2, f3, f4, f5;   // Pointers to f_i functions

//---------------------- Objects parametes -------------------------

// Object 0
double m0;        // Mass [Ms]
double s0;        // Spin [rad day⁻¹]
double o0;        // Obliquity [rad]
double radius0;   // Radius [AU]
double k20;       // Love coefficient
double dt0;       // Delta t (obliquity) [seg]
double k2dt0;     // Love coefficient * Delta t (obliquity) [seg]
double alpha0;    // Alpha (inertia)

static double k0; // Tidal force magnitude
double c0;        // Moment of inertia
double q0;        // Tidal dissipation factor
double k0cte;     // Tidal force magnitude constant (k0 * a⁻¹/²)

// Object 1
double m1;        // Mass [Ms]
double a1;        // Semimajor axis [AU]
double e1;        // Eccentricity
static double n1; // Mean movement [rad seg⁻¹]
double s1;        // Spin [rad day⁻¹]
double o1;        // Obliquity [rad]
double radius1;   // Radius [AU]
double k21;       // Love coefficient
double dt1;       // Delta t (obliquity) [seg]
double k2dt1;     // Love coefficient * Delta t (obliquity) [seg]
double alpha1;    // Alpha (inertia)

static double k1; // Tidal force magnitude
double c1;        // Moment of inertia
double q1;        // Tidal dissipation factor
double k1cte;     // Tidal force magnitude constant (k1 * a⁻¹/²)

// ---------------------- Run -------------------------
static double t;  // Time


