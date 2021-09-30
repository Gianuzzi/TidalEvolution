//---------------------- Constants ---------------------------------
// INTEGRATION
static const unsigned int N = 6; // Amount of parameters
static const size_t y_size = N * sizeof(double); // Size (bytes) of parameters

// Functions
# define SQUARE(x) (x * x)

// Physics
// # define G 39.478f // AU³ yr⁻² Ms⁻¹
# define PI 3.1415927 // Pi
# define TWO_PI 6.2831853 // Two Pi
# define G SQUARE(0.01720209895) // AU³ days⁻² Ms⁻¹
# define KM2AU(x) (x / 1.496E8)
# define MJ2MS(x) (x * 9.54e-4)
# define MT2MS(x) (x * 3.04e-6)
# define DAY2YR(x) (x / 365.2563)
# define MU (G * (m0 + m1))

// Implicit integrations methods
# define MAX_ITER 10000 // Max iter
# define MIN_ERR 1e-9 // Min err

// User defined constants
unsigned int n_iter; // Amount of iterations
static double dt; // TimeStep [days]
static double t0; // Initial itegration time
static double tf; // Final itegration time
static double Logt; // Cte for ratio (t_[i+1] / t_[i])
static int n_points; // Amount of total output points

# define m1prime (m0 * m1) / (m0 + m1)
// # define m1prima (Body0.mass*Body1.mass)/(Body0.mass+Body1.mass)



//---------------------- Dynamics ---------------------------------
// Time calculations
double start_time; // Initial cpu time
double elapsed; // Total cpu time

FILE *fp; // File with solution array


//---------------------- Function types ----------------------------
typedef double (*derivator)(double, double *);
typedef double (*integrator)(double, double *, double, int, derivator);


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
double *y; // Params pointer [a1, e1, spin1, oblic1, spin0, oblic0]
double t0_old; // Time checkpoint for logarithmic scale output

//---------------------------------------------------------------
/// STRUCT TYPE

// mass, a, e, i, M, w, O, n, spin, epsilon, radius, alpha, k2, dt
typedef struct Object
{
    const double mass;
    double a;
    double e;
    double i;
    double M;
    double w;
    double O;
    double n;
    double spin;
    double epsilon;
    const double radius;
    const double alpha;
    const double k2;
    double dt;
} Object, *Object_ptr;

Object Body0;
Object Body1;

Object_ptr B0_ptr = &Body0;
Object_ptr B1_ptr = &Body1;

// struct TidAll
// {
//     struct object Body0;
//     struct object Body1;
// } TidAll, *TidAll_ptr;
