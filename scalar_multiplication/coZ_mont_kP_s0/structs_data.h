#ifndef STRUCTS_DATA
#define STRUCTS_DATA


//~ IMPORTANT : We take 'phi = 1 << WORD_SIZE'
#define WORD_SIZE 64
#define POLY_DEG 5
#define NB_COEFF 6
#define NB_ADD_MAX 38

#define RHO_LOG2 48  // rho = 1 << RHO_LOG2.

typedef __int128 int128;

//~ representatives of polynomials P0 and P1, used for conversion into the AMNS
static int64_t poly_P0[NB_COEFF] = {-6277362495496, 523027163910, -9584902612820, -2137107975958, 18657495145588, 7868340501930};
static int64_t poly_P1[NB_COEFF] = {-29678146215139, -9150175267373, 4805433369377, 18589431413630, 9601662501410, 5304347000407};

//~ representatives of polynomials Pi, for i=2,...,n-1
static int64_t polys_P[(NB_COEFF - 2)][NB_COEFF];

static mpz_t modul_p;
static mpz_t gama_pow[POLY_DEG];

#endif

