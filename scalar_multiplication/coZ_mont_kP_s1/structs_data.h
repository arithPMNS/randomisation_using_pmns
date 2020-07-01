#ifndef STRUCTS_DATA
#define STRUCTS_DATA


//~ IMPORTANT : We take 'phi = 1 << WORD_SIZE'
#define WORD_SIZE 64
#define POLY_DEG 5
#define NB_COEFF 6
#define NB_ADD_MAX 18

#define RAND_UP_BOUND 7

#define RHO_LOG2 50  // rho = 1 << RHO_LOG2.

typedef __int128 int128;

//~ representatives of polynomials P0 and P1, used for conversion into the AMNS
static int64_t poly_P0[NB_COEFF] = {14341040004112, 12823019907127, 17071019866162, -283208545568, -12399058895610, -4969842846402};
static int64_t poly_P1[NB_COEFF] = {19688664322726, 24755155911171, 19406089238075, 6677515902099, -14139160744811, -5933139979368};

//~ representatives of polynomials Pi, for i=2,...,n-1
static int64_t polys_P[(NB_COEFF - 2)][NB_COEFF];

static mpz_t modul_p;
static mpz_t gama_pow[POLY_DEG];

//~~~~~~~~~~~~~~~~~~~~~~~ random generator stuff ~~~~~~~~~~~~~~~~~~

#define SEEDEXPANDER_MAX_LENGTH             4294967295


AES_XOF_struct aes_seedexpander;

//~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#endif

