#ifndef STRUCTS_DATA
#define STRUCTS_DATA


//~ IMPORTANT : We take 'phi = 1 << WORD_SIZE'
#define WORD_SIZE 64
#define POLY_DEG 5
#define NB_COEFF 6
#define NB_ADD_MAX 5

#define RAND_UP_BOUND 7

#define RHO_LOG2 50  // rho = 1 << RHO_LOG2.

typedef __int128 int128;

//~ representatives of polynomials P0 and P1, used for conversion into the AMNS
static int64_t poly_P0[NB_COEFF] = {1667982572806, 1612493844509, -8764554677309, -12793029824619, -8885025994580, -1464133454158};
static int64_t poly_P1[NB_COEFF] = {-7471337341872, 2362808291243, -5068509920966, -8752072807476, -6109478364357, 2031617687120};

//~ representatives of polynomials Pi, for i=2,...,n-1
static int64_t polys_P[(NB_COEFF - 2)][NB_COEFF];

static mpz_t modul_p;
static mpz_t gama_pow[POLY_DEG];

//~~~~~~~~~~~~~~~~~~~~~~~ random generator stuff ~~~~~~~~~~~~~~~~~~

#define SEEDEXPANDER_MAX_LENGTH             4294967295


AES_XOF_struct aes_seedexpander;

//~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#endif

