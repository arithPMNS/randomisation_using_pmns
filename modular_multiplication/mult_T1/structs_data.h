#ifndef STRUCTS_DATA
#define STRUCTS_DATA


//~ IMPORTANT : We take 'phi = 1 << WORD_SIZE'
#define WORD_SIZE 64
#define POLY_DEG 5
#define NB_COEFF 6
#define NB_ADD_MAX 4

#define RAND_UP_BOUND 7

#define RHO_LOG2 50  // rho = 1 << RHO_LOG2.

#define BETA_LOG2 47  // beta = 1 << BETA_LOG2, for conversion

typedef __int128 int128;

//~ representations of polynomials Pi, used for conversion into the AMNS
//~ Note: each Pi is a representation of ((1 << BETA_LOG2)^i * phi^2)
int64_t polys_P[NB_COEFF][NB_COEFF] = {
	{1234317683957, 4003804707729, 7226159628765, 890089753249, 2908109985955, 2399745195435},
	{3346495884018, 6951736482028, 3692172826646, -828086939081, 4623697724695, 1055922684655},
	{1924932799176, 9684348167139, 445321944744, 4336436773335, 7511079625150, 2140082962089},
	{2939433148074, 1145712344437, 2206485039221, 6315832225441, 1540592708651, 2651300698217},
	{10834769720226, 2672495267424, 6878085499697, 5974187531900, 4585109675848, 1783165917386},
	{5308476785205, 6322768853035, -2513504629157, 3149487371550, 4832150635938, 2654267372012}};

mpz_t modul_p;
mpz_t gama_pow[POLY_DEG];

//~~~~~~~~~~~~~~~~~~~~~~~ random generator stuff ~~~~~~~~~~~~~~~~~~

#define SEEDEXPANDER_MAX_LENGTH             4294967295


AES_XOF_struct aes_seedexpander;

//~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#endif

