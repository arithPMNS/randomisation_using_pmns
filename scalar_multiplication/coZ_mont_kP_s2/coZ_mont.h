#ifndef CO_Z_MONT
#define CO_Z_MONT 


#define SCAL_SIZE 256 // 'order_of_P' bit-size


//~ note : if not part of a 'CoZ_Points' object, then the Z coordinate is equal to 1, i.e. affines coordinates
typedef struct Inter_Point{
	int64_t X[NB_COEFF];
	int64_t Y[NB_COEFF];
} Inter_Point;

typedef struct CoZ_Points{
	Inter_Point R[2];
	int64_t Z[NB_COEFF];
} CoZ_Points;

typedef struct JacPoint{
	int64_t X[NB_COEFF];
	int64_t Y[NB_COEFF];
	int64_t Z[NB_COEFF];
} JacPoint;

typedef struct affPoint{
	mpz_t x;
	mpz_t y;
} affPoint;


mpz_t const_ONE; // will be used to compute (random) representation of 1 in the AMNS; note : 1 is 'phi' in the AMNS


//~ curve parameters
mpz_t crv_a;
mpz_t crv_b;
int64_t curve_a[NB_COEFF];
int64_t curve_b[NB_COEFF];


//~ base point
affPoint aff_P;
mpz_t order_of_P;


void init_curve_and_base_point(int64_t *rand_pol);

void free_point_data();


int is_on_curve_aff(affPoint *ap);

int is_on_curve_jac(JacPoint *jp);

void print_point_in_aff(affPoint *ap);

void print_point_in_jac(JacPoint *jp);


void init_affPoint(affPoint *ap);

void clear_affPoint(affPoint *ap);


void from_aff_to_jacOneZ_toMontDomain(Inter_Point *jp, affPoint *ap, int64_t *rand_pol);

void from_aff_to_jac_toMontDomain(JacPoint *jp, affPoint *ap, int64_t *rand_pol);

void from_jac_to_aff_outOfMontDomain(affPoint *ap, JacPoint *jp);


//~ -------------------------------------------------------------------------------------

//~ IMPORTANT : two exact_coeffs_reduction are done here to reduce the required number of 'free additions' from 12 to 3
//~ Action: (R, S) <- (2*P, P), with the same Z
void dblu_with_free_adds__fixed_randPoly(int64_t *Z_out, Inter_Point *R, Inter_Point *S, Inter_Point *P, int64_t *rand_pol);

//~ Action: (R, S) = (P+Q, P) <- (P, Q), with the same R.Z = S.Z = Z_out, where P.Z = Q.Z = Z_in
//~ Note: the sequence of instructions allows inputs and outputs to be the same.
void zaddu_with_free_adds__fixed_randPoly(Inter_Point *R, Inter_Point *S, Inter_Point *P, Inter_Point *Q, int64_t *Z_in, int64_t *Z_out, int64_t *rand_pol); 

//~ Action: (R, S) = (P+Q, P-Q) <- (P, Q), with the same R.Z = S.Z = Z_out, where P.Z = Q.Z = Z_in
//~ Note: the sequence of instructions allows inputs and outputs to be the same.
void zaddc_with_free_adds__fixed_randPoly(Inter_Point *R, Inter_Point *S, Inter_Point *P, Inter_Point *Q, int64_t *Z_in, int64_t *Z_out, int64_t *rand_pol); 

void coZ_montgomery_scal_mult_with_free_adds__fixed_randPoly(JacPoint *rop_Q, Inter_Point *op_P, char *scal_bits, int64_t *rand_pol);



#endif





















