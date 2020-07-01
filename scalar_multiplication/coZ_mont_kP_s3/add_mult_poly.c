#include "add_mult_poly.h"


void add_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] + pb[j];
}

//~ computes : pa + 2.pb
void double_add_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] + 2*pb[j];
}

void sub_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] - pb[j];
}

void neg_poly(int64_t *rop, int64_t *op){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = -op[j];
}

//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.
void scalar_mult_poly(int64_t *rop, int64_t *op, int64_t scalar){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = scalar * op[j];
}

//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.
void double_poly_coeffs(int64_t *rop, int64_t *op){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = op[j] << 1;
}

//~ assumes 'nb_pos' and/or coeffs of 'op' small enough to avoid an overflow.
void lshift_poly_coeffs(int64_t *rop, int64_t *op, int nb_pos){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = op[j] << nb_pos;
}

//~ computes: J = Z.M[E]
//~ WARNING: 'rand_zero' and 'rand_pol' must be different
void compute_rand_zero(int64_t *rand_zero, int64_t *rand_pol){

	rand_zero[0] = rand_pol[1] * -14477579451936 + rand_pol[2] * -3417025446624 + rand_pol[0] * -958314944803 + rand_pol[4] * 15312819876900 + rand_pol[5] * -4383068741696 + rand_pol[3] * 14245031534996;
	rand_zero[1] = rand_pol[2] * -14477579451936 + rand_pol[3] * -3417025446624 + rand_pol[1] * -958314944803 + rand_pol[5] * 15312819876900 + rand_pol[0] * 1095767185424 + rand_pol[4] * 14245031534996;
	rand_zero[2] = rand_pol[4] * -3417025446624 + rand_pol[3] * -14477579451936 + rand_pol[2] * -958314944803 + rand_pol[0] * -3828204969225 + rand_pol[1] * 1095767185424 + rand_pol[5] * 14245031534996;
	rand_zero[3] = rand_pol[5] * -3417025446624 + rand_pol[4] * -14477579451936 + rand_pol[3] * -958314944803 + rand_pol[0] * -3561257883749 + rand_pol[1] * -3828204969225 + rand_pol[2] * 1095767185424;
	rand_zero[4] = rand_pol[5] * -14477579451936 + rand_pol[4] * -958314944803 + rand_pol[1] * -3561257883749 + rand_pol[2] * -3828204969225 + rand_pol[3] * 1095767185424 + rand_pol[0] * 854256361656;
	rand_zero[5] = rand_pol[5] * -958314944803 + rand_pol[2] * -3561257883749 + rand_pol[0] * 3619394862984 + rand_pol[3] * -3828204969225 + rand_pol[4] * 1095767185424 + rand_pol[1] * 854256361656;
}

//~ random polynomial generator.
//~ IMPORTANT : we assume that : NB_COEFF = 6, RAND_UP_BOUND = 7
void gen_rand_pol(int64_t *rand_pol){
	
	unsigned char random_data[NB_COEFF];
	
	seedexpander(&aes_seedexpander, random_data, NB_COEFF);
	
	
	//~ only for: NB_COEFF = 6, RAND_UP_BOUND = 7 = 2^3 - 1
	rand_pol[0] = (random_data[0] & RAND_UP_BOUND) - ((random_data[0] >> 3) & RAND_UP_BOUND);
	rand_pol[1] = (random_data[1] & RAND_UP_BOUND) - ((random_data[1] >> 3) & RAND_UP_BOUND);
	rand_pol[2] = (random_data[2] & RAND_UP_BOUND) - ((random_data[2] >> 3) & RAND_UP_BOUND);
	rand_pol[3] = (random_data[3] & RAND_UP_BOUND) - ((random_data[3] >> 3) & RAND_UP_BOUND);
	rand_pol[4] = (random_data[4] & RAND_UP_BOUND) - ((random_data[4] >> 3) & RAND_UP_BOUND);
	rand_pol[5] = (random_data[5] & RAND_UP_BOUND) - ((random_data[5] >> 3) & RAND_UP_BOUND);
}

//~ Computes pa(X)*pb(X) mod(E)
void mult_mod_poly(int64_t *rop, int64_t *pa, int64_t *pb, int64_t *rand_zero){

	int64_t rand_op[NB_COEFF];
	int128 tmp_prod_result[NB_COEFF];

	add_poly(rand_op, pb, rand_zero); // computes: B' = B + J

	//~ computation of: A.B'[E]
	tmp_prod_result[0] = (int128)pa[0] * rand_op[0] - (((int128)pa[1] * rand_op[5] + (int128)pa[2] * rand_op[4] + (int128)pa[3] * rand_op[3] + (int128)pa[4] * rand_op[2] + (int128)pa[5] * rand_op[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * rand_op[1] + (int128)pa[1] * rand_op[0] - (((int128)pa[2] * rand_op[5] + (int128)pa[3] * rand_op[4] + (int128)pa[4] * rand_op[3] + (int128)pa[5] * rand_op[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * rand_op[2] + (int128)pa[1] * rand_op[1] + (int128)pa[2] * rand_op[0] - (((int128)pa[3] * rand_op[5] + (int128)pa[4] * rand_op[4] + (int128)pa[5] * rand_op[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * rand_op[3] + (int128)pa[1] * rand_op[2] + (int128)pa[2] * rand_op[1] + (int128)pa[3] * rand_op[0] - (((int128)pa[4] * rand_op[5] + (int128)pa[5] * rand_op[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * rand_op[4] + (int128)pa[1] * rand_op[3] + (int128)pa[2] * rand_op[2] + (int128)pa[3] * rand_op[1] + (int128)pa[4] * rand_op[0] - (((int128)pa[5] * rand_op[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * rand_op[5] + (int128)pa[1] * rand_op[4] + (int128)pa[2] * rand_op[3] + (int128)pa[3] * rand_op[2] + (int128)pa[4] * rand_op[1] + (int128)pa[5] * rand_op[0];

	internal_reduction(rop, tmp_prod_result);
	
	double_add_poly(rop, rop, rand_zero); //~last step: + 2.J
}

//~ Computes pa(X)^2 mod(E)
void square_mod_poly(int64_t *rop, int64_t *pa, int64_t *rand_zero){

	int64_t rand_op[NB_COEFF];
	int128 tmp_prod_result[NB_COEFF];

	add_poly(rand_op, pa, rand_zero); // computes: A' = A + J

	//~ computation of: A'^2[E]
	tmp_prod_result[0] = (int128)rand_op[0] * rand_op[0] - (((((int128)rand_op[4] * rand_op[2] + (int128)rand_op[5] * rand_op[1]) << 1) + (int128)rand_op[3] * rand_op[3]) << 2);
	tmp_prod_result[1] = (((int128)rand_op[1] * rand_op[0]) << 1) - (((int128)rand_op[4] * rand_op[3] + (int128)rand_op[5] * rand_op[2]) << 3);
	tmp_prod_result[2] = (((int128)rand_op[2] * rand_op[0]) << 1) + (int128)rand_op[1] * rand_op[1] - (((((int128)rand_op[5] * rand_op[3]) << 1) + (int128)rand_op[4] * rand_op[4]) << 2);
	tmp_prod_result[3] = (((int128)rand_op[2] * rand_op[1] + (int128)rand_op[3] * rand_op[0]) << 1) - (((int128)rand_op[5] * rand_op[4]) << 3);
	tmp_prod_result[4] = (((int128)rand_op[3] * rand_op[1] + (int128)rand_op[4] * rand_op[0]) << 1) + (int128)rand_op[2] * rand_op[2] - (((int128)rand_op[5] * rand_op[5]) << 2);
	tmp_prod_result[5] = (((int128)rand_op[3] * rand_op[2] + (int128)rand_op[4] * rand_op[1] + (int128)rand_op[5] * rand_op[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
	
	double_add_poly(rop, rop, rand_zero); //~last step: + 2.J
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmpQ_0, tmpQ_1, tmpQ_2, tmpQ_3, tmpQ_4, tmpQ_5;
	int128 tmpZero_0, tmpZero_1, tmpZero_2, tmpZero_3, tmpZero_4, tmpZero_5;

	//~ computation of : op*neg_inv_ri_rep_coeff mod(E, mont_phi)
	tmpQ_0 = (uint64_t)op[0] * 4801487534992572291UL + (uint64_t)op[3] * 1583920016881280228UL + (uint64_t)op[4] * 6798934983241754948UL + (uint64_t)op[5] * 8867761565262434256UL + (uint64_t)op[1] * 14814638768499915992UL + (uint64_t)op[2] * 17861618070575694308UL;
	tmpQ_1 = (uint64_t)op[1] * 4801487534992572291UL + (uint64_t)op[3] * 17861618070575694308UL + (uint64_t)op[5] * 6798934983241754948UL + (uint64_t)op[0] * 2394745627111779340UL + (uint64_t)op[2] * 14814638768499915992UL + (uint64_t)op[4] * 1583920016881280228UL;
	tmpQ_2 = (uint64_t)op[2] * 4801487534992572291UL + (uint64_t)op[4] * 17861618070575694308UL + (uint64_t)op[1] * 2394745627111779340UL + (uint64_t)op[0] * 7523638291044337071UL + (uint64_t)op[3] * 14814638768499915992UL + (uint64_t)op[5] * 1583920016881280228UL;
	tmpQ_3 = (uint64_t)op[3] * 4801487534992572291UL + (uint64_t)op[5] * 17861618070575694308UL + (uint64_t)op[0] * 8827392032634455751UL + (uint64_t)op[2] * 2394745627111779340UL + (uint64_t)op[1] * 7523638291044337071UL + (uint64_t)op[4] * 14814638768499915992UL;
	tmpQ_4 = (uint64_t)op[4] * 4801487534992572291UL + (uint64_t)op[1] * 8827392032634455751UL + (uint64_t)op[0] * 4757967519210852231UL + (uint64_t)op[3] * 2394745627111779340UL + (uint64_t)op[2] * 7523638291044337071UL + (uint64_t)op[5] * 14814638768499915992UL;
	tmpQ_5 = (uint64_t)op[5] * 4801487534992572291UL + (uint64_t)op[2] * 8827392032634455751UL + (uint64_t)op[0] * 10131398363157184714UL + (uint64_t)op[1] * 4757967519210852231UL + (uint64_t)op[4] * 2394745627111779340UL + (uint64_t)op[3] * 7523638291044337071UL;

	//~ computation of : tmp_q*red_int_coeff mod(E)
	tmpZero_0 = (int128)tmpQ_1 * -14477579451936L + (int128)tmpQ_2 * -3417025446624L + (int128)tmpQ_0 * -958314944803L + (int128)tmpQ_4 * 15312819876900L + (int128)tmpQ_5 * -4383068741696L + (int128)tmpQ_3 * 14245031534996L;
	tmpZero_1 = (int128)tmpQ_2 * -14477579451936L + (int128)tmpQ_3 * -3417025446624L + (int128)tmpQ_1 * -958314944803L + (int128)tmpQ_5 * 15312819876900L + (int128)tmpQ_0 * 1095767185424L + (int128)tmpQ_4 * 14245031534996L;
	tmpZero_2 = (int128)tmpQ_4 * -3417025446624L + (int128)tmpQ_3 * -14477579451936L + (int128)tmpQ_2 * -958314944803L + (int128)tmpQ_0 * -3828204969225L + (int128)tmpQ_1 * 1095767185424L + (int128)tmpQ_5 * 14245031534996L;
	tmpZero_3 = (int128)tmpQ_5 * -3417025446624L + (int128)tmpQ_4 * -14477579451936L + (int128)tmpQ_3 * -958314944803L + (int128)tmpQ_0 * -3561257883749L + (int128)tmpQ_1 * -3828204969225L + (int128)tmpQ_2 * 1095767185424L;
	tmpZero_4 = (int128)tmpQ_5 * -14477579451936L + (int128)tmpQ_4 * -958314944803L + (int128)tmpQ_1 * -3561257883749L + (int128)tmpQ_2 * -3828204969225L + (int128)tmpQ_3 * 1095767185424L + (int128)tmpQ_0 * 854256361656L;
	tmpZero_5 = (int128)tmpQ_5 * -958314944803L + (int128)tmpQ_2 * -3561257883749L + (int128)tmpQ_0 * 3619394862984L + (int128)tmpQ_3 * -3828204969225L + (int128)tmpQ_4 * 1095767185424L + (int128)tmpQ_1 * 854256361656L;

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmpZero_0) >> WORD_SIZE;
	rop[1] = (op[1] + tmpZero_1) >> WORD_SIZE;
	rop[2] = (op[2] + tmpZero_2) >> WORD_SIZE;
	rop[3] = (op[3] + tmpZero_3) >> WORD_SIZE;
	rop[4] = (op[4] + tmpZero_4) >> WORD_SIZE;
	rop[5] = (op[5] + tmpZero_5) >> WORD_SIZE;
}

void exact_coeffs_reduction(int64_t *rop, int64_t *op){
	int i;
	int128 tmp[NB_COEFF];

	for(i=0; i<NB_COEFF; i++)
		tmp[i] = (int128) op[i];

	internal_reduction(rop, tmp);
	
	tmp[0] = (int128)rop[0] * poly_P0[0] - (((int128)rop[1] * poly_P0[5] + (int128)rop[2] * poly_P0[4] + (int128)rop[3] * poly_P0[3] + (int128)rop[4] * poly_P0[2] + (int128)rop[5] * poly_P0[1]) << 2);
	tmp[1] = (int128)rop[0] * poly_P0[1] + (int128)rop[1] * poly_P0[0] - (((int128)rop[2] * poly_P0[5] + (int128)rop[3] * poly_P0[4] + (int128)rop[4] * poly_P0[3] + (int128)rop[5] * poly_P0[2]) << 2);
	tmp[2] = (int128)rop[0] * poly_P0[2] + (int128)rop[1] * poly_P0[1] + (int128)rop[2] * poly_P0[0] - (((int128)rop[3] * poly_P0[5] + (int128)rop[4] * poly_P0[4] + (int128)rop[5] * poly_P0[3]) << 2);
	tmp[3] = (int128)rop[0] * poly_P0[3] + (int128)rop[1] * poly_P0[2] + (int128)rop[2] * poly_P0[1] + (int128)rop[3] * poly_P0[0] - (((int128)rop[4] * poly_P0[5] + (int128)rop[5] * poly_P0[4]) << 2);
	tmp[4] = (int128)rop[0] * poly_P0[4] + (int128)rop[1] * poly_P0[3] + (int128)rop[2] * poly_P0[2] + (int128)rop[3] * poly_P0[1] + (int128)rop[4] * poly_P0[0] - (((int128)rop[5] * poly_P0[5]) << 2);
	tmp[5] = (int128)rop[0] * poly_P0[5] + (int128)rop[1] * poly_P0[4] + (int128)rop[2] * poly_P0[3] + (int128)rop[3] * poly_P0[2] + (int128)rop[4] * poly_P0[1] + (int128)rop[5] * poly_P0[0];
	
	internal_reduction(rop, tmp);
}


















