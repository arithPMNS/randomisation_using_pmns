#include "add_mult_poly.h"


void add_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] + pb[j];
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

	rand_zero[0] = rand_pol[0] * -2243091601121 + rand_pol[4] * 3489445074728 + rand_pol[5] * 9255697449168 + rand_pol[3] * 14115396788464 + rand_pol[1] * -6855691867408 + rand_pol[2] * 16076437253240;
	rand_zero[1] = rand_pol[1] * -2243091601121 + rand_pol[5] * 3489445074728 + rand_pol[4] * 14115396788464 + rand_pol[2] * -6855691867408 + rand_pol[0] * -2313924362292 + rand_pol[3] * 16076437253240;
	rand_zero[2] = rand_pol[2] * -2243091601121 + rand_pol[0] * -872361268682 + rand_pol[5] * 14115396788464 + rand_pol[3] * -6855691867408 + rand_pol[1] * -2313924362292 + rand_pol[4] * 16076437253240;
	rand_zero[3] = rand_pol[3] * -2243091601121 + rand_pol[1] * -872361268682 + rand_pol[4] * -6855691867408 + rand_pol[2] * -2313924362292 + rand_pol[5] * 16076437253240 + rand_pol[0] * -3528849197116;
	rand_zero[4] = rand_pol[4] * -2243091601121 + rand_pol[2] * -872361268682 + rand_pol[5] * -6855691867408 + rand_pol[3] * -2313924362292 + rand_pol[1] * -3528849197116 + rand_pol[0] * -4019109313310;
	rand_zero[5] = rand_pol[5] * -2243091601121 + rand_pol[0] * 1713922966852 + rand_pol[3] * -872361268682 + rand_pol[4] * -2313924362292 + rand_pol[2] * -3528849197116 + rand_pol[1] * -4019109313310;
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
void mult_mod_poly(int64_t *rop, int64_t *pa, int64_t *pb){

	int128 tmp_prod_result[NB_COEFF];

	//~ computation of: A.B'[E]
	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(E)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	//~ computation of: A'^2[E]
	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmpQ_0, tmpQ_1, tmpQ_2, tmpQ_3, tmpQ_4, tmpQ_5;
	int128 tmpZero_0, tmpZero_1, tmpZero_2, tmpZero_3, tmpZero_4, tmpZero_5;

	//~ computation of : op*neg_inv_ri_rep_coeff mod(E, mont_phi)
	tmpQ_0 = (uint64_t)op[0] * 10132241992947132833UL + (uint64_t)op[4] * 5185227616001927592UL + (uint64_t)op[3] * 10594776064724140208UL + (uint64_t)op[1] * 12999443846141959472UL + (uint64_t)op[2] * 11738032480074340968UL + (uint64_t)op[5] * 1805543253268868560UL;
	tmpQ_1 = (uint64_t)op[1] * 10132241992947132833UL + (uint64_t)op[5] * 5185227616001927592UL + (uint64_t)op[4] * 10594776064724140208UL + (uint64_t)op[0] * 8771986223537558668UL + (uint64_t)op[2] * 12999443846141959472UL + (uint64_t)op[3] * 11738032480074340968UL;
	tmpQ_2 = (uint64_t)op[2] * 10132241992947132833UL + (uint64_t)op[4] * 11738032480074340968UL + (uint64_t)op[5] * 10594776064724140208UL + (uint64_t)op[1] * 8771986223537558668UL + (uint64_t)op[3] * 12999443846141959472UL + (uint64_t)op[0] * 17150437169709069718UL;
	tmpQ_3 = (uint64_t)op[3] * 10132241992947132833UL + (uint64_t)op[5] * 11738032480074340968UL + (uint64_t)op[2] * 8771986223537558668UL + (uint64_t)op[4] * 12999443846141959472UL + (uint64_t)op[0] * 6574678020673740756UL + (uint64_t)op[1] * 17150437169709069718UL;
	tmpQ_4 = (uint64_t)op[4] * 10132241992947132833UL + (uint64_t)op[0] * 10900549935263578470UL + (uint64_t)op[3] * 8771986223537558668UL + (uint64_t)op[5] * 12999443846141959472UL + (uint64_t)op[1] * 6574678020673740756UL + (uint64_t)op[2] * 17150437169709069718UL;
	tmpQ_5 = (uint64_t)op[5] * 10132241992947132833UL + (uint64_t)op[1] * 10900549935263578470UL + (uint64_t)op[4] * 8771986223537558668UL + (uint64_t)op[2] * 6574678020673740756UL + (uint64_t)op[3] * 17150437169709069718UL + (uint64_t)op[0] * 1361825056891898036UL;

	//~ computation of : tmp_q*red_int_coeff mod(E)
	tmpZero_0 = (int128)tmpQ_0 * -2243091601121L + (int128)tmpQ_4 * 3489445074728L + (int128)tmpQ_5 * 9255697449168L + (int128)tmpQ_3 * 14115396788464L + (int128)tmpQ_1 * -6855691867408L + (int128)tmpQ_2 * 16076437253240L;
	tmpZero_1 = (int128)tmpQ_1 * -2243091601121L + (int128)tmpQ_5 * 3489445074728L + (int128)tmpQ_4 * 14115396788464L + (int128)tmpQ_2 * -6855691867408L + (int128)tmpQ_0 * -2313924362292L + (int128)tmpQ_3 * 16076437253240L;
	tmpZero_2 = (int128)tmpQ_2 * -2243091601121L + (int128)tmpQ_0 * -872361268682L + (int128)tmpQ_5 * 14115396788464L + (int128)tmpQ_3 * -6855691867408L + (int128)tmpQ_1 * -2313924362292L + (int128)tmpQ_4 * 16076437253240L;
	tmpZero_3 = (int128)tmpQ_3 * -2243091601121L + (int128)tmpQ_1 * -872361268682L + (int128)tmpQ_4 * -6855691867408L + (int128)tmpQ_2 * -2313924362292L + (int128)tmpQ_5 * 16076437253240L + (int128)tmpQ_0 * -3528849197116L;
	tmpZero_4 = (int128)tmpQ_4 * -2243091601121L + (int128)tmpQ_2 * -872361268682L + (int128)tmpQ_5 * -6855691867408L + (int128)tmpQ_3 * -2313924362292L + (int128)tmpQ_1 * -3528849197116L + (int128)tmpQ_0 * -4019109313310L;
	tmpZero_5 = (int128)tmpQ_5 * -2243091601121L + (int128)tmpQ_0 * 1713922966852L + (int128)tmpQ_3 * -872361268682L + (int128)tmpQ_4 * -2313924362292L + (int128)tmpQ_2 * -3528849197116L + (int128)tmpQ_1 * -4019109313310L;

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


















