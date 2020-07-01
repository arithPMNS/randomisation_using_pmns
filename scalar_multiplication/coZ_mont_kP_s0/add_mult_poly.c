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

//~ Computes pa(X)*pb(X) mod(E)
void mult_mod_poly(int64_t *rop, int64_t *pa, int64_t *pb){

	int128 tmp_prod_result[NB_COEFF];

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
	tmpQ_0 = (uint64_t)op[3] * 98746733395996164UL + (uint64_t)op[4] * 2403407688631628396UL + (uint64_t)op[2] * 3259566405501951084UL + (uint64_t)op[0] * 12636426123386989239UL + (uint64_t)op[1] * 11077915543176631628UL + (uint64_t)op[5] * 6345162433761097372UL;
	tmpQ_1 = (uint64_t)op[4] * 98746733395996164UL + (uint64_t)op[2] * 11077915543176631628UL + (uint64_t)op[5] * 2403407688631628396UL + (uint64_t)op[3] * 3259566405501951084UL + (uint64_t)op[1] * 12636426123386989239UL + (uint64_t)op[0] * 12248767446841889369UL;
	tmpQ_2 = (uint64_t)op[5] * 98746733395996164UL + (uint64_t)op[0] * 4010834096269480805UL + (uint64_t)op[4] * 3259566405501951084UL + (uint64_t)op[3] * 11077915543176631628UL + (uint64_t)op[2] * 12636426123386989239UL + (uint64_t)op[1] * 12248767446841889369UL;
	tmpQ_3 = (uint64_t)op[1] * 4010834096269480805UL + (uint64_t)op[5] * 3259566405501951084UL + (uint64_t)op[4] * 11077915543176631628UL + (uint64_t)op[3] * 12636426123386989239UL + (uint64_t)op[2] * 12248767446841889369UL + (uint64_t)op[0] * 18422057390360552575UL;
	tmpQ_4 = (uint64_t)op[2] * 4010834096269480805UL + (uint64_t)op[5] * 11077915543176631628UL + (uint64_t)op[4] * 12636426123386989239UL + (uint64_t)op[3] * 12248767446841889369UL + (uint64_t)op[1] * 18422057390360552575UL + (uint64_t)op[0] * 13020166453906675941UL;
	tmpQ_5 = (uint64_t)op[1] * 13020166453906675941UL + (uint64_t)op[3] * 4010834096269480805UL + (uint64_t)op[0] * 15677265187915393709UL + (uint64_t)op[5] * 12636426123386989239UL + (uint64_t)op[4] * 12248767446841889369UL + (uint64_t)op[2] * 18422057390360552575UL;

	//~ computation of : tmp_q*red_int_coeff mod(E)
	tmpZero_0 = (int128)tmpQ_3 * -11878367806152L + (int128)tmpQ_0 * 4063330874669L + (int128)tmpQ_2 * -17513256816176L + (int128)tmpQ_4 * -8547885628664L + (int128)tmpQ_5 * -4216937025428L + (int128)tmpQ_1 * 16826431529368L;
	tmpZero_1 = (int128)tmpQ_0 * 1054234256357L + (int128)tmpQ_4 * -11878367806152L + (int128)tmpQ_1 * 4063330874669L + (int128)tmpQ_3 * -17513256816176L + (int128)tmpQ_5 * -8547885628664L + (int128)tmpQ_2 * 16826431529368L;
	tmpZero_2 = (int128)tmpQ_1 * 1054234256357L + (int128)tmpQ_5 * -11878367806152L + (int128)tmpQ_2 * 4063330874669L + (int128)tmpQ_4 * -17513256816176L + (int128)tmpQ_3 * 16826431529368L + (int128)tmpQ_0 * 2136971407166L;
	tmpZero_3 = (int128)tmpQ_2 * 1054234256357L + (int128)tmpQ_3 * 4063330874669L + (int128)tmpQ_5 * -17513256816176L + (int128)tmpQ_0 * 2969591951538L + (int128)tmpQ_4 * 16826431529368L + (int128)tmpQ_1 * 2136971407166L;
	tmpZero_4 = (int128)tmpQ_3 * 1054234256357L + (int128)tmpQ_0 * 4378314204044L + (int128)tmpQ_4 * 4063330874669L + (int128)tmpQ_1 * 2969591951538L + (int128)tmpQ_5 * 16826431529368L + (int128)tmpQ_2 * 2136971407166L;
	tmpZero_5 = (int128)tmpQ_4 * 1054234256357L + (int128)tmpQ_0 * -4206607882342L + (int128)tmpQ_1 * 4378314204044L + (int128)tmpQ_5 * 4063330874669L + (int128)tmpQ_2 * 2969591951538L + (int128)tmpQ_3 * 2136971407166L;

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


















