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

//~ computes : pa + 2.pb
void double_add_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] + 2*pb[j];
}

//~ computes : pa - 2.pb
void double_sub_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] - 2*pb[j];
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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(E)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmpQ_0, tmpQ_1, tmpQ_2, tmpQ_3, tmpQ_4, tmpQ_5;
	int128 tmpZero_0, tmpZero_1, tmpZero_2, tmpZero_3, tmpZero_4, tmpZero_5;

	//~ computation of : op*neg_inv_ri_rep_coeff mod(E, mont_phi)
	tmpQ_0 = (uint64_t)op[3] * 9066109803895389729UL + (uint64_t)op[4] * 17228072911253529543UL + (uint64_t)op[5] * 1878471730737693744UL + (uint64_t)op[2] * 11834913054877804533UL + (uint64_t)op[1] * 9815045943720945279UL + (uint64_t)op[0] * 1791998369291743285UL;
	tmpQ_1 = (uint64_t)op[4] * 9066109803895389729UL + (uint64_t)op[5] * 17228072911253529543UL + (uint64_t)op[0] * 626157243579231248UL + (uint64_t)op[3] * 11834913054877804533UL + (uint64_t)op[2] * 9815045943720945279UL + (uint64_t)op[1] * 1791998369291743285UL;
	tmpQ_2 = (uint64_t)op[5] * 9066109803895389729UL + (uint64_t)op[0] * 5742690970417843181UL + (uint64_t)op[1] * 626157243579231248UL + (uint64_t)op[4] * 11834913054877804533UL + (uint64_t)op[2] * 1791998369291743285UL + (uint64_t)op[3] * 9815045943720945279UL;
	tmpQ_3 = (uint64_t)op[0] * 3022036601298463243UL + (uint64_t)op[1] * 5742690970417843181UL + (uint64_t)op[2] * 626157243579231248UL + (uint64_t)op[5] * 11834913054877804533UL + (uint64_t)op[3] * 1791998369291743285UL + (uint64_t)op[4] * 9815045943720945279UL;
	tmpQ_4 = (uint64_t)op[0] * 3944971018292601511UL + (uint64_t)op[1] * 3022036601298463243UL + (uint64_t)op[2] * 5742690970417843181UL + (uint64_t)op[3] * 626157243579231248UL + (uint64_t)op[4] * 1791998369291743285UL + (uint64_t)op[5] * 9815045943720945279UL;
	tmpQ_5 = (uint64_t)op[1] * 3944971018292601511UL + (uint64_t)op[2] * 3022036601298463243UL + (uint64_t)op[3] * 5742690970417843181UL + (uint64_t)op[4] * 626157243579231248UL + (uint64_t)op[0] * 3271681981240315093UL + (uint64_t)op[5] * 1791998369291743285UL;

	//~ computation of : tmp_q*red_int_coeff mod(E)
	tmpZero_0 = (int128)tmpQ_3 * -5091835506183L + (int128)tmpQ_0 * 2587093525133L + (int128)tmpQ_5 * -2514491652669L + (int128)tmpQ_1 * 4697248739220L + (int128)tmpQ_2 * 369410497077L + (int128)tmpQ_4 * 8830871768445L;
	tmpZero_1 = (int128)tmpQ_4 * -5091835506183L + (int128)tmpQ_1 * 2587093525133L + (int128)tmpQ_2 * 4697248739220L + (int128)tmpQ_3 * 369410497077L + (int128)tmpQ_5 * 8830871768445L + (int128)tmpQ_0 * -838163884223L;
	tmpZero_2 = (int128)tmpQ_5 * -5091835506183L + (int128)tmpQ_2 * 2587093525133L + (int128)tmpQ_3 * 4697248739220L + (int128)tmpQ_4 * 369410497077L + (int128)tmpQ_0 * 2943623922815L + (int128)tmpQ_1 * -838163884223L;
	tmpZero_3 = (int128)tmpQ_3 * 2587093525133L + (int128)tmpQ_0 * -1697278502061L + (int128)tmpQ_4 * 4697248739220L + (int128)tmpQ_5 * 369410497077L + (int128)tmpQ_1 * 2943623922815L + (int128)tmpQ_2 * -838163884223L;
	tmpZero_4 = (int128)tmpQ_0 * 123136832359L + (int128)tmpQ_4 * 2587093525133L + (int128)tmpQ_1 * -1697278502061L + (int128)tmpQ_5 * 4697248739220L + (int128)tmpQ_3 * -838163884223L + (int128)tmpQ_2 * 2943623922815L;
	tmpZero_5 = (int128)tmpQ_1 * 123136832359L + (int128)tmpQ_5 * 2587093525133L + (int128)tmpQ_2 * -1697278502061L + (int128)tmpQ_3 * 2943623922815L + (int128)tmpQ_0 * 1565749579740L + (int128)tmpQ_4 * -838163884223L;

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

	tmp[0] = (int128)rop[0] * polys_P[0][0] + (((int128)rop[1] * polys_P[0][5] + (int128)rop[2] * polys_P[0][4] + (int128)rop[3] * polys_P[0][3] + (int128)rop[4] * polys_P[0][2] + (int128)rop[5] * polys_P[0][1]) * 3);
	tmp[1] = (int128)rop[0] * polys_P[0][1] + (int128)rop[1] * polys_P[0][0] + (((int128)rop[2] * polys_P[0][5] + (int128)rop[3] * polys_P[0][4] + (int128)rop[4] * polys_P[0][3] + (int128)rop[5] * polys_P[0][2]) * 3);
	tmp[2] = (int128)rop[0] * polys_P[0][2] + (int128)rop[1] * polys_P[0][1] + (int128)rop[2] * polys_P[0][0] + (((int128)rop[3] * polys_P[0][5] + (int128)rop[4] * polys_P[0][4] + (int128)rop[5] * polys_P[0][3]) * 3);
	tmp[3] = (int128)rop[0] * polys_P[0][3] + (int128)rop[1] * polys_P[0][2] + (int128)rop[2] * polys_P[0][1] + (int128)rop[3] * polys_P[0][0] + (((int128)rop[4] * polys_P[0][5] + (int128)rop[5] * polys_P[0][4]) * 3);
	tmp[4] = (int128)rop[0] * polys_P[0][4] + (int128)rop[1] * polys_P[0][3] + (int128)rop[2] * polys_P[0][2] + (int128)rop[3] * polys_P[0][1] + (int128)rop[4] * polys_P[0][0] + (((int128)rop[5] * polys_P[0][5]) * 3);
	tmp[5] = (int128)rop[0] * polys_P[0][5] + (int128)rop[1] * polys_P[0][4] + (int128)rop[2] * polys_P[0][3] + (int128)rop[3] * polys_P[0][2] + (int128)rop[4] * polys_P[0][1] + (int128)rop[5] * polys_P[0][0];

	internal_reduction(rop, tmp);
}

