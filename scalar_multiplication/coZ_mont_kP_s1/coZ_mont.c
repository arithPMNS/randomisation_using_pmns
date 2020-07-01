#include "coZ_mont.h" 

void init_curve_and_base_point(){
	mpz_inits (order_of_P, crv_a, crv_b, const_ONE, NULL);
	init_affPoint(&aff_P);
	
	//~ curve equation is : y^2 = x^3 + 5 = x^3 + ax + b
	mpz_set_ui(crv_b, 5);
	from_int_to_amns(curve_a, crv_a); // note : crv_a is set to 0 with the initialisation above; note: conv to mont_domain is done here 
	from_int_to_amns(curve_b, crv_b);  // note: conv to mont_domain is done here 
	
	mpz_set_ui(const_ONE, 1);
	
	mpz_set_str (aff_P.x, "33287780256326021305866951798235204535566731496774587847840493300116154562335", 10);
	mpz_set_str (aff_P.y, "57200006855776545486303128794809149920862104667428222900231750752046655528108", 10);
	mpz_set_str (order_of_P, "115792089237316195423570985008687907852920869663551026044856920202296625078603", 10);
}

void free_point_data(){
	mpz_clears (order_of_P, crv_a, crv_b, const_ONE, NULL);
	clear_affPoint(&aff_P);
}


void init_affPoint(affPoint *ap){
	mpz_init (ap->x);
	mpz_init (ap->y);
}

void clear_affPoint(affPoint *ap){
	mpz_clear (ap->x);
	mpz_clear (ap->y);
}


void print_point_in_jac(JacPoint *jp){
	printf("X : "); print_element(jp->X); printf("\n");
	printf("Y : "); print_element(jp->Y); printf("\n");
	printf("Z : "); print_element(jp->Z); printf("\n");
}


void print_point_in_aff(affPoint *ap){
	gmp_printf("x : %Zd\n", ap->x);
	gmp_printf("y : %Zd\n", ap->y);
}

//~ note: jp->Z is assume to be one
void from_aff_to_jacOneZ_toMontDomain(Inter_Point *jp, affPoint *ap){
	from_int_to_amns(jp->X, ap->x);
	from_int_to_amns(jp->Y, ap->y);
}


void from_aff_to_jac_toMontDomain(JacPoint *jp, affPoint *ap){
	from_int_to_amns(jp->X, ap->x);
	from_int_to_amns(jp->Y, ap->y);
	from_int_to_amns(jp->Z, const_ONE); // computation of a (random) representation of 'phi'
}

//~ IMPORTANT : assumes elements of 'ap' already initialised
void from_jac_to_aff_outOfMontDomain(affPoint *ap, JacPoint *jp){
	
	mpz_t tz1, tz2, jpx, jpy, jpz;
	mpz_inits (tz1, tz2, jpx, jpy, jpz, NULL);
	
	from_amns_to_int(jpx, jp->X);
	from_amns_to_int(jpy, jp->Y);
	from_amns_to_int(jpz, jp->Z);
	
	mpz_invert (tz1, jpz, modul_p); // computation of 1/Z
	
	mpz_mul (tz2, tz1, tz1);
	mpz_mod (tz2, tz2, modul_p); // computation of 1/Z^2
	
	mpz_mul (ap->x, jpx, tz2);
	mpz_mod (ap->x, ap->x, modul_p);
	
	mpz_mul (ap->y, jpy, tz2);
	mpz_mod (ap->y, ap->y, modul_p);
	mpz_mul (ap->y, ap->y, tz1);
	mpz_mod (ap->y, ap->y, modul_p);
	
	mpz_clears (tz1, tz2, jpx, jpy, jpz, NULL);
}


int is_on_curve_aff(affPoint *ap){
	int rep;
	mpz_t l_op, r_op, tmp1, tmp2;
	mpz_inits (l_op, r_op, tmp1, tmp2, NULL);
	
	
	mpz_mul (l_op, ap->y, ap->y);
	mpz_mod (l_op, l_op, modul_p);
	
	
	mpz_set(r_op, crv_b);
	
	mpz_mul (tmp1, crv_a, ap->x);
	mpz_mod (tmp1, tmp1, modul_p);
	
	mpz_mul (tmp2, ap->x, ap->x);
	mpz_mod (tmp2, tmp2, modul_p);
	mpz_mul (tmp2, tmp2, ap->x);
	mpz_mod (tmp2, tmp2, modul_p);
	
	mpz_add(r_op, r_op, tmp1);
	mpz_add(r_op, r_op, tmp2);
	mpz_mod (r_op, r_op, modul_p);
	
	
	rep = mpz_cmp(l_op, r_op);
	
	mpz_clears (l_op, r_op, tmp1, tmp2, NULL);
	
	return (rep == 0);
}


int is_on_curve_jac(JacPoint *jp){
	
	int rep;
	affPoint tmp;
	init_affPoint(&tmp);
	
	from_jac_to_aff_outOfMontDomain(&tmp, jp);
	
	rep = is_on_curve_aff(&tmp);
	
	clear_affPoint(&tmp);
	
	return rep;
}


//~ ------------------------------------------------------------

//~ IMPORTANT : two exact_coeffs_reduction are done here to reduce the required number of 'free additions' from 12 to 3
//~ NOTE : it is assumed that 'P->Z' is 1
//~ Action: (R, S) <- (2*P, P), with the same Z
void dblu_with_free_adds(int64_t *Z_out, Inter_Point *R, Inter_Point *S, Inter_Point *P){
	
	int64_t B[NB_COEFF];
	int64_t E[NB_COEFF];
	int64_t L[NB_COEFF]; 
	
	square_mod_poly(B, P->X);
	square_mod_poly(E, P->Y);
	square_mod_poly(L, E);
	
	lshift_poly_coeffs(S->Y, L, 3);
	exact_coeffs_reduction(S->Y, S->Y); // first
	
	add_poly(E, E, P->X); 
	square_mod_poly(E, E);
	sub_poly(E, E, B); 
	sub_poly(E, E, L); 
	double_poly_coeffs(S->X, E);
	exact_coeffs_reduction(S->X, S->X); // second
	
	scalar_mult_poly(B, B, 3);
	add_poly(B, B, curve_a); 
	
	double_poly_coeffs(E, S->X);
	square_mod_poly(R->X, B);
	sub_poly(R->X, R->X, E); 
	
	sub_poly(L, S->X, R->X); 
	mult_mod_poly(R->Y, L, B);
	sub_poly(R->Y, R->Y, S->Y); 
	
	double_poly_coeffs(Z_out, P->Y);
}

//~ Action: (R, S) = (P+Q, P) <- (P, Q), with the same R.Z = S.Z = Z_out, where P.Z = Q.Z = Z_in
//~ Note: the sequence of instructions allows inputs and outputs to be the same.
void zaddu_with_free_adds(Inter_Point *R, Inter_Point *S, Inter_Point *P, Inter_Point *Q, int64_t *Z_in, int64_t *Z_out){
	
	int64_t C[NB_COEFF];
	int64_t D[NB_COEFF];
	int64_t W[NB_COEFF];
	
	sub_poly(C, P->X, Q->X); 
	
	mult_mod_poly(Z_out, Z_in, C);
	
	square_mod_poly(C, C);
	mult_mod_poly(S->X, P->X, C); // W1
	mult_mod_poly(W, Q->X, C);    // W2
	
	sub_poly(C, P->Y, Q->Y);
	square_mod_poly(D, C);
	
	sub_poly(R->X, D, S->X);
	sub_poly(R->X, R->X, W);
	
	sub_poly(W, S->X, W);
	mult_mod_poly(S->Y, P->Y, W);
	
	sub_poly(D, S->X, R->X);
	mult_mod_poly(R->Y, D, C);
	sub_poly(R->Y, R->Y, S->Y);
}


//~ Action: (R, S) = (P+Q, P-Q) <- (P, Q), with the same R.Z = S.Z = Z_out, where P.Z = Q.Z = Z_in
//~ Note: the sequence of instructions allows inputs and outputs to be the same.
void zaddc_with_free_adds(Inter_Point *R, Inter_Point *S, Inter_Point *P, Inter_Point *Q, int64_t *Z_in, int64_t *Z_out){
	
	int64_t C[NB_COEFF];
	int64_t D1[NB_COEFF];
	int64_t D2[NB_COEFF];
	int64_t W1[NB_COEFF];
	int64_t W2[NB_COEFF];
	
	sub_poly(C, P->X, Q->X); 
	
	mult_mod_poly(Z_out, Z_in, C);
	
	square_mod_poly(C, C);
	mult_mod_poly(W1, P->X, C);
	mult_mod_poly(W2, Q->X, C);
	
	sub_poly(C, P->Y, Q->Y);
	square_mod_poly(D1, C);
	
	sub_poly(R->X, D1, W1);
	sub_poly(R->X, R->X, W2); // X3
	
	sub_poly(D1, W1, W2);
	mult_mod_poly(D1, P->Y, D1); // A1
	
	
	add_poly(D2, P->Y, Q->Y);  // Y1 + Y2
	
	sub_poly(R->Y, W1, R->X);
	mult_mod_poly(R->Y, R->Y, C);
	sub_poly(R->Y, R->Y, D1);     // Y3
	
	square_mod_poly(S->X, D2); 
	sub_poly(S->X, S->X, W1);
	sub_poly(S->X, S->X, W2);   // X4
	
	sub_poly(S->Y, W1, S->X);
	mult_mod_poly(S->Y, S->Y, D2);
	sub_poly(S->Y, S->Y, D1);     // Y3
} 



//~ Here, it is assumed that 'op->Z' is 1, i.e. affines coordinates
void coZ_montgomery_scal_mult_with_free_adds(JacPoint *rop_Q, Inter_Point *op_P, char *scal_bits){
	
	int i,b;
	CoZ_Points tmp_coZ_points;
	
	
	dblu_with_free_adds(tmp_coZ_points.Z, &(tmp_coZ_points.R[1]), &(tmp_coZ_points.R[0]), op_P);
	
	for(i=(SCAL_SIZE-2); i>=0; i--){
		b = scal_bits[i];
		
		zaddc_with_free_adds(&(tmp_coZ_points.R[1-b]), &(tmp_coZ_points.R[b]), &(tmp_coZ_points.R[b]), &(tmp_coZ_points.R[1-b]), tmp_coZ_points.Z, tmp_coZ_points.Z);
		
		zaddu_with_free_adds(&(tmp_coZ_points.R[b]), &(tmp_coZ_points.R[1-b]), &(tmp_coZ_points.R[1-b]), &(tmp_coZ_points.R[b]), tmp_coZ_points.Z, tmp_coZ_points.Z);
	}
	
	copy_poly(rop_Q->X, tmp_coZ_points.R[0].X);
	copy_poly(rop_Q->Y, tmp_coZ_points.R[0].Y);
	copy_poly(rop_Q->Z, tmp_coZ_points.Z);
	//~ copy_poly(rop_Q->X, tmp_coZ_points.R[1].X);
	//~ copy_poly(rop_Q->Y, tmp_coZ_points.R[1].Y);
	//~ copy_poly(rop_Q->Z, tmp_coZ_points.Z);
}



//~ ------------------------------------------------------------

//~ IMPORTANT : two exact_coeffs_reduction are done here to reduce the required number of 'free additions' from 12 to 3
//~ NOTE : it is assumed that 'P->Z' is 1
//~ Action: (R, S) <- (2*P, P), with the same Z
void dblu_without_free_adds(int64_t *Z_out, Inter_Point *R, Inter_Point *S, Inter_Point *P){
	
	int64_t B[NB_COEFF];
	int64_t E[NB_COEFF];
	int64_t L[NB_COEFF]; 
	
	square_mod_poly(B, P->X);
	square_mod_poly(E, P->Y);
	square_mod_poly(L, E);
	
	lshift_poly_coeffs(S->Y, L, 3);
	exact_coeffs_reduction(S->Y, S->Y);  // <-------
	
	add_poly(E, E, P->X);
	exact_coeffs_reduction(E, E);  // <-------
	square_mod_poly(E, E);
	sub_poly(E, E, B); 
	sub_poly(E, E, L); 
	double_poly_coeffs(S->X, E);
	exact_coeffs_reduction(S->X, S->X);  // <-------
	
	scalar_mult_poly(B, B, 3);
	add_poly(B, B, curve_a);
	exact_coeffs_reduction(B, B);  // <-------
	
	double_poly_coeffs(E, S->X);
	exact_coeffs_reduction(E, E);  // <-------
	square_mod_poly(R->X, B);
	sub_poly(R->X, R->X, E);
	exact_coeffs_reduction(R->X, R->X);  // <-------
	
	sub_poly(L, S->X, R->X);
	exact_coeffs_reduction(L, L);  // <-------
	mult_mod_poly(R->Y, L, B);
	sub_poly(R->Y, R->Y, S->Y);
	exact_coeffs_reduction(R->Y, R->Y);  // <-------
	
	double_poly_coeffs(Z_out, P->Y);
	exact_coeffs_reduction(Z_out, Z_out);  // <-------
}


//~ Action: (R, S) = (P+Q, P) <- (P, Q), with the same R.Z = S.Z = Z_out, where P.Z = Q.Z = Z_in
//~ Note: the sequence of instructions allows inputs and outputs to be the same.
void zaddu_without_free_adds(Inter_Point *R, Inter_Point *S, Inter_Point *P, Inter_Point *Q, int64_t *Z_in, int64_t *Z_out){
	
	int64_t C[NB_COEFF];
	int64_t D[NB_COEFF];
	int64_t W[NB_COEFF];
	
	sub_poly(C, P->X, Q->X);
	exact_coeffs_reduction(C, C);  // <-------
	
	mult_mod_poly(Z_out, Z_in, C);
	
	square_mod_poly(C, C);
	mult_mod_poly(S->X, P->X, C); // W1
	mult_mod_poly(W, Q->X, C);    // W2
	
	sub_poly(C, P->Y, Q->Y);
	exact_coeffs_reduction(C, C);  // <-------
	square_mod_poly(D, C);
	
	sub_poly(R->X, D, S->X);
	sub_poly(R->X, R->X, W);
	exact_coeffs_reduction(R->X, R->X);  // <-------
	
	sub_poly(W, S->X, W);
	exact_coeffs_reduction(W, W);  // <-------
	mult_mod_poly(S->Y, P->Y, W);
	
	sub_poly(D, S->X, R->X);
	exact_coeffs_reduction(D, D);  // <-------
	mult_mod_poly(R->Y, D, C);
	sub_poly(R->Y, R->Y, S->Y);
	exact_coeffs_reduction(R->Y, R->Y);  // <-------
}

//~ Action: (R, S) = (P+Q, P-Q) <- (P, Q), with the same R.Z = S.Z = Z_out, where P.Z = Q.Z = Z_in
//~ Note: the sequence of instructions allows inputs and outputs to be the same.
void zaddc_without_free_adds(Inter_Point *R, Inter_Point *S, Inter_Point *P, Inter_Point *Q, int64_t *Z_in, int64_t *Z_out){
	
	int64_t C[NB_COEFF];
	int64_t D1[NB_COEFF];
	int64_t D2[NB_COEFF];
	int64_t W1[NB_COEFF];
	int64_t W2[NB_COEFF];
	
	sub_poly(C, P->X, Q->X);
	exact_coeffs_reduction(C, C);  			// <-------
	
	mult_mod_poly(Z_out, Z_in, C);
	
	square_mod_poly(C, C);
	mult_mod_poly(W1, P->X, C);
	mult_mod_poly(W2, Q->X, C);
	
	sub_poly(C, P->Y, Q->Y);
	exact_coeffs_reduction(C, C);  			// <-------
	square_mod_poly(D1, C);
	
	sub_poly(R->X, D1, W1);
	sub_poly(R->X, R->X, W2); // X3
	exact_coeffs_reduction(R->X, R->X);		// <-------
	
	sub_poly(D1, W1, W2);
	exact_coeffs_reduction(D1, D1);			// <-------
	mult_mod_poly(D1, P->Y, D1); // A1
	
	
	add_poly(D2, P->Y, Q->Y);  // Y1 + Y2
	exact_coeffs_reduction(D2, D2);			// <-------
	
	sub_poly(R->Y, W1, R->X);
	exact_coeffs_reduction(R->Y, R->Y);			// <-------
	mult_mod_poly(R->Y, R->Y, C);
	sub_poly(R->Y, R->Y, D1);     // Y3
	exact_coeffs_reduction(R->Y, R->Y);			// <-------
	
	square_mod_poly(S->X, D2); 
	sub_poly(S->X, S->X, W1);
	sub_poly(S->X, S->X, W2);   // X4
	exact_coeffs_reduction(S->X, S->X);			// <-------
	
	sub_poly(S->Y, W1, S->X);
	exact_coeffs_reduction(S->Y, S->Y);			// <-------
	mult_mod_poly(S->Y, S->Y, D2);
	sub_poly(S->Y, S->Y, D1);     // Y3
	exact_coeffs_reduction(S->Y, S->Y);			// <-------
} 

//~ Here, it is assumed that 'op->Z' is 1, i.e. affines coordinates
void coZ_montgomery_scal_mult_without_free_adds(JacPoint *rop_Q, Inter_Point *op_P, char *scal_bits){

	int i,b;
	CoZ_Points tmp_coZ_points;
	
	
	dblu_without_free_adds(tmp_coZ_points.Z, &(tmp_coZ_points.R[1]), &(tmp_coZ_points.R[0]), op_P);
	
	for(i=(SCAL_SIZE-2); i>=0; i--){
		b = scal_bits[i];
		
		zaddc_without_free_adds(&(tmp_coZ_points.R[1-b]), &(tmp_coZ_points.R[b]), &(tmp_coZ_points.R[b]), &(tmp_coZ_points.R[1-b]), tmp_coZ_points.Z, tmp_coZ_points.Z);
		
		zaddu_without_free_adds(&(tmp_coZ_points.R[b]), &(tmp_coZ_points.R[1-b]), &(tmp_coZ_points.R[1-b]), &(tmp_coZ_points.R[b]), tmp_coZ_points.Z, tmp_coZ_points.Z);
	}
	
	copy_poly(rop_Q->X, tmp_coZ_points.R[0].X);
	copy_poly(rop_Q->Y, tmp_coZ_points.R[0].Y);
	copy_poly(rop_Q->Z, tmp_coZ_points.Z);
	//~ copy_poly(rop_Q->X, tmp_coZ_points.R[1].X);
	//~ copy_poly(rop_Q->Y, tmp_coZ_points.R[1].Y);
	//~ copy_poly(rop_Q->Z, tmp_coZ_points.Z);
}



































