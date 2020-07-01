#include "amns_init.h"


void init_data(){
	int i;
	for(i=0; i<POLY_DEG; i++)
		mpz_init (gama_pow[i]);

	mpz_init (modul_p);


	mpz_set_str (modul_p, "115792089237316195423570985008687907853269984665640564039457584007913129638397", 10);

	mpz_set_str (gama_pow[0], "62725036072198853858007603808615797342193930046000307761413395459424990168884", 10);
	for(i=1; i<POLY_DEG; i++){
		mpz_mul (gama_pow[i], gama_pow[i-1], gama_pow[0]);
		mpz_mod (gama_pow[i], gama_pow[i], modul_p);
	}

	//~ IMPORTANT : initialisations above must be done before those below.
	compute_Pi_polys();
}


//~ computes representatives of of polynomials Pi
void compute_Pi_polys(){
	int i, l;
	int64_t tmp_poly[NB_COEFF];
	
	//~ computation of a representative of 'phi*rho'
	from_mont_domain(tmp_poly, poly_P1);

	l = NB_COEFF - 2;
	if (l > 0){
		mult_mod_poly(polys_P[0], poly_P1, tmp_poly);
		
		for(i=1; i<l; i++)
			mult_mod_poly(polys_P[i], polys_P[i-1], tmp_poly);
	}
}


void free_data(){
	int i;
	for(i=0; i<POLY_DEG; i++)
		mpz_clear (gama_pow[i]);

	mpz_clear (modul_p);
}

