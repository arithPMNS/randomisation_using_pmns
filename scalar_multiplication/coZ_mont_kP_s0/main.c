#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <gmp.h>

#include "intel_measurement_stuff.c"

#include "structs_data.h"
#include "add_mult_poly.c"
#include "useful_functs.c"
#include "amns_init.c"
#include "coZ_mont.c"


//~ Compilation command : gcc -Wall -O3 main.c -o main -lgmp
//~ Execution command : ./main

//~ Important : polynomials representations form is P(X) = a0 + ... + an.X^n = (a0, ..., an).


int main(void){
	
	srand(time(NULL));
	
	unsigned long seed = time(NULL);
	gmp_randstate_t r;
	gmp_randinit_default(r);
	gmp_randseed_ui(r, seed);
	
	mpz_t scal_k[NTEST];
	for (int i=0; i<NTEST; i++)
		mpz_init (scal_k[i]);

	char scal_k_bits[NTEST][SCAL_SIZE];
	
	JacPoint jac_Q;
	affPoint aff_Q; // to represent the result in affine coordinates
	Inter_Point jac_P;
	
	
	init_data(); // for the AMNS
	init_curve_and_base_point(); // for the curve
	init_affPoint(&aff_Q);
	
	//~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	uint64_t cycles[NTEST*NSAMPLES]={0};

	unsigned long long timermin , timermax, meanTimer1min =0, meanTimer1max =0, t1,t2, diff_t;
	unsigned long long *statTimer1 ;
	
	
	//~ Génération des scalaires
	for(int i=0;i<NTEST;i++){
		mpz_urandomm(scal_k[i], r, order_of_P); // à faire après initialisation de la courbe
		while(mpz_tstbit(scal_k[i], (SCAL_SIZE-1)) == 0)
			mpz_urandomm(scal_k[i], r, order_of_P); // on veut un scalaire de taille toujours SCAL_SIZE bits, utile pour le Co-Z Échelle de Montgomrey (voir algo). IMPORTANT : s'assurer que 'order_of_P' le permet, sinon risque de boucle infinie ici ....
		
		//~ pour éviter de 'polluer' la multiplication scalaire en utilisant GNU MP
		for (int j = 0; j<SCAL_SIZE; j++)
			scal_k_bits[i][j] = mpz_tstbit(scal_k[i], j);
	}
	
	// cache memory heating
	for(int i=0;i<NTEST;i++){
		from_aff_to_jacOneZ_toMontDomain(&jac_P, &aff_P); // conversion step
		coZ_montgomery_scal_mult_with_free_adds(&jac_Q, &jac_P, scal_k_bits[i]);
	}
	
	//~ -------------------- ECSM with 'delta' big enough ----------------------------------------------
	
	// timing
	for(int i=0;i<NSAMPLES;i++){
		
		// Génération d'un jeu de paramètres aléatoires : génération des scalaires
		for(int j=0;j<NTEST;j++){
			mpz_urandomm(scal_k[j], r, order_of_P); 
			while(mpz_tstbit(scal_k[j], (SCAL_SIZE-1)) == 0)
				mpz_urandomm(scal_k[j], r, order_of_P); 
			for (int l = 0; l<SCAL_SIZE; l++) //~ pour éviter de 'polluer' la multiplication scalaire en utilisant GNU MP
				scal_k_bits[j][l] = mpz_tstbit(scal_k[j], l);
		}
		
		timermin = (unsigned long long int)0x1<<63;
		timermax = 0;
		for(int j=0;j<NTEST;j++){
			t1 = cpucyclesStart();
			
			// Appel des fonctions à mesurer, avec le jeu de paramètre précédant
			from_aff_to_jacOneZ_toMontDomain(&jac_P, &aff_P); // conversion step
			coZ_montgomery_scal_mult_with_free_adds(&jac_Q, &jac_P, scal_k_bits[j]);
			
			t2 = cpucyclesStop();
			if (t2 < t1){
				diff_t = 18446744073709551615ULL-t1;
				diff_t = t2+diff_t+1;
			}
			else
				diff_t = t2-t1;
			if(timermin>diff_t) timermin = diff_t;
			else if(timermax < diff_t) timermax = diff_t;
			cycles[i*NTEST+j]=diff_t;
		}

		meanTimer1min += timermin;
		meanTimer1max += timermax;
	}
	statTimer1   = quartiles(cycles,NTEST*NSAMPLES);

	

	//~ -------------------- results and timings ----------------------------------------------------
	
	printf("\n");
	print_point_in_aff(&aff_P);
	if (is_on_curve_aff(&aff_P))
		printf("\nbase point : OK\n\n");
	else
		printf("\nbase point : NOT OK\n\n");
	/*print_point_in_jac(&jac_P);
	if (is_on_curve_jac(&jac_P))
		printf("jac base point : OK\n\n");
	else
		printf("jac base point : NOT OK\n\n");*/

	printf("-------------------------- result --------------------------- \n\n");
	
	from_jac_to_aff_outOfMontDomain(&aff_Q, &jac_Q);
	print_point_in_jac(&jac_Q); printf("\n");
	print_point_in_aff(&aff_Q); printf("\n");
	if (is_on_curve_jac(&jac_Q))
		printf("result point : OK\n\n");
	else
		printf("result point : NOT OK\n\n");
	
	
	printf("\nscal_mult: min : %lld, max : %lld,  Q1 : %lld, Q2 : %lld, Q3 : %lld CPU cycles\n", meanTimer1min/NSAMPLES, meanTimer1max/NSAMPLES, statTimer1[0], statTimer1[1], statTimer1[2]);
	
	

	for (int i=0; i<NTEST; i++)
		mpz_clear (scal_k[i]);
	gmp_randclear(r);

	clear_affPoint(&aff_Q);
	free_point_data(); // for the curve
	free_data(); // for the AMNS
	
	return 0;
}


























