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


//~ Compilation command : gcc -Wall -O3 main.c -o main -lgmp 
//~ Execution command : ./main

//~ Important : polynomials representations form is P(X) = a0 + ... + an.X^n = (a0, ..., an).


int main(void){
	
	srand(time(NULL));
	
	unsigned long seed = time(NULL);
	gmp_randstate_t r;
	gmp_randinit_default(r);
	gmp_randseed_ui(r, seed);
	
	mpz_t A[NTEST], B[NTEST];
	for (int i=0; i<NTEST; i++){
		mpz_init (A[i]);
		mpz_init (B[i]);
	}
	
	int64_t pa[NTEST][NB_COEFF];
	int64_t pb[NTEST][NB_COEFF];
	int64_t pc[NTEST][NB_COEFF];


	init_data();
	
	//~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	uint64_t cycles[NTEST*NSAMPLES]={0};

	unsigned long long timermin , timermax, meanTimer1min =0, meanTimer1max =0, t1,t2, diff_t;
	unsigned long long *statTimer1 ;
	
	
	// cache memory heating
	for(int i=0;i<NTEST;i++){
		mpz_urandomm(A[i], r, modul_p);
		mpz_urandomm(B[i], r, modul_p);
		
		from_int_to_amns(pa[i], A[i]);
		from_int_to_amns(pb[i], B[i]);
	}
	for(int i=0;i<NTEST;i++)
		mult_mod_poly(pc[i], pa[i], pb[i]);
	
	//~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// timing
	for(int i=0;i<NSAMPLES;i++){
		
		// Génération d'un jeu de paramètres aléatoires
		for(int j=0;j<NTEST;j++){
		mpz_urandomm(A[j], r, modul_p);
		mpz_urandomm(B[j], r, modul_p);
		
		from_int_to_amns(pa[j], A[j]);
		from_int_to_amns(pb[j], B[j]);
		}
		
		timermin = (unsigned long long int)0x1<<63;
		timermax = 0;
		for(int j=0;j<NTEST;j++){
			t1 = cpucyclesStart();
			
			// Appel de la fonction à mesurer, avec le jeu de paramètre précédant
			mult_mod_poly(pc[j], pa[j], pb[j]);
			
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

	printf("\nmult_mod_poly: min : %lld, max : %lld,  Q1 : %lld, Q2 : %lld, Q3 : %lld CPU cycles\n", meanTimer1min/NSAMPLES, meanTimer1max/NSAMPLES, statTimer1[0], statTimer1[1], statTimer1[2]);
	

	for (int i=0; i<NTEST; i++){
		mpz_clear (A[i]);
		mpz_clear (B[i]);
	}
	gmp_randclear(r);

	free_data();
	
	return 0;
}

