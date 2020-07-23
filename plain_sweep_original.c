/*
 *  Simple simulation of a sweep in space
 *
 *
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_randist.h>

// run by typing : gcc plain_sweep_original.c -o plain_sweep -lgsl
// ./plain_sweep L N s m tfinal ranseed



//Todo : when beneficial mutation is fixed in 10% of the demes, save n[i] and t. Get the right most deme with a beneficial allele, choose a random individual among them, and give a neutral mutation and track it.
// When the neutral mutation goes extinct, set n[i] and t to the saved values and start tracking again. Save n[i] and n_neutral[i]'s. 
// Later, there should be one more input value (argv) that sets how far out from the center of mass (i_cm = sum(n[i] * i) / N ) to put in a neutral mutation.

// Global variables
const gsl_rng_type * T;
gsl_rng * R;

void next_gen(unsigned int n[], double s, double mig, int L, unsigned int N) {

    double x[L]; // frequencies
    double xmig[L]; // after migration
    int i;

	// Initialize frequencies:
	for (i = 0; i < L; i++){
		x[i] = (double)n[i] / N;
	}
	// Migration:
	for (i = 1; i < L - 1; i++){
		xmig[i] = x[i] + mig * (0.5 * (x[i-1] + x[i+1]) - x[i]);
	}
	xmig[0] = x[0] + mig * 0.5 * (x[1] - x[0]);
	xmig[L-1] = x[L-1] + mig * 0.5 * (x[L-2] - x[L-1]);

    // Sampling and selection within demes
    for (i = 0; i < L; i++) {
		n[i] = gsl_ran_binomial(R, xmig[i] + s * xmig[i] * (1 - xmig[i]), N);
	}

}

/*function roundint to round a double to an integer*/
int roundint(double x) {
	int result;
	
	result = floor(x);
	if (x - result >= 0.5) result++;
	
	return result;
} /*roundint*/



int main(int argc, char *argv[]) {
    long SEED;
    double mig, s; // migration rate, selection
    unsigned int N, L, t, tfinal; // deme size, number of demes, time, max number of generations
    FILE *datafile, *paramfile;
    char *outfile = malloc(1000);
    char *outfile1 = malloc(1000);
    unsigned int i, j, ntot;
    //Initialize variables:
    if (argc != 7) {
        printf("Usage: L N s m tfinal ranseed\n");
        return 1; 
    }
    
    j = 1;
    L = (unsigned int) roundint(atof(argv[j++]));
    N = (unsigned int) roundint(atof(argv[j++]));
    s = atof(argv[j++]);
    mig = atof(argv[j++]);
    tfinal = (unsigned int) roundint(atof(argv[j++]));
    SEED = atof(argv[j++]);
    sprintf(outfile, "L=%u_N=%u_s=%f_m=%f_tfinal=%u", L, N, s, mig, tfinal);
	// outfile = argv[j++];
    outfile = strcat(outfile, ".txt");
	strcpy(outfile1,outfile);

    // Print out variables to parameter file:
    paramfile = fopen(strcat(outfile1, "_params.txt"), "w");
    fprintf(paramfile, "L = %u\nN = %u\ns = %f\nm = %f\ntfinal = %u\nseed = %ld\n", L, N, s, mig, tfinal, SEED);
    fclose(paramfile);
	
    // gsl random setup:
    gsl_rng_env_setup();
    T = gsl_rng_default;
    R = gsl_rng_alloc (T);
    gsl_rng_set(R,SEED);


    // Initialize population:
    unsigned int n[L];
    // Leftmost demes fixed for sweeping allele
    // Fill enough demes so that it is very unlikely to die out:
    if (1 + 10 / (N * s) > L / 2){
    	printf("Warning: meta-population is too small for selection to be strong.\n");
    	return 1;
    }

    for (i = 0; i < L; i++){
        n[i] = 0;
    }

    while (n[0] == 0){
        for (i = 0; i < L; i++){
            n[i] = 0;
        }
        n[0] = 1;
        //Open the datafile for recording:
        //datafile = fopen(strcat(outfile, ".txt"), "w");
        datafile = fopen(outfile, "w");
        // Run until alleles fixes (or goes extinct):
        for (t = 0; t < tfinal; t++) {

            // Record the status of the population and check for fixation:
            ntot = 0;
            for (i = 0; i < L; i++){
                fprintf(datafile, " %d", n[i]);
                ntot += n[i];
            }
            fprintf(datafile, "\n");

            // Stop if one of the alleles is fixed or extinct:
            if ((ntot == 0) || (ntot == N * L)) {break;}
            
            // Evolve for one generation
            next_gen(n, s, mig, L, N);

        }
            
        fclose(datafile);
    }




	if (t == tfinal){
		printf("Simulation finished without fixation.\n");
		return 1;
	}
	
	// if (n[0] == 0){
	// 	printf("Allele went extinct.\n");
	// 	return 2;
	// }


    return 0;
}
