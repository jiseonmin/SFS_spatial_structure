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

// run by typing : gcc plain_sweep.c -o plain_sweep -lgsl
// ./plain_sweep L N s m tc n_sim dist_edge ranseed outfile



//Todo : when beneficial mutation is fixed in 10% of the demes, save n[i] and t. Get the right most deme with a beneficial allele, choose a random individual among them, and give a neutral mutation and track it.
// When the neutral mutation goes extinct, set n[i] and t to the saved values and start tracking again. Save n[i] and n_neutral[i]'s. 
// Later, there should be one more input value (argv) that sets how far out from the center of mass (i_cm = sum(n[i] * i) / N ) to put in a neutral mutation.


// n0[i] = In i-th deme, # of individuals with a beneficial mutation but no neutral marker
// n1[i] = beneficial mutation + neutral marker.

// Global variables
const gsl_rng_type * T;
gsl_rng * R;

void next_gen(unsigned int n0[], unsigned int n1[], double s, double mig, int L, unsigned int N) {

    double x0[L]; // frequencies
    double x0mig[L]; // after migration
    double x1[L];
    double x1mig[L]; 
    int i;
    unsigned int n[L];

	// Initialize frequencies:
	for (i = 0; i < L; i++){
		x0[i] = (double)n0[i] / N;
        x1[i] = (double)n1[i] / N;
	}
	// Migration:
	for (i = 1; i < L - 1; i++){
		x0mig[i] = x0[i] + mig * (0.5 * (x0[i-1] + x0[i+1]) - x0[i]);
		x1mig[i] = x1[i] + mig * (0.5 * (x1[i-1] + x1[i+1]) - x1[i]);
	}
	x0mig[0] = x0[0] + mig * 0.5 * (x0[1] - x0[0]);
	x0mig[L-1] = x0[L-1] + mig * 0.5 * (x0[L-2] - x0[L-1]);
	x1mig[0] = x1[0] + mig * 0.5 * (x1[1] - x1[0]);
	x1mig[L-1] = x1[L-1] + mig * 0.5 * (x1[L-2] - x1[L-1]);

    // Sampling and selection within demes
    for (i = 0; i < L; i++) {
		n[i] = gsl_ran_binomial(R, (x0mig[i] + x1mig[i]) + s * (x0mig[i] + x1mig[i]) * (1 - x0mig[i] - x1mig[i]), N);
        n0[i] = gsl_ran_binomial(R, x0mig[i] / (x0mig[i] + x1mig[i]), n[i]);
        n1[i] = n[i] - n0[i];
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
    unsigned int N, L, t, tc; // deme size, number of demes, time, number of generation when we add a neutral mutation
    int dist_edge; // distance from the front to seed the neutral mutation
    FILE *datafile, *paramfile, *lineagefile;
    char *outfile;
    char *outfile1 = malloc(1000);
    char *outfile2 = malloc(1000);

    unsigned int i, j, ntot, i_end, i_end_copy, NSIM, nsim, n1tot, n1tot2, ntottc;
    //Initialize variables:
    if (argc != 10) {
        printf("Usage: L N s m tc n_sim dist_edge ranseed outfile\n");
        return 1; 
    }
    
    j = 1;
    L = (unsigned int) roundint(atof(argv[j++]));
    N = (unsigned int) roundint(atof(argv[j++]));
    s = atof(argv[j++]);
    mig = atof(argv[j++]);
    tc = (unsigned int) roundint(atof(argv[j++]));
    NSIM = atof(argv[j++]);
    dist_edge = atof(argv[j++]);
    SEED = atof(argv[j++]);
	outfile = argv[j++];
	strcpy(outfile1,outfile);
	strcpy(outfile2,outfile);


    // Print out variables to parameter file:
    paramfile = fopen(strcat(outfile1, "_params.txt"), "w");
    fprintf(paramfile, 
    "L = %u\nN = %u\ns = %f\nm = %f\ntc = %u\nseed = %ld\nnsim = %u\ndistance from the front =%d"
    , L, N, s, mig, tc, SEED, NSIM, dist_edge);
    fclose(paramfile);
	
    // gsl random setup:
    gsl_rng_env_setup();
    T = gsl_rng_default;
    R = gsl_rng_alloc (T);
    gsl_rng_set(R,SEED);


    // Initialize population:
    unsigned int n0[L];
    unsigned int n1[L];

    // We will save n0 and n1 at time t_c and add a neutral mutation. If the neutral mutation goes extinct, we go back to t = t_c and try again.
    unsigned int n0_tc[L];
    unsigned int n1_tc[L];

    // Leftmost demes fixed for sweeping allele
    // Fill enough demes so that it is very unlikely to die out:
    if (1 + 10 / (N * s) > L / 2){
    	printf("Warning: meta-population is too small for selection to be strong.\n");
    	return 1;
    }
    for (i = 0; i < 1 + 10 / (N * s); i++){
    	n0[i] = N;
        n1[i] = 0;
    }
    // All other demes fixed for ancestral allele:
    for (i = (int)(1 + 10 / (N * s)); i < L; i++) {
        n0[i] = 0;
        n1[i] = 0;
    }
    
    //Open the datafile for recording:
    datafile = fopen(strcat(outfile, ".txt"), "w");

    // Run until tc 
    t = 0;
    for (t = 0; t < tc; t++) {
		// Record the status of the population and check for fixation:

		ntot = 0;
        n1tot = 0;
        i_end = -1;
		for (i = 0; i < L; i++){
			fprintf(datafile, " %d", n0[i]);
			ntot += n0[i] + n1[i];
            n1tot += n1[i];
            if (n0[i] > 0) {i_end++;}
		}
		fprintf(datafile, "\n");
        
        // Stop if one of the alleles is fixed or extinct:
		if ((ntot == 0) || (ntot == N * L)) {break;}
		
		// Evolve for one generation
        next_gen(n0, n1, s, mig, L, N);

    }
    // save n0 and n1 at t = tc
    for (i = 0; i < L; i ++){
        n0_tc[i] = n0[i];
        n1_tc[i] = n1[i];
    }
    ntottc = ntot;
    fclose(datafile);

    
	if (n0[0] == 0){
		printf("Allele went extinct.\n");

		return 2;
    }
    
    else {
        lineagefile = fopen(strcat(outfile2, "_lineages.txt"), "w");
        n1_tc[i_end + dist_edge] = 1;
        n0_tc[i_end + dist_edge] -= 1;
        
        for (nsim = 0; nsim < NSIM + 1; nsim ++){
            fprintf(lineagefile, " %d", nsim);
            fprintf(lineagefile, "\n");
            n1tot = 1;
            for (i = 0; i < L; i ++){
                n0[i] = n0_tc[i];
                n1[i] = n1_tc[i];
            }
            t = tc;
            ntot = ntottc;
            // Track the lineage until it goes extinct or beneficial mutation fixes. 
            while ((n1tot > 0) && (ntot < N * L)){
                t++;
                ntot = 0;
                n1tot2 = 0;
                for (i = 0; i < L; i++){
                    // fprintf(lineagefile, " %d", n1[i]);
                    ntot += n0[i] + n1[i];
                    n1tot2 += n1[i];
                }
                n1tot = n1tot2;
                fprintf(lineagefile, " %d", n1tot);
                                
                // Evolve for one generation
                next_gen(n0, n1, s, mig, L, N);

                }
            
            fprintf(lineagefile, "\n");



        }
        fclose(lineagefile);

    }

	


    return 0;
}
