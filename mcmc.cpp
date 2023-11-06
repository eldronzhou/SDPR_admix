#include "mcmc.h"
#include <cmath>
#include <random>
#include "gsl/gsl_randist.h"
#include "gsl/gsl_statistics_double.h"
#include "gsl/gsl_blas.h"
#include <fstream>
#include <sstream>
#include <string.h>

using std::log; using std::exp; using std::sqrt;
using std::cout; using std::endl;

#define square(x) ((x)*(x))

void init_state(Dat *dat, MCMC_state *state) {
    state->n_cluster = state->M_cluster[0] + state->M_cluster[1] + \
		       state->M_cluster[2] + state->M_cluster[3] + state->M_cluster[4];
    state->beta1 = (double *) calloc(dat->n_snp, sizeof(double));
    state->beta2 = (double *) calloc(dat->n_snp, sizeof(double));
    state->residual = (double *) malloc(sizeof(double)*dat->n_ind);
    state->G = (double *) calloc(dat->n_ind, sizeof(double));
    

    dat->y = (double *) malloc(sizeof(double)*dat->n_ind);
    for (size_t i=0; i<dat->n_ind; i++) {
	dat->y[i] = dat->pheno[i];
    }

    for (size_t i=0; i<dat->n_ind; i++) {
	state->residual[i] = dat->pheno[i];
    }

    state->assgn = (int *) malloc(sizeof(int)*dat->n_snp);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(1, state->n_cluster);
    
    //std::ifstream infile("../effect_size/sim_1.txt");
    for (size_t i=0; i<dat->n_snp; i++) {
	state->assgn[i] = dist(gen);
	//infile >> state->assgn[i];
    }

    state->suffstats = (int *) calloc(state->n_cluster, sizeof(int));
    for (size_t i=0; i<state->n_cluster; i++) {
	if (i > dat->n_snp) {
	    break;
	}
	state->suffstats[state->assgn[i]]++;
    }
    state->sumsq = (double *) calloc(state->n_cluster, sizeof(double));
  
    state->cluster_var = (double *) malloc(sizeof(double)*state->n_cluster);
    for (size_t i=1; i<state->n_cluster; i++) {
	std::gamma_distribution<> rgamma(state->suffstats[i]/2.0+state->a0k, \
		1.0/state->b0k);
	state->cluster_var[i] = 1.0/rgamma(gen);
    }

    state->pi = (double *) malloc(sizeof(double)*state->n_cluster);
   
    int population[5] = {1, state->M_cluster[1]+1, \
 	    state->M_cluster[1]+state->M_cluster[2]+1, \
	    state->M_cluster[1]+state->M_cluster[2]+state->M_cluster[3]+1, \
            state->M_cluster[1]+state->M_cluster[2]+state->M_cluster[3]+state->M_cluster[4]+1}; 
    
    for (size_t j=1; j<5; j++) {
	for (int i=population[j-1]; i<population[j]; i++) {
	    state->pi[i] = 1.0/state->M_cluster[j];
	}
    }
    
    state->log_p = (double *) malloc(sizeof(double)*state->n_cluster);
    state->V = (double *) malloc(sizeof(double)*state->n_cluster);
    state->aj1 = (double *) malloc(sizeof(double)*state->n_cluster);
    state->aj2 = (double *) malloc(sizeof(double)*state->n_cluster);
    state->mu_j1 = (double *) malloc(sizeof(double)*state->n_cluster);
    state->mu_j2 = (double *) malloc(sizeof(double)*state->n_cluster);
    state->cj = (double *) malloc(sizeof(double)*state->n_cluster);
}

void destroy_state(MCMC_state *state) {
    free(state->beta1);
    free(state->beta2);
    free(state->residual);
    free(state->G);
    free(state->assgn);
    free(state->suffstats);
    free(state->sumsq);
    free(state->cluster_var);
    free(state->pi);
    free(state->log_p);
    free(state->V);
    free(state->aj1);
    free(state->aj2);
    free(state->mu_j1);
    free(state->mu_j2);
    free(state->cj);
}

void mcmc(Dat *dat, std::string out_path, int iter, int burn, double maf, double rho_0) {
    MCMC_state state;
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    std::random_device rd;
    std::mt19937 gen(rd());
   
    init_state(dat, &state);
    state.rho = rho_0;

    double *res_beta1 = (double *) calloc(dat->n_snp, sizeof(double));
    double *res_beta2 = (double *) calloc(dat->n_snp, sizeof(double));

    // set up OLS env
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(dat->n_ind, dat->n_cov);
    gsl_vector *y = gsl_vector_alloc(dat->n_ind);
    gsl_matrix *W = gsl_matrix_alloc(dat->n_ind, dat->n_cov);
    if (dat->n_cov != 0) {
	for (size_t i=0; i<dat->n_ind; i++) {
	    for (size_t j=0; j<dat->n_cov; j++) {
		gsl_matrix_set(W, i, j, dat->covar[j][i]);
	    }
	}
    }
    gsl_matrix *cov = gsl_matrix_alloc(dat->n_cov, dat->n_cov);
    gsl_vector *alpha = gsl_vector_alloc(dat->n_cov);
    double chisq = 0;

    for (size_t n_mcmc=0; n_mcmc<iter; n_mcmc++) {
	int n_total = iter - 1 -  burn;

	// OLS
	if (dat->n_cov != 0) {
	    for (size_t i=0; i<dat->n_ind; i++) {
		gsl_vector_set(y, i, dat->y[i] - state.eta * state.G[i]);
	    }
	    gsl_multifit_linear(W, y, alpha, cov, &chisq, work);
	    
	    gsl_blas_dgemv(CblasNoTrans, 1.0, W, alpha, 0, y);

	    // update pheno
	    for (size_t i=0; i<dat->n_ind; i++) {
		dat->pheno[i] = dat->y[i] - gsl_vector_get(y, i);
		state.residual[i] = dat->pheno[i] - state.eta * state.G[i];
	    }
	}

	//std::ifstream infile1("../effect_size/Scene1/AFR/sim_1.txt");
	//std::ifstream infile2("../effect_size/Scene1/EUR/sim_1.txt");

	for (size_t i=0; i<dat->n_snp; i++) {
	    /*if (dat->maf1[i] < maf || dat->maf2[i] < maf) {
		state.assgn[i] = 0;
		continue;
	    }*/

	    // sample z_j
	    // compute residual related terms
	    double bj1 = 0, bj2 = 0;
	    
	    for (size_t j=0; j<dat->n_ind; j++) {
		bj1 += state.residual[j] * dat->geno1[i][j];
		bj2 += state.residual[j] * dat->geno2[i][j];
	    }

	    bj1 += state.eta * ( dat->geno1_sq[i] * state.beta1[i] + dat->geno12_prod[i] * state.beta2[i] );
	    bj2 += state.eta * ( dat->geno2_sq[i] * state.beta2[i] + dat->geno12_prod[i] * state.beta1[i] );
	    bj1 *= state.eta / state.sigmae2;
	    bj2 *= state.eta / state.sigmae2;
	    
	    // 0 effects for ancestry 1 & 2
	    state.log_p[0] = log(state.pi_pop[0] + 1e-40);
	    
	    // ancestry 1 specific terms
	    for (int k=1; k<state.M_cluster[1]+1; k++) {
		state.log_p[k] = -.5 * log(square(state.eta) * dat->geno1_sq[i] * state.cluster_var[k] / state.sigmae2 + 1) \
		                  + .5 * square(bj1) / (square(state.eta) * dat->geno1_sq[i] / \
					 state.sigmae2 + 1.0 / state.cluster_var[k] ) \
				  + log(state.pi_pop[1] * state.pi[k] + 1e-40); 
	    }
	    
	    // ancestry 2 specific terms
	    for (int k=state.M_cluster[1]+1; k<state.M_cluster[1]+state.M_cluster[2]+1; k++) {
		state.log_p[k] = -.5 * log(square(state.eta) * dat->geno2_sq[i] * state.cluster_var[k] / state.sigmae2 + 1) \
				 + .5 * square(bj2) / (square(state.eta) * dat->geno2_sq[i] / \
					 state.sigmae2 + 1.0 / state.cluster_var[k] ) \
				 + log(state.pi_pop[2] * state.pi[k] + 1e-40);
	    }
	     
	    // ancestry 1 & 2 shared terms
	    for (int k=state.M_cluster[1]+state.M_cluster[2]+1; k<state.n_cluster; k++) {
		
		double rho = 0, pi_pop = 0;
		if (k < state.M_cluster[1]+state.M_cluster[2]+state.M_cluster[3]+1) {
		    // correlated effect sizes
		    rho = state.rho;
		    pi_pop = state.pi_pop[3];
		}
		else {
		    // same effect sizes 
		    rho = 0.9999;
		    pi_pop = state.pi_pop[4];
		}

		state.aj1[k] = dat->geno1_sq[i] / (2*state.sigmae2) * square(state.eta) \
			       + .5/( (1 - square(rho)) * state.cluster_var[k]);

		state.aj2[k] = dat->geno2_sq[i] / (2*state.sigmae2) * square(state.eta) \
			       + .5/( (1 - square(rho)) * state.cluster_var[k]);

		state.cj[k] = - dat->geno12_prod[i] / state.sigmae2 * square(state.eta) \
			       + rho / ( (1 - square(rho)) * state.cluster_var[k]);

		state.mu_j1[k] = ( 2 * state.aj2[k] * bj1 + state.cj[k] * bj2 ) / \
				 ( 4 * state.aj1[k] * state.aj2[k] - square(state.cj[k]) );
		
		state.mu_j2[k] = ( 2 * state.aj1[k] * bj2 + state.cj[k] * bj1 ) / \
				 ( 4 * state.aj1[k] * state.aj2[k] - square(state.cj[k]) );
		
		state.log_p[k] = -.5 * log( 4 * state.aj1[k] * state.aj2[k] - square(state.cj[k]) ) - \
				 log( state.cluster_var[k] ) - \
				 .5 * log( 1 - square(rho) ) + \
				 state.aj1[k] * square(state.mu_j1[k]) + \
				 state.aj2[k] * square(state.mu_j2[k]) - \
				 state.cj[k] * state.mu_j1[k] * state.mu_j2[k] + \
				 log(pi_pop * state.pi[k] + 1e-40);
	    }

	    double max = state.log_p[0], log_exp_sum = 0;
	    for (size_t k=1; k<state.n_cluster; k++) {
		// only keep ancestry 2 for ancestry 2 specific variants
		if (dat->maf1[i] < maf) {
		    if (k>=state.M_cluster[1]+1 && k<state.M_cluster[1]+state.M_cluster[2]+1) {
			if (max < state.log_p[k]) {
			    max = state.log_p[k];
			}
		    }
		}
		// only keep ancestry 1 for ancestry 1 specific variants
		else if (dat->maf2[i] < maf) {
		    if (k>=1 && k<state.M_cluster[1]+1) {
			if (max < state.log_p[k]) {
			    max = state.log_p[k];
			}
		    }
		}
		else {
		    if (max < state.log_p[k]) {
			max = state.log_p[k];
		    }
		}
	    }

	    log_exp_sum = 0;
	    for (size_t k=0; k<state.n_cluster; k++) {
		// only keep ancestry 2 for ancestry 2 specific variants
		if (dat->maf1[i] < maf) {
		    if (k>=state.M_cluster[1]+1 && k<state.M_cluster[1]+state.M_cluster[2]+1) {
			log_exp_sum += exp(state.log_p[k] - max);
		    }
		}
		// only keep ancestry 1 for ancestry 1 specific variants
		else if (dat->maf2[i] < maf) {
		    if (k>=1 && k<state.M_cluster[1]+1) {
			log_exp_sum += exp(state.log_p[k] - max);
		    }
		}
		else {
		    log_exp_sum += exp(state.log_p[k] - max);
		}
	    }
	    log_exp_sum = max + log(log_exp_sum);

	    state.assgn[i] = state.n_cluster - 1;
	    double rnd = gsl_rng_uniform(r);
	    for (size_t k=0; k<state.n_cluster-1; k++) {
		if (dat->maf1[i] < maf) {
		    // only keep ancestry 2 for ancestry 2 specific variants
		    state.assgn[i] = 0;
		    break;
		    if (k>=1 && k<state.M_cluster[1]+1) {
			continue;
		    }
		    else if (k >= state.M_cluster[1]+state.M_cluster[2]+1) {
			state.assgn[i] = state.M_cluster[1]+state.M_cluster[2];
			break;
		    }
		}
		else if (dat->maf2[i] < maf) {
		    // only keep ancestry 1 for ancestry 1 specific variants
		    state.assgn[i] = 0;
		    break;
		    if (k>=state.M_cluster[1]+1 && k<state.M_cluster[1]+state.M_cluster[2]+1) {
			continue;
		    }
		    else if (k >= state.M_cluster[1]+state.M_cluster[2]+1) {
			state.assgn[i] = state.M_cluster[1];
			break;
		    }
		}
		rnd -= exp(state.log_p[k] - log_exp_sum);
		if (rnd < 0) {
		    state.assgn[i] = k;
		    break;
		}
	    }
	    
	    double beta1_old = state.beta1[i];
	    double beta2_old = state.beta2[i];

	    //double tmp;
	    std::string id;

	    if (state.assgn[i] == 0) {
		state.beta1[i] = 0;
		state.beta2[i] = 0;
		//infile1 >> id >> tmp;
		//infile2 >> id >> tmp;
	    }
	    else if (state.assgn[i] >= 1 && state.assgn[i] < state.M_cluster[1]+1) {
		if (dat->geno1_sq[i] != 0) {
		    double C = dat->geno1_sq[i] * square(state.eta) / state.sigmae2 + \
			         1.0 / state.cluster_var[state.assgn[i]];

		    state.beta1[i] = 1 / sqrt(C) *  gsl_ran_gaussian(r, 1.0) \
				     + bj1 / C;
		    state.beta2[i] = 0; 
		}
		else {
		    state.beta1[i] = 0;
		    state.beta2[i] = 0;
		}
		//infile1 >> id >> tmp;
		//infile2 >> id >> tmp;
	    }
	    else if (state.assgn[i] >= state.M_cluster[1]+1 && state.assgn[i] < state.M_cluster[1]+state.M_cluster[2]+1) {
		if (dat->geno2_sq[i] != 0) {
		    double C = dat->geno2_sq[i] * square(state.eta) / state.sigmae2 + \
			       1.0 / state.cluster_var[state.assgn[i]];

		    state.beta1[i] = 0;
		    state.beta2[i] = 1 / sqrt(C) *  gsl_ran_gaussian(r, 1.0) \
				     + bj2 / C;
		}
		else {
		    state.beta1[i] = 0;
		    state.beta2[i] = 0;
		}
		//infile1 >> id >> tmp;
		//infile2 >> id >> tmp;
	    }
	    else {
		double sigma_x = sqrt(2*state.aj2[state.assgn[i]] / \
			(4 * state.aj1[state.assgn[i]] * state.aj2[state.assgn[i]] \
			 - square(state.cj[state.assgn[i]])));
		double sigma_y = sqrt(2*state.aj1[state.assgn[i]] / \
			(4 * state.aj1[state.assgn[i]] * state.aj2[state.assgn[i]] \
			 - square(state.cj[state.assgn[i]])));
		double rho = state.cj[state.assgn[i]] / \
			     (2 * sqrt(state.aj1[state.assgn[i]] * state.aj2[state.assgn[i]]));

		gsl_ran_bivariate_gaussian(r, sigma_x, sigma_y, rho, &state.beta1[i], &state.beta2[i]);
		
		state.beta1[i] += state.mu_j1[state.assgn[i]];
		state.beta2[i] += state.mu_j2[state.assgn[i]];

		//infile1 >> id >> tmp;
		//infile2 >> id >> tmp;
	    }

	    /*std::string id;
	    infile1 >> id >> state.beta1[i];
	    infile2 >> id >> state.beta2[i];*/

	    // update genetic effects and residual	    
	    for (size_t j=0; j<dat->n_ind; j++) {
		state.G[j] += dat->geno1[i][j] * (state.beta1[i] - beta1_old) + \
		       dat->geno2[i][j] * (state.beta2[i] - beta2_old);

		state.residual[j] = dat->pheno[j] - state.eta * state.G[j];
	    }
	}
	
	// update suffstats
	for (size_t k=0; k<state.n_cluster; k++) {
	    state.suffstats[k] = 0;
	    state.sumsq[k] = 0;
	}
	
	for (size_t i=0; i<dat->n_snp; i++) {
	    state.suffstats[state.assgn[i]]++;
	    if (state.assgn[i] < state.M_cluster[1]+1) {
		state.sumsq[state.assgn[i]] += square(state.beta1[i]) / 2.0;
	    } 
	    else if (state.assgn[i] >= state.M_cluster[1]+1 && state.assgn[i] < state.M_cluster[1]+state.M_cluster[2]+1) {
		state.sumsq[state.assgn[i]] += square(state.beta2[i]) / 2.0;
	    }
	    else {
		double rho = 0;
		           
		if (state.assgn[i] < state.M_cluster[1]+state.M_cluster[2]+state.M_cluster[3]+1) {
		    rho = state.rho;
		}
		else {
		    rho = 0.9999;
		}
		state.sumsq[state.assgn[i]] += ( square(state.beta1[i]) + square(state.beta2[i]) - \
			2 * rho * state.beta1[i] * state.beta2[i] ) / (2 * (1 - square(rho)) );
	    }
	}

	// sample cluster_var
	for (size_t k=1; k<state.n_cluster; k++) {
	    double a = state.suffstats[k] / 2.0 + state.a0k;
	    double b = 1.0 / (state.sumsq[k] + state.b0k);
	    state.cluster_var[k] = 1.0 / gsl_ran_gamma(r, a, b);
	}
	
	// sample V
	int population[5] = {1, state.M_cluster[1]+1, \
	                     state.M_cluster[1]+state.M_cluster[2]+1, \
	                     state.M_cluster[1]+state.M_cluster[2]+state.M_cluster[3]+1, \
			     state.M_cluster[1]+state.M_cluster[2]+state.M_cluster[3]+state.M_cluster[4]+1}; 
	
	double *a = (double *) calloc(state.n_cluster, sizeof(double));

	a[(size_t)state.n_cluster-2] = state.suffstats[(size_t)state.n_cluster-1];
	for (int k=state.n_cluster-3; k>0; k--) {
	    if (k == state.M_cluster[1] || k == state.M_cluster[1]+state.M_cluster[2] || k == state.M_cluster[1]+state.M_cluster[2]+state.M_cluster[3]) {
		a[k] = 0;
	    }
	    else { 
		a[k] = state.suffstats[k+1] + a[k+1];
	    }
	}

	for (size_t j=1; j<5; j++) {
	    for (int k=population[j-1]; k<population[j]; k++) {
		state.V[k] = gsl_ran_beta(r, 1 + state.suffstats[k], \
			    state.alpha[j-1] + a[k]);
	    }
	}
	state.V[state.M_cluster[1]] = 1;
	state.V[state.M_cluster[1]+state.M_cluster[2]] = 1;
	state.V[state.M_cluster[1]+state.M_cluster[2]+state.M_cluster[3]] = 1;
	state.V[state.M_cluster[1]+state.M_cluster[2]+state.M_cluster[3]+state.M_cluster[4]] = 1;


	// update p
	a[1] = 1 - state.V[1]; 
	state.pi[1] = state.V[1];
	for (int k=2; k<state.M_cluster[1]+1; k++) {
	    a[k] = a[k-1] * (1 - state.V[k]);
	    state.pi[k] = a[k-1] * state.V[k];
	    if (state.V[k] == 1) {
		for (int m=k+1; m<state.M_cluster[1]+1; m++) {
		    a[m] = 0;
		    state.pi[m] = 0;
		}
		break;
	    }
	}
	if (state.pi[state.M_cluster[1]] < 0) {
	    state.pi[state.M_cluster[1]] = 0;
	}

	a[state.M_cluster[1]+1] = 1 - state.V[state.M_cluster[1]+1];
	state.pi[state.M_cluster[1]+1] = state.V[state.M_cluster[1]+1];
	for (int k=state.M_cluster[1]+2; k<state.M_cluster[1]+state.M_cluster[2]+1; k++) {
	    a[k] = a[k-1] * (1 - state.V[k]);
	    state.pi[k] = a[k-1] * state.V[k];
	    if (state.V[k] == 1) {
		for (int m=k+1; m<state.M_cluster[1]+state.M_cluster[2]+1; m++) {
		    a[m] = 0;
		    state.pi[m] = 0;
		}
		break;
	    }
	}
	if (state.pi[state.M_cluster[1]+state.M_cluster[2]] < 0) {
	    state.pi[state.M_cluster[1]+state.M_cluster[2]] = 0;
	}

	a[state.M_cluster[1]+state.M_cluster[2]+1] = 1 - state.V[state.M_cluster[1]+state.M_cluster[2]+1];
	state.pi[state.M_cluster[1]+state.M_cluster[2]+1] = state.V[state.M_cluster[1]+state.M_cluster[2]+1];
	for (int k=state.M_cluster[1]+state.M_cluster[2]+2; k<state.M_cluster[1]+state.M_cluster[2]+state.M_cluster[3]+1; k++) {
	    a[k] = a[k-1] * (1 - state.V[k]);
	    state.pi[k] = a[k-1] * state.V[k];
	    if (state.V[k] == 1) {
		for (int m=k+1; m<state.n_cluster; m++) {
		    a[m] = 0;
		    state.pi[m] = 0;
		}
		break;
	    }
	}
	if (state.pi[state.M_cluster[1]+state.M_cluster[2]+state.M_cluster[3]] < 0) {
	    state.pi[state.M_cluster[1]+state.M_cluster[2]+state.M_cluster[3]] = 0;
	}

	a[state.M_cluster[1]+state.M_cluster[2]+state.M_cluster[3]+1] = 1 - state.V[state.M_cluster[1]+state.M_cluster[2]+state.M_cluster[3]+1];
	state.pi[state.M_cluster[1]+state.M_cluster[2]+state.M_cluster[3]+1] = state.V[state.M_cluster[1]+state.M_cluster[2]+state.M_cluster[3]+1];
	for (int k=state.M_cluster[1]+state.M_cluster[2]+state.M_cluster[3]+2; k<state.n_cluster; k++) {
	    a[k] = a[k-1] * (1 - state.V[k]);
	    state.pi[k] = a[k-1] * state.V[k];
	    if (state.V[k] == 1) {
		for (int m=k+1; m<state.n_cluster; m++) {
		    a[m] = 0;
		    state.pi[m] = 0;
		}
		break;
	    }
	}
	if (state.pi[(size_t)state.n_cluster-1] < 0) {
	    state.pi[(size_t)state.n_cluster-1] = 0;
	}

	free(a);

	// sample pi_pop
	int null_anc1 = 0, null_anc2 = 0;        
	int nonull_anc1 = 0, nonull_anc2 = 0;
	
	double alpha[5] = {1, 1, 1, 1, 1};
	alpha[0] += state.suffstats[0];
	for (size_t j=1; j<5; j++) {
	    for (int k=population[j-1]; k<population[j]; k++) {
		alpha[j] += state.suffstats[k];
	    }
	}
	alpha[0] = alpha[0] - null_anc1 - null_anc2 ;
	alpha[1] -= nonull_anc1;
	alpha[2] -= nonull_anc2;
	gsl_ran_dirichlet(r, 5, alpha, state.pi_pop);
	//state.pi_pop[1] = 0; state.pi_pop[2] = 0;

	// sample alpha
	double sum = 0, m = 0;
	for (int k=1; k<state.M_cluster[1]+1; k++) {
	    if (state.V[k] != 0) {
		sum += log(1 - state.V[k]);
		m++;
	    }
	}
	state.alpha[0] = gsl_ran_gamma(r, .1+m-1, 1.0/(.1-sum));

	sum = 0, m = 0;
	for (int k=state.M_cluster[1]+1; k<state.M_cluster[1]+state.M_cluster[2]+1; k++) {
	    if (state.V[k] != 0) {
		sum += log(1 - state.V[k]);
		m++;
	    }
	}
	state.alpha[1] = gsl_ran_gamma(r, .1+m-1, 1.0/(.1-sum));

	sum = 0, m = 0;
	for (int k=state.M_cluster[1]+state.M_cluster[2]+1; k<state.M_cluster[1]+state.M_cluster[2]+state.M_cluster[3]+1; k++) {
	    if (state.V[k] != 0) {
		sum += log(1 - state.V[k]);
		m++;
	    }
	}
	state.alpha[2] = gsl_ran_gamma(r, .1+m-1, 1.0/(.1-sum)); 

	sum = 0, m = 0;
	for (int k=state.M_cluster[1]+state.M_cluster[2]+state.M_cluster[3]+1; k<state.n_cluster; k++) {
	    if (state.V[k] != 0) {
		sum += log(1 - state.V[k]);
		m++;
	    }
	}
	state.alpha[3] = gsl_ran_gamma(r, .1+m-1, 1.0/(.1-sum));

	// sample eta and compute h2
	double G2 = 0, YG_prod = 0;
	for (size_t j=0; j<dat->n_ind; j++) {
	    YG_prod += state.G[j] * dat->pheno[j];
	    G2 += square(state.G[j]);
	}
	state.eta = gsl_ran_ugaussian(r) * sqrt(1.0/(G2 + 1e-6)) + \
		    YG_prod / (G2 + 1e-6);

	// avoid numeric issues if eta is too small
	if (state.eta < 1e-20) {
	    for (size_t j=0; j<state.n_cluster; j++) {
		state.cluster_var[j] *= square(state.eta);
	    }
	    for (size_t i=0; i<dat->n_snp; i++) {
		state.beta1[i] *= state.eta;
		state.beta2[i] *= state.eta;
	    }
	    for (size_t j=0; j<dat->n_ind; j++) {
		state.G[j] *= state.eta;
	    }
	    state.eta = 1;
	}

	double residual2 = 0;
	for (size_t j=0; j<dat->n_ind; j++) {
	    state.residual[j] = dat->pheno[j] - state.eta * state.G[j];
	    residual2 += square(state.residual[j]);
	}

	// sample sigmae2
	double ae = dat->n_ind / 2.0 + 0.1;
	double be = 1 / (residual2 / 2.0 + 0.1);
	state.sigmae2 = 1.0 / gsl_ran_gamma(r, ae, be);
	
	// output h2
	double h2 = square(state.eta) * gsl_stats_variance(state.G, 1, dat->n_ind) / gsl_stats_variance(dat->pheno, 1, dat->n_ind);

	// record the beta
	if (n_mcmc > burn) {
	    for (size_t i=0; i<dat->n_snp; i++) {
		res_beta1[i] += state.eta * state.beta1[i] / n_total;
		res_beta2[i] += state.eta * state.beta2[i] / n_total;
	    }
	}

	double max_beta1 = 0, max_beta2 = 0;
	for (size_t i=0; i<dat->n_snp; i++) {
	    if (state.eta * state.beta1[i] > max_beta1) {
		max_beta1 = state.eta * state.beta1[i];
	    }
	    if (state.eta * state.beta2[i] > max_beta2) {
		max_beta2 = state.eta * state.beta2[i];
	    }
	}

	if (n_mcmc % 1 == 0) {
	    std::ostringstream eta_out;
	    eta_out << state.eta;
	    std::string strobj = eta_out.str();

	cout << std::to_string(h2) << " " << std::to_string(state.pi_pop[0]) \
	    << " " << std::to_string(state.pi_pop[1]) << " " << \
	    std::to_string(state.pi_pop[2]) << " " << std::to_string(state.pi_pop[3]) \
	    << " " << std::to_string(state.pi_pop[4]) << " max beta1: " << std::to_string(max_beta1) \
	    <<" max beta2: " << std::to_string(max_beta2) << endl; \
	}
    }   
    
    destroy_state(&state);
    gsl_rng_free(r);

    std::ofstream out(out_path);
    for (size_t i=0; i<dat->n_snp; i++) {
	out << dat->chr[i] << "\t" << dat->pos[i] << "\t" << \
	    dat->id[i] << "\t" << dat->ref[i] << "\t" << \
	    dat->alt[i] << "\t" << res_beta1[i] << \
	    "\t" << res_beta2[i] << endl;
    }
    out.close();
    free(res_beta1);
    free(res_beta2);

    gsl_vector_free(y);
    gsl_vector_free(alpha);
    gsl_matrix_free(W);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free(work);
}

/*int main(int argc, char *argv[]) {
    Dat dat;

    double rho = 0;

    std::string pheno_path, geno1_path, \
	geno2_path, vcf_path, msp_path, out_path, covar_path;

    int i = 1, iter = 1000, burn = 500;
    while (i < argc) {
	if (strcmp(argv[i], "-pheno") == 0) {
	    pheno_path = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-geno1") == 0) {
	    geno1_path = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-geno2") == 0) {
	    geno2_path = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-vcf") == 0) {
	    vcf_path = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-msp") == 0) {
	    msp_path = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-rho") == 0) {
	    rho = std::stoi(argv[i+1]);
	    i += 2;
	}
	else if (strcmp(argv[i], "-iter") == 0) {
	    iter = std::stoi(argv[i+1]);
	    i +=2;
	}
	else if (strcmp(argv[i], "-burn") == 0) {
	    burn = std::stoi(argv[i+1]);
	    i += 2;
	}
	else if (strcmp(argv[i], "-rho") == 0) {
            rho = std::stod(argv[i+1]);
            i += 2;
        }
	else if (strcmp(argv[i], "-out") == 0) {
	    out_path = argv[i+1];
	    i += 2;
	}
	else if (strcmp(argv[i], "-covar") == 0) {
	    covar_path = argv[i+1];
	    i += 2;
	}
	else {
	    cout << "Invalid option: " << argv[i] << endl;
	    return 0;
	}
    }

    get_size_vcf(pheno_path.c_str(), vcf_path.c_str(), &dat);
    
    read_lanc(vcf_path.c_str(), msp_path.c_str(), &dat);

    read_pheno(pheno_path.c_str(), &dat);

    if (!covar_path.empty()) {
	read_cov(covar_path.c_str(), &dat);
    }

    //linear(&dat, out_path.c_str());
   
    double maf = 0;

    check_maf(&dat, maf);
    
    prod_geno(&dat);

    maf = 0.05;
    mcmc(&dat, out_path.c_str(), iter, burn, maf, rho);

    for (size_t i=0; i<dat.n_snp; i++) {
	free(dat.geno1[i]);
	free(dat.geno2[i]);
    }
    free(dat.geno1);
    free(dat.geno2);
    if (!covar_path.empty()) {
	free(dat.covar);
    }
    free(dat.maf1);
    free(dat.maf2);
    free(dat.n_anc1);
    free(dat.n_anc2);
    free(dat.y);
    free(dat.pheno);
    free(dat.geno1_sq);
    free(dat.geno2_sq);
    free(dat.geno12_prod);
    return 0;
}*/
