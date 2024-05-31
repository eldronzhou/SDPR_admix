#include "parse_gen.h"
#include "regress.h"

typedef struct {
    double eta = 1.0;
    double rho = 0.93;
    //double sigmae2 = 0.14;
    double sigmae2 = 1.0;
    double a0k = 0.5, b0k = 0.5;
    int M_cluster[4] = {1, 500, 500, 1000};
    double n_cluster;
    double *beta1;
    double *beta2;
    double *residual;
    double *G;
    int *assgn;
    double *cluster_var;
    double *pi;
    double *log_p;
    double *V;
    int *suffstats;
    double *sumsq;
    //double pi_pop[4] = {0.8, 0.05, 0.05, 0.1};
    double pi_pop[4] = {0.25, 0.25, 0.25, 0.25};
    double alpha[3] = {1.0, 1.0, 1.0};
    double *aj1;
    double *aj2;
    double *cj;
    double *mu_j1;
    double *mu_j2;
    double bj1 = 0;
    double bj2 = 0;
} MCMC_state;

typedef struct {
    double *beta1;
    double *beta2;
    double h2;
} MCMC_samples;

void mcmc(Dat *dat, std::string out_path, int iter, int burn, double maf, double rho, int thread);
