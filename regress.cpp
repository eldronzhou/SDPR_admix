#include "regress.h"

using std::cout; using std::endl; 

void linear(Dat *dat, std::string out_path) {
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(dat->n_ind, dat->n_cov+2);

    gsl_vector *y = gsl_vector_alloc(dat->n_ind);

    for (size_t i=0; i<dat->n_ind; i++) {
	gsl_vector_set(y, i, dat->pheno[i]);
    }


    gsl_matrix *X = gsl_matrix_alloc(dat->n_ind, dat->n_cov+2);

    for (size_t i=0; i<dat->n_ind; i++) {
	for (size_t j=0; j<dat->n_cov; j++) {
	    gsl_matrix_set(X, i, j, dat->covar[j][i]);
	}
    }

    gsl_matrix *cov = gsl_matrix_alloc(dat->n_cov+2, dat->n_cov+2);
    gsl_vector *beta = gsl_vector_alloc(dat->n_cov+2);
    double chisq = 0;
    double eff1, se1, eff2, se2;

    std::ofstream out(out_path);
    for (size_t i=0; i<dat->n_snp; i++) {
	for (size_t j=0; j<dat->n_ind; j++) {
	    gsl_matrix_set(X, j, dat->n_cov, dat->geno1[i][j]);
	    gsl_matrix_set(X, j, dat->n_cov+1, dat->geno2[i][j]);
	}
	gsl_multifit_linear(X, y, beta, cov, &chisq, work);
	eff1 = gsl_vector_get(beta, dat->n_cov);
	se1 = gsl_matrix_get(cov, dat->n_cov, dat->n_cov);
	eff2 = gsl_vector_get(beta, dat->n_cov+1);
	se2 = gsl_matrix_get(cov, dat->n_cov+1, dat->n_cov+1);
	out << eff1 << "\t" << se1 << "\t" << eff2 << "\t" << se2 << endl;
    }

    out.close();

    gsl_vector_free(y);
    gsl_vector_free(beta);
    gsl_matrix_free(X);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free(work);
}

