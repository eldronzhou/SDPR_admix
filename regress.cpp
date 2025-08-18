#include "regress.h"
#include "parse_gen.h"
#include <cstring>
#include "gsl/gsl_randist.h"
#include "gsl/gsl_statistics_double.h"
#include "gsl/gsl_blas.h"
using std::cout; using std::endl;

using std::cout; using std::endl; 

void linear(Dat *dat, std::string out_path, int thread) {
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
	out << dat->id[i] << "\t" << dat->ref[i] << "\t" << dat->alt[i] << "\t" << eff1 << "\t" << se1 << "\t" << eff2 << "\t" << se2 << endl;
    }

    out.close();

    gsl_vector_free(y);
    gsl_vector_free(beta);
    gsl_matrix_free(X);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free(work);
}

/*int main(int argc, char *argv[]) {
    Dat dat;

    std::string pheno_path, vcf_path, msp_path, out_path, covar_path;

    int anc = 0;

    int i = 1, thread = 1;
    while (i < argc) {
        if (strcmp(argv[i], "-pheno") == 0) {
            pheno_path = argv[i+1];
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
	else if (strcmp(argv[i], "-anc") == 0) {
            anc = std::stoi(argv[i+1]);
            i += 1;
        }
        else if (strcmp(argv[i], "-out") == 0) {
            out_path = argv[i+1];
            i += 2;
        }
        else if (strcmp(argv[i], "-covar") == 0) {
            covar_path = argv[i+1];
            i += 2;
        }
	else if (strcmp(argv[i], "-thread") == 0) {
            thread = std::stoi(argv[i+1]);
            i += 2;
        }
    }

    get_size_vcf(pheno_path.c_str(), vcf_path.c_str(), &dat);

    read_lanc(vcf_path.c_str(), msp_path.c_str(), anc, &dat);

    read_pheno(pheno_path.c_str(), &dat);

    if (!covar_path.empty()) {
        read_cov(covar_path.c_str(), &dat);
    }
    linear(&dat, out_path.c_str(), thread);
}*/

