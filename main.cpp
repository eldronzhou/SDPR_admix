#include <iostream>
#include "mcmc.h"
#include <string.h>
#include "time.h"


using std::cout; using std::endl;
using std::string;

void print_use() {
    cout << "Usage: SDPR_admix -options" << endl << endl
    << "-iter (optional) number of iterations for MCMC. Default is 1000." << endl << endl
    << "-burn (optional) number of burn-in for MCMC. Default is 200." << endl << endl	
    << "-pheno (required) path to the phenotype file. The phenotype will be read from the 3rd column of the specified space- or tab-delimited file. There is no header and NA value can be included." << endl << endl
    << "-vcf (required) path to the phased genotype file." << endl << endl
    << "-msp (required) path to the directory containing Rfmix2 solved local ancestry files." << endl << endl
    << "-covar (optional) path to the covariate file. Covariates will be reading from the first column. There is no header for the covariate file." << endl << endl
    << "-out (required) path to the output file." << endl << endl
    << "-rho (required) Cross-ancestry genetic correlation." << endl << endl	
    << "-h print the options." << endl << endl;
}

int main(int argc, char *argv[]) {
    Dat dat;

    if (argc == 1) {
	print_use();
	return 0;
    }

    double rho = 0;

    std::string pheno_path, geno1_path, \
        geno2_path, vcf_path, msp_path, out_path, covar_path;

    int i = 1, iter = 1000, burn = 500;
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
	else if (strcmp(argv[i], "-h") == 0) {
	    print_use();
	    return 0;
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
}
