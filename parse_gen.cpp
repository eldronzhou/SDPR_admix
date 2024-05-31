#include "parse_gen.h"
#include <assert.h>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <zlib.h>
#include <cstring>

using std::cout; using std::endl; using std::ifstream;
using std::string; using std::getline;

// get n_snp and n_ind
void get_size_vcf(const string &pheno_path, const string &geno_path, Dat *dat) {
    size_t n_ind = 0;
    size_t n_invalid = 0;
    size_t n_snp = 0;

    string line;
    ifstream infile1(pheno_path.c_str());
    string id;
    string y;
    size_t i = 0;
    while (infile1 >> id >> id >> y) {
	try {
	    std::stod(y); 
	    dat->ind_idx.push_back(i);
	    n_ind++;
	}
	catch (std::invalid_argument&) {
	    n_invalid++;
	}
	i++;
    }
    dat->n_ind = n_ind;
    cout << "Warning: " + std::to_string(n_invalid) + \
	" individuals with invalid phenotypes." << endl;

    gzFile infile2 = gzopen(geno_path.c_str(), "rb");

    char buffer[4096];

    while (gzgets(infile2, buffer, 4096)) {
	// Remove the newline character if it exists
        size_t length = strlen(buffer);
        if (length > 0 && buffer[length - 1] == '\n') {
            buffer[length - 1] = '\0';
        }
	
	line = buffer;
	if (line.find("##") == 0) {
	    continue;
	}
	else if (line.find("#") == 0) {
	    continue;
	}
	else {
	    n_snp++;
	}
    }

    gzclose(infile2);

    /*ifstream infile2(geno_path.c_str());
    while(getline(infile2, line)) {
	if (line.find("##") == 0) {
	    continue;
	}
	else if (line.find("#") == 0) {
	    continue;
	}
	else {
	    n_snp++;	
	}
    }*/

    dat->n_snp = n_snp;
    cout << "In total " + std::to_string(n_snp) + " SNPs and " \
	+ std::to_string(n_ind) + " individuals to be readed." << endl;
}

void read_pheno(const std::string &pheno_path, Dat *dat) {
    ifstream infile(pheno_path.c_str());

    cout << "Reading phenotype file from: " + pheno_path + "." << endl;

    string id;
    string y;
    double *pheno = (double *) malloc(dat->n_ind*sizeof(double));

    size_t i = 0, idx = 0;
    while (infile >> id >> id >> y) {
	if (i == dat->ind_idx[idx]) {
	    pheno[idx] = stod(y);
	    idx++;
	}
	i++;
    }
    dat->pheno = pheno;

    cout << "Readed phenotype from " + std::to_string(idx) + " individuals." << endl;
}

void read_cov(const std::string &cov_path,  Dat *dat) {
    ifstream infile(cov_path.c_str());

    cout << "Reading covariate file from: " + cov_path + "." << endl;

    size_t n_cov = 0, i = 0, idx = 0;

    string line, token;

    while (getline(infile, line)) {
	n_cov = 0;
	std::istringstream iss(line);
	if (i == 0) {
	    while (getline(iss, token, '\t')) {   
		n_cov++;
	    }
	    cout << "Reading " + std::to_string(n_cov) + " covariates." << endl;
	    dat->n_cov = n_cov;
	    dat->covar = (double **) malloc(n_cov*sizeof(double *));
	    for (size_t k=0; k<n_cov; k++) {
		dat->covar[k] = (double *) malloc(dat->n_ind*sizeof(double));
	    }
	    n_cov = 0;
	    iss.clear();
	    iss.str(line);
	}
	if (i == dat->ind_idx[idx]) {
	    while (getline(iss, token, '\t')) {
		dat->covar[n_cov][idx] = stod(token);
		n_cov++;
	    }
	    idx++;
	}
	i++;
    }
    infile.close();
}

void read_lanc(const std::string &vcf_path,  const std::string &msp_path, Dat *dat) {
    
    string line1, line2;
    string token1, token2;

    ifstream mspfile(msp_path.c_str());
    cout << "Reading RFmix msp file from: " + msp_path + "." << endl;
    
    //ifstream infile(vcf_path.c_str());
    //cout << "Reading VCF file from: " + vcf_path + "." << endl;
    gzFile infile = gzopen(vcf_path.c_str(), "rb");
    char buffer[4096];

    // skip first two lines of msp file
    getline(mspfile, line1);
    getline(mspfile, line1);

    // skip the header of vcf file
    int n_ind = 0;
    while (gzgets(infile, buffer, 4096)) {
	// Remove the newline character if it exists
        size_t length = strlen(buffer);
        if (length > 0 && buffer[length - 1] == '\n') {
            buffer[length - 1] = '\0';
        }
	line2 = buffer;

	if (line2.find("##") == 0) {
	    continue;
	}
	else if (line2.find("#") == 0) {
	    int idx_2 = 0;
	    std::istringstream iss2(line2);
	    while (getline(iss2, token2, '\t')) {
		idx_2++;
	    }
	    n_ind = idx_2 - 9;
	}
	else {
	    break;
	}
    }

    int *hap_lanc = (int *) malloc(2*n_ind*sizeof(int));
    unsigned spos = 0, epos = 0, pos = 0;
    int chr_vcf = 1, chr_msp = 1;
    
    // read the pos from the first line of vcf file
    std::istringstream iss2(line2);
    int idx2 = 0;
    for (; idx2<9; idx2++) {
	getline(iss2, token2, '\t');
	if (idx2 == 0) {
	    chr_vcf = std::stoi(token2);
	}
	if (idx2 == 1) {
	    pos = std::stoul(token2);
	}
    }

    size_t idx_snp = 0;
    dat->geno1 = (double **) malloc(dat->n_snp*sizeof(double*));
    dat->geno2 = (double **) malloc(dat->n_snp*sizeof(double*));
    dat->n_anc1 = (double *) calloc(dat->n_snp, sizeof(double));
    dat->n_anc2 = (double *) calloc(dat->n_snp, sizeof(double));
    for (size_t i=0; i<dat->n_snp; i++) {
	dat->geno1[i] = (double *) calloc(dat->n_ind, sizeof(double));
	dat->geno2[i] = (double *) calloc(dat->n_ind, sizeof(double));
	if (!dat->geno1[i] || !dat->geno2[i]) {
	    cout << "Error: memory allocation failed for" + \
		std::to_string(i) + " th SNP.";
	    exit(EXIT_FAILURE);
	}
    }

    while (getline(mspfile, line1)) {
	// read msp file
	std::istringstream iss1(line1);
	for (int idx1=0; idx1<2*n_ind+6; idx1++) {
	    getline(iss1, token1, '\t'); 
	    if (idx1 == 0) {
		chr_msp = std::stoi(token1);
	    }
	    else if (idx1 == 1) {
		spos = std::stoul(token1);
	    }
	    else if (idx1 == 2) {
		epos = std::stoul(token1);
	    }
	    else if (idx1 >= 6) {
		hap_lanc[idx1-6] = std::stoi(token1);
		if (hap_lanc[idx1-6] != 0 && hap_lanc[idx1-6] != 1) {
		    cout << "RFmix field must be either 0 or 1." << endl;
		    return;
		}
	    }
	}

	if ((chr_vcf != chr_msp) || (pos != spos)) {
	    cout << "Inconsistent starting position: chr_vcf: " + std::to_string(chr_vcf) + \
		" chr_msp: " + std::to_string(chr_msp) + " pos: " + \
		std::to_string(pos) + " spos: " + std::to_string(spos) << endl;
	    exit(EXIT_FAILURE);
	} 
	     
	// read vcf file
	while ((chr_vcf == chr_msp && pos >= spos && pos < epos) || idx_snp == dat->n_snp-1) {

	    if (idx_snp == dat->n_snp-1) {
		assert(chr_vcf == chr_msp && pos == epos);
	    }

	    // reset the stream
	    std::istringstream iss2(line2);

	    size_t k = 0;
	    for (idx2=0; idx2<n_ind+9; idx2++) {
		getline(iss2, token2, '\t');

		if (idx2 == 0) {
		    dat->chr.push_back(token2);
		}

		if (idx2 == 1) {
		    dat->pos.push_back(token2);
		}

		if (idx2 == 2) {
		    dat->id.push_back(token2);
		}

		if (idx2 == 3) {
		    dat->ref.push_back(token2);
		}

		if (idx2 == 4) {
		    dat->alt.push_back(token2);
		}

		if (idx2 >= 9) {

		    // individuals kept in analysis
		    if (idx2-9 != dat->ind_idx[k]) {
			continue;
		    }

		    // check phasing
		    if (token2[1] != '|') {
			cout << "Genotype must be phased." << endl;
			return;
		    }
		    
		    // read the genotype
		    if (token2[0] == '.' || token2[3] == '.') {
			cout << "Missing genotype not supported yet." << endl;
			return;
		    }
		    if (hap_lanc[2*(idx2-9)] == 0) {
			dat->geno1[idx_snp][k] += std::stod(&token2[0]);
			dat->n_anc1[idx_snp]++;
		    }
		    else {
			dat->geno2[idx_snp][k] += std::stod(&token2[0]);
			dat->n_anc2[idx_snp]++;
		    }

		    if (hap_lanc[2*(idx2-9)+1] == 0) {
			dat->geno1[idx_snp][k] += std::stod(&token2[2]);
			dat->n_anc1[idx_snp]++;
		    }
		    else  {
			dat->geno2[idx_snp][k] += std::stod(&token2[2]);
			dat->n_anc2[idx_snp]++;
		    }
		    k++;
		}
	    }
	    
	    // read the next line of vcf and update pos
	    assert(idx2 == n_ind+9);
	    
	    if (gzgets(infile, buffer, 4096)) {
		// Remove the newline character if it exists
		size_t length = strlen(buffer);
		if (length > 0 && buffer[length - 1] == '\n') {
		    buffer[length - 1] = '\0';
		}

		line2 = buffer;
		iss2.clear();
		iss2.str(line2);
		for (idx2=0; idx2<9; idx2++) {
		    getline(iss2, token2, '\t');
		    if (idx2 == 0) {
			chr_vcf = std::stoi(token2);
		    }
		    if (idx2 == 1) {
			pos = std::stoul(token2);
		    }
		}
		idx_snp++;
	    }
	    else {
		break;
	    }
	}	
    }
    
    cout << "Readed " << std::to_string(idx_snp+1) << \
	" SNPs from " << std::to_string(dat->n_ind) << " individuals." << endl;
    free(hap_lanc);
    gzclose(infile);

    // write to the output
    /*std::ofstream out1("./X_afr.txt");
    std::ofstream out2("./X_eur.txt");
    for (size_t i=0; i<dat->n_snp; i++) {
	for (size_t j=0; j<dat->n_ind; j++) {
	    if (j != dat->n_ind-1) {
		out1 << dat->geno1[i][j] << " ";
		out2 << dat->geno2[i][j] << " ";
	    }
	    else {
		out1 << dat->geno1[i][j] << endl;
		out2 << dat->geno2[i][j] << endl;
	    }
	}
    }
    out1.close();
    out2.close();*/
}

void prod_geno(Dat *dat) {
    cout << "Precalculating genotype related products." << endl;

    dat->geno1_sq = (double *) calloc(dat->n_snp, sizeof(double));
    dat->geno2_sq = (double *) calloc(dat->n_snp, sizeof(double));
    dat->geno12_prod = (double *) calloc(dat->n_snp, sizeof(double));

    if (!dat->geno1_sq || !dat->geno2_sq || !dat->geno12_prod) {
	cout << "Error: memory allocation failed for genotype related products." << endl;
	exit(EXIT_FAILURE);
    }

    for (size_t i=0; i<dat->n_snp; i++) {
	for (size_t j=0; j<dat->n_ind; j++) {
	    dat->geno1_sq[i] += dat->geno1[i][j]*dat->geno1[i][j];	    
	}
    }

    for (size_t i=0; i<dat->n_snp; i++) {
	for (size_t j=0; j<dat->n_ind; j++) {
	    dat->geno2_sq[i] += dat->geno2[i][j]*dat->geno2[i][j];
	}
    }

    for (size_t i=0; i<dat->n_snp; i++) {
	for (size_t j=0; j<dat->n_ind; j++) {
	    dat->geno12_prod[i] += dat->geno1[i][j]*dat->geno2[i][j];
	}
    }

    cout << "Finished calculating genotype related products." << endl;
}

void check_maf(Dat *dat, double maf) {
    size_t n_bad = 0;
    std::vector<size_t> ok_idx;
    std::vector<std::string> chr;
    std::vector<std::string> id;
    std::vector<std::string> pos;
    std::vector<std::string> ref;
    std::vector<std::string> alt;

    dat->maf1 = (double *) malloc(dat->n_snp*sizeof(double));
    dat->maf2 = (double *) malloc(dat->n_snp*sizeof(double));
    
    for (size_t i=0; i<dat->n_snp; i++) {
	double mean1 = 0, mean2 = 0;
	for (size_t j=0; j<dat->n_ind; j++) {
	    mean1 += dat->geno1[i][j]; 
	    mean2 += dat->geno2[i][j]; 	
	}
	mean1 /= (dat->n_anc1[i]); 	
	mean2 /= (dat->n_anc2[i]);
	if (mean1 > .5) {
	    mean1 = 1 - mean1;
	}
	if (mean2 > .5) {
	    mean2 = 1 - mean2;
	}
	dat->maf1[i] = mean1;
	dat->maf2[i] = mean2;

	if (mean1 < maf || mean1 > 1-maf || mean2 < maf || mean2 > 1-maf) {
	    n_bad++;
	    continue;
	}
	ok_idx.push_back(i);
    }
    
    size_t n_snp = dat->n_snp - n_bad;

    double **geno1 = (double **) malloc(n_snp*sizeof(double*));
    double **geno2 = (double **) malloc(n_snp*sizeof(double*));

    size_t k = 0;
    for (size_t i=0; i<dat->n_snp; i++) {
	if (i == ok_idx[k]) {
	    geno1[k] = dat->geno1[i];
	    geno2[k] = dat->geno2[i];
	    chr.push_back(dat->chr[i]);
	    id.push_back(dat->id[i]);
	    pos.push_back(dat->pos[i]);
	    ref.push_back(dat->ref[i]);
	    alt.push_back(dat->alt[i]);
	    k++;
	}
	else {
	    free(dat->geno1[i]);
	    free(dat->geno2[i]);
	}
    }
    free(dat->geno1);
    free(dat->geno2);
    dat->geno1 = geno1;
    dat->geno2 = geno2;
    dat->n_snp = n_snp;
    dat->chr = chr;
    dat->id = id;
    dat->pos = pos;
    dat->ref = ref;
    dat->alt = alt;
    cout << "Excluded " + std::to_string(n_bad) << " SNPs due to MAF filter." << endl;
}




