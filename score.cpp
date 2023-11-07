#include "score.h"
#include <assert.h>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <unordered_map>
#include <regex>
#include <zlib.h>

using std::cout; using std::endl; using std::ifstream;
using std::string; using std::getline; using std::unordered_map;


void print_use() {
    cout << "Usage: score -options" << endl << endl
    << "-score (required) path to the weight file output by SDPR_admix. The format is." << endl << endl
    << "-vcf (required) path to the phased genotype file." << endl << endl
    << "-msp (required) path to the directory containing Rfmix2 solved local ancestry files." << endl << endl
    << "-out (required) path to the output file." << endl << endl
    << "-h print the options." << endl << endl;
}

void read_score(const std::string &score_path, unordered_map<string, Dat*> &dat_dict) {
    ifstream infile(score_path.c_str());

    if (!infile) {
	throw std::runtime_error("Error: cannot open "
		"the score file: " + score_path);
    }

    string id, A1, A2, line;
    double beta1, beta2;
    int chr;
    unsigned pos;

    int n_snp = 0;

    while (getline(infile, line, '\n')) {
	std::istringstream ss(line);
	ss >> chr >> pos >> id >> A1 >> A2 >> beta1 >> beta2;
	Dat *dat = new Dat;
	dat->chr = chr;
	dat->pos = pos;
	dat->A1 = A1;
	dat->A2 = A2;
	dat->beta1 = beta1;
	dat->beta2 = beta2;
	if (!dat_dict.insert(std::pair<string, \
		    Dat*>(id, dat)).second) {
	    throw std::runtime_error("Error: duplicate SNP found "
		    "in ref snpInfo file: " + id);
	}
	n_snp++;
    }
    cout << "Readed " << std::to_string(n_snp) << " SNPs from the score file." << endl;
}

void get_size_vcf(const std::string &vcf_path, PRS_dat* prs_dat) {
    size_t n_snp = 0;
    size_t n_ind = 0;
    string line, token;
    
    gzFile infile2 = gzopen(vcf_path.c_str(), "rb");
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
	    int idx = 0;
	    std::istringstream iss(line);
	    while (getline(iss, token, '\t')) {
		if (idx >= 9) {
		    prs_dat->id.push_back(token);
		}
		idx++;
	    }
	    n_ind = idx - 9;
	}
	else {
	    n_snp++;
	}
    }
    prs_dat->n_snp = n_snp;
    prs_dat->n_ind = n_ind;
    gzclose(infile2);
    cout << "In total " + std::to_string(n_snp) + " SNPs and " + \
	std::to_string(n_ind) + " to be readed." << endl;
}

void score(const std::string &vcf_path,  const std::string &msp_path, unordered_map<string, Dat*> &dat_dict, PRS_dat* prs_dat) {
    
    string line1, line2;
    string token1, token2;

    ifstream mspfile(msp_path.c_str());
    cout << "Reading RFmix msp file from: " + msp_path + "." << endl;
    
    
    gzFile infile2 = gzopen(vcf_path.c_str(), "rb");
    char buffer[4096];
    cout << "Reading VCF file from: " + vcf_path + "." << endl;

    // skip first two lines of msp file
    getline(mspfile, line1);
    getline(mspfile, line1);

    // skip the header of vcf file
    int n_ind = prs_dat->n_ind;
    while (gzgets(infile2, buffer, 4096)) {
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
	    continue;
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
    int n_snp = 0;

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
	string A1, A2;
	double beta1 = 1.0, beta2 = 1.0;
	unordered_map<string, Dat*>::iterator idx_dict;
	
	while ((chr_vcf == chr_msp && pos >= spos && pos < epos) || idx_snp == prs_dat->n_snp-1) {

	    if (idx_snp == prs_dat->n_snp-1) {
		assert(chr_vcf == chr_msp && pos == epos);
	    }

	    // reset the stream
	    std::istringstream iss2(line2);

	    size_t k = 0;

	    for (idx2=0; idx2<n_ind+9; idx2++) {
		getline(iss2, token2, '\t');

		// find the SNP
		if (idx2 == 2) {
		    idx_dict = dat_dict.find(token2);

		    if (idx_dict == dat_dict.end()) {
			idx2 = n_ind+9;
			break;
		    }
		}

		// check A1 and A2
		if (idx2 == 3) {
		    A1 = token2;
		}
		if (idx2 == 4) {
		    A2 = token2;
		}

		if (idx2 == 5) {
		    if (A1 == idx_dict->second->A1 && A2 == idx_dict->second->A2) {
			beta1 = idx_dict->second->beta1;
			beta2 = idx_dict->second->beta2;
			n_snp++;
		    }
		    else if (A2 == idx_dict->second->A1 && A1 == idx_dict->second->A2) {
			beta1 = -1.0 * idx_dict->second->beta1;
			beta2 = -1.0 * idx_dict->second->beta2;
			n_snp++;
		    }
		    else {
			idx2 = n_ind+9;
			break;
		    }
		}
		
		if (idx2 >= 9) {

		    // check phasing
		    if (token2[1] != '|') {
			cout << "Genotype must be phased." << endl;
			exit(1);
		    }
		    
		    // read the genotype
		    if (token2[0] == '.' || token2[3] == '.') {
			cout << "Missing genotype not supported yet." << endl;
			exit(1);
		    }
		    if (hap_lanc[2*(idx2-9)] == 0) {
			prs_dat->geno1[k] += beta1 * std::stod(&token2[0]); 
		    }
		    else {
			prs_dat->geno2[k] += beta2 * std::stod(&token2[0]);
		    }

		    if (hap_lanc[2*(idx2-9)+1] == 0) {
			prs_dat->geno1[k] += beta1 * std::stod(&token2[2]);
		    }
		    else  {
			prs_dat->geno2[k] += beta2 * std::stod(&token2[2]);
		    }
		    k++;
		}
	    }
	    
	    // read the next line of vcf and update pos
	    assert(idx2 == n_ind+9);
	    
	    if (gzgets(infile2, buffer, 4096)) {
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
    
    cout << "Calculated PRS for " << std::to_string(n_snp) << " SNPs." << endl;
    free(hap_lanc);
}


int main(int argc, char *argv[]) {
        Dat dat;

	std::string score_path, vcf_path, msp_path, out_path;

	if (argc == 1) {
	    print_use();
	    return 0;
	}

	int i = 1;
		 
	while (i < argc) {
	    if (strcmp(argv[i], "-score") == 0) {
		score_path = argv[i+1];
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
	    else if (strcmp(argv[i], "-out") == 0) {
		out_path = argv[i+1];
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

	unordered_map<string, Dat*> dat_dict;
	read_score(score_path.c_str(), dat_dict);

	PRS_dat prs_dat;

	string vcf_chr, msp_chr;
	
	if (vcf_path.find('#') < vcf_path.length()) {
	    for (size_t i=1; i<23; i++) {
		vcf_chr = std::regex_replace(vcf_path, std::regex("#"), std::to_string(i));
		msp_chr = std::regex_replace(msp_path, std::regex("#"), std::to_string(i));
		cout << "Working on Chr" << std::to_string(i) << "." << endl;
		get_size_vcf(vcf_chr.c_str(), &prs_dat);
		if (i == 1) {
		    prs_dat.geno1 = (double *) calloc(prs_dat.n_ind, sizeof(double));
		    prs_dat.geno2 = (double *) calloc(prs_dat.n_ind, sizeof(double));
		}
		score(vcf_chr.c_str(), msp_chr.c_str(), dat_dict, &prs_dat);
	    }
	}
	else {
	    cout << "Working on " << vcf_path << endl;
	    get_size_vcf(vcf_path.c_str(), &prs_dat);
	    prs_dat.geno1 = (double *) calloc(prs_dat.n_ind, sizeof(double));
	    prs_dat.geno2 = (double *) calloc(prs_dat.n_ind, sizeof(double));
	    score(vcf_path.c_str(), msp_path.c_str(), dat_dict, &prs_dat);
	}
	
	unordered_map<string, Dat*>::iterator it;
	for (it=dat_dict.begin(); it != dat_dict.end(); it++) {
	    delete(it->second);
	}

	std::ofstream out(out_path);
	for (size_t i=0; i<prs_dat.n_ind; i++) {
	    out << prs_dat.id[i] << "\t" << prs_dat.geno1[i] << "\t" << prs_dat.geno2[i] << endl; 
	}
	out.close();

	free(prs_dat.geno1);
	free(prs_dat.geno2);

	return 0;
}

