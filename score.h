#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>

typedef struct {
    int chr;
    unsigned pos;
    std::string A1;
    std::string A2;
    double beta1;
    double beta2;
} Dat;

typedef struct {
    double *geno1;
    double *geno2;
    size_t n_snp;
    size_t n_ind;
    std::vector<std::string> id;
} PRS_dat;

void score(const std::string &vcf_path,  const std::string &msp_path, std::unordered_map<std::string, Dat*> &dat_dict);

