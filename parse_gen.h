#include <iostream>
#include <fstream>
#include <vector>

#ifndef PARSE_GEN_H
#define PARSE_GEN_H

typedef struct {
    double **geno1;
    double **geno2;
    double *pheno;
    double *y;
    std::vector<size_t> ind_idx;
    size_t n_snp;
    size_t n_ind;
    size_t n_cov = 0;
    double *geno1_sq;
    double *geno2_sq;
    double *geno12_prod;
    double **covar;
    double **proj;
    std::vector<std::string> chr;
    std::vector<std::string> id;
    std::vector<std::string> pos;
    std::vector<std::string> ref;
    std::vector<std::string> alt;
    double *maf1;
    double *maf2;
    double *n_anc1;
    double *n_anc2;
} Dat;

#endif

void get_size_vcf(const std::string &pheno_path, const std::string &geno_path, Dat *dat);

void read_geno(const std::string &geno_path, Dat *dat, size_t pop);

void read_pheno(const std::string &pheno_path, Dat *dat);

void read_cov(const std::string &cov_path,  Dat *dat);

void read_lanc(const std::string &vcf_path,  const std::string &msp_path, Dat *dat);

void check_maf(Dat *dat, double maf);

void prod_geno(Dat *dat);
