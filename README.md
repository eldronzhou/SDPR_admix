# SDPR_admix
A statistical method to calculate PRS in admixed population. SDPR_admix is still in the early stage of testing. If you have encounbtered a problem, please contact geyu.zhou@yale.edu.

## Installation

To install SDPR, you need to first download the repo:

```
git clone https://github.com/eldronzhou/SDPR_admix.git
```

## Quick start

SDPR_admix can be run from the command line. To see the full list of options, please type

```bash
./SDPR_admix -h
```
Below are the required options.

-vcf (required): Path to the phased genotype file in the vcf format.
-msp (required): Path to path to the directory containing Rfmix2 solved local ancestry files.
-pheno (required): path to the phenotype file. The phenotype will be read from the 3rd column of the specified space- or tab-delimited file. There is no header and NA value can be included.
-covar (required): path to the covariate file. The format is same as the phenotype file (no header, FID and IID in first two columns, covariates in remaining columns).
-out (required): Path to the output file containing estimated effect sizes.
-rho (required): Trans-ethnic genetic correlation output by PopCorn between 0 and 1. Default is 0.8. 

## Running SDPR_admix

```bash
./SDPR_admix -vcf test/chr22_train.vcf -vcf chr22_train.vcf -msp test/chr22_train.msp.tsv -pheno test/train.pheno -covar test/covar.txt -out test/res.txt
```
The output has the following format:



## Derive the PRS
