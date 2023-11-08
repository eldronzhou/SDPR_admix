# SDPR_admix
A statistical method to calculate PRS in admixed population. SDPR_admix is still in the early stage of testing. If you have encounbtered a problem, please contact geyu.zhou@yale.edu.

## Installation

To install SDPR, you need to first download the repo:

```
git clone https://github.com/eldronzhou/SDPR_admix.git
```

You can then compile SDPR by running `make`. If there is run time error that the shared library "libgsl.so.0" not found, you can fix it by typing

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:SDPR_admix_dir/gsl/lib
```


## Quick start

SDPR_admix can be run from the command line. To see the full list of options, please type

```bash
./SDPR_admix -h
```
Below are the required options.

- vcf (required): Path to the phased genotype file in the gzipped vcf format.
- msp (required): Path to path to the directory containing RFMix2 solved local ancestry files.
- pheno (required): path to the phenotype file. The phenotype will be read from the 3rd column of the specified space- or tab-delimited file. There is no header and NA value can be included.
- covar (required): path to the covariate file. The format is same as the phenotype file (no header, FID and IID in first two columns, covariates in remaining columns).
- out (required): Path to the output file containing estimated effect sizes.
- rho (required): Trans-ethnic genetic correlation output by PopCorn between 0 and 1. Default is 0.8. 

## Running SDPR_admix

```bash
./SDPR_admix -vcf test/chr22_train.vcf.gz -msp test/chr22_train.msp.tsv -pheno test/train.pheno -covar test/covar.txt -out test/res.txt
```
The output file has the following format:

```
1       833068  rs12562034      G       A       0       0
1       843942  rs4040617       A       G       0       0
...
```
where the columns are chromsome, position, variant ID, effect allele, non-effect allele, effect sizes corrsponding to population 0 in RFMix2 file and effect sizes corrsponding to population 1 in the RFMix2 file.

## Derive the PRS

The following command can be used to calculate ancestry-aware PRS for the admixed population.

```bash
./score -vcf test/chr22_train.vcf.gz -msp test/chr22_train.msp.tsv -score test/res.txt -out test/test.profile
```
If you have bgzipped vcf and RFMix2 local ancestry files for chr1-22 with the name `test/chr[1-22]_train.vcf.gz`, then you can use the following command for iterative calculation over all chromsomes:

```bash
./score -vcf test/chr#_train.vcf.gz -msp test/chr#_train.msp.tsv -score test/res.txt -out test/test.profile
```
