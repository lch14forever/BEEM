# BEEM

<img src="logo.png" height="200" align="right" />

BEEM is an approach to infer models for microbial community dynamics based on metagenomic sequencing data (16S or shotgun-metagenomics). It is based on the commonly used [generalized Lotka-Volterra modelling](https://en.wikipedia.org/wiki/Generalized_Lotka–Volterra_equation) (gLVM) framework. BEEM uses an iterative EM algorithm to simultaneously infer scaling factors (microbial biomass) and model parameters (microbial growth rate and interaction terms) from **longitudinal** data and can thus work directly with the relative abundance values that are obtained with metagenomic sequencing.

Note: BEEM stands for **B**iomass **E**stimation and model inference with an **E**xpectation **M**aximization algorithm. We have now extended the BEEM framework to be able to work with cross-sectional data (BEEM-static, check out our R package [here](https://github.com/CSB5/BEEM-static)).

## Dependencies

BEEM was written in R (>=3.3.1) and requires the following packages: 
 - foreach
 - doMC: this currently only works on MacOS or LinuxOS
 - lokern
 - pspline
 - monomvn

You can install BEEM as an R package using devtools

```r
devtools::install_github('csb5/beem')
```

## Input data

The input files for BEEM should have the same format as described in the manual for [MDSINE](https://bitbucket.org/MDSINE/mdsine/). The following two files are required by BEEM:

### OTU table

This should be a tab-delimited text file whose first row has the sample IDs and the first column has the OTU IDs (or taxonomic annotations). Each row should then contain the relative abundance of one OTU across all samples and each column should contain the relative abundances of all OTUs in that sample. 

### Metadata

The metadata file should be a tab-delimited text file with the following columns:
```
sampleID    isIncluded    subjectID    measurementID
```
 - `sampleID`: sample IDs matching the first row of the OTU table
 - `isIncluded`: whether the sample should be included in the analysis (1-include, 0-exclude)
 - `subjectID`: indicator for which biological replicate the sample belongs to
 - `measurementID`: time in standardized units from the start of the experiment

### Sample data

We have provided several sample input files that were also analyzed in our manuscript.

#### Data from [Props et. al. (2016)](https://www.nature.com/articles/ismej2016117)

 - OTU count table: `vignettes/props_et_al_analysis/counts.sel.txt`
 - Metadata: `vignettes/props_et_al_analysis/metadata.sel.txt`

#### Data from [Gibbons et. al. (2017)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005364)

 - OTU count table: `vignettes/gibbons_et_al_analysis/{DA,DB,M3,F4}.counts.txt`
 - Metadata: `vignettes/gibbons_et_al_analysis/{DA,DB,M3,F4}.metadata.txt`

## Usage

### Basic Usage (R commands)

```r
## Load functions
library(beem)
## Read inputs
counts <- read.table('counts.txt', head=F, row.names=1)
metadata <- read.table('metadata.txt', head=T)
## Run BEEM
res <- EM(dat=input, meta=metadata)
## Estimate parameters
biomass <- biomassFromEM(res)
write.table(biomass, 'biomass.txt', col.names=F, row.names=F, quote=F)
gLVparameters <- paramFromEM(res, counts, metadata)
write.table(gLVparameters, 'gLVparameters.txt', col.names=T, row.names=F, sep='\t' , quote=F)
```
### Output format

BEEM estimated parameters is an R `data.frame` (a table) with the following columns in order:
 
 - `parameter_type`: `growth_rate` or `interaction`
 - `source_taxon`: source taxon for interaction (`NA` if `parameter_type` is `growth_rate`)
 - `target_taxon`: target taxon for interaction or growth rate
 - `value`: parameter value 
 - `significance`: confidence level of the inferred interaction (only meaningful for interactions)
 
### Analyses in the manuscript

The commands for reproducing the analysis reportd in the manuscript are presented as jupyter notebooks: (1) [notebook on a demo of the gLVM simulation](https://github.com/CSB5/BEEM/blob/master/vignettes/simulation.ipynb), (2) [notebook for Props et. al.](https://github.com/CSB5/BEEM/blob/master/vignettes/props_et_al.ipynb) and (3) [notebook for Gibbons et. al.](https://github.com/CSB5/BEEM/blob/master/vignettes/gibbons_et_al.ipynb).

## Citation
C Li, K R Chng, J S Kwah, T V Av-Shalom, L Tucker-Kellogg & N Nagarajan. (2019). An expectation-maximization algorithm enables accurate ecological modeling using longitudinal metagenome sequencing data. [*Microbiome*](https://rdcu.be/bPl3T).

## Contact
Please direct any questions or feedback to Chenhao Li (cli40@mgh.harvard.edu) and Niranjan Nagarajan (nagarajann@gis.a-star.edu.sg).
