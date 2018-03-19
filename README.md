# BEEM
BEEM is an approach to infer models for microbial community dynamics based on metagenomic sequencing data (16S or shotgun-metagenomics). It is based on the commonly used [generalized Lotka-Volterra modelling](https://en.wikipedia.org/wiki/Generalized_Lotka–Volterra_equation) (gLVM) framework. BEEM uses an iterative EM-like algorithm to simultaneously infer scaling factors (microbial biomass) and model parameters (microbial growth rate and interaction terms) and can thus work directly with the relative abundance values that are obtained with metagenomic sequencing. A preprint describing this work will be posted on bioRxiv soon.

Note: BEEM stands for **B**iomass **E**stimation and model inference with an **E**xpectation **M**aximization-like algorithm. 

## Dependencies

BEEM was written in R (>3.3.1) and requires the following packages: 
 - foreach
 - doMC
 - lokern
 - pspline
 - monomvn

BEEM scripts can be loaded with the following command in R:
```r
source('path/to/this/repo/emFunctions.r')
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
 - `measurementID`: timepoint for the sample

### Sample data

We have provided several sample input files that were also analyzed in our manuscript.

#### Data from [Props et. al. (2016)](https://www.nature.com/articles/ismej2016117)

 - OTU count `table: isme_analysis/counts.sel.txt`
 - Metadata: `isme_analysis/metadata.sel.txt`

#### Data from [Gibbons et. al. (2017)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005364)

 - OTU count table: `time_series_analysis/{DA,DB,M3,F4}.counts.txt`
 - Metadata: `time_series_analysis/{DA,DB,M3,F4}.metadata.txt`

## Usage

### Basic Usage (R commands)

```r
## Load functions
source("emFunctions.r")
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

The commands for reproducing the analysis reportd in the manuscript are presented as two jupyter notebooks: (1) [notebook for Props et. al.](https://github.com/CSB5/BEEM/blob/master/isme.ipynb) and (2) [notebook for Gibbons et. al.](https://github.com/CSB5/BEEM/blob/master/time_series_meta.ipynb).

## Citation
C Li, L Tucker-Kellogg & N Nagarajan. (2018). System	Biology	Modeling	with Compositional Microbiome	Data	Reveals Personalized	Gut	Microbial	Dynamics	and	Keystone	Species. *Submitted*.

## Contact
Please direct any questions or feedback to Chenhao Li (lich@gis.a-star.edu.sg) and Niranjan Nagarajan (nagarajann@gis.a-star.edu.sg).
