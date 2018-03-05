# BEEM
BEEM stands for *B*iomass *E*stimation and model inference with an *E*xpectation *M*aximization like algorithm. BEEM implemented an iterative EM-like algorithm to solve the total biomass (microbiome load) from high-throughput environmental sequencing data (16S targeted or shotgun-metagenomics) based on the generalized Lotka-Volterra Model (gLVM). BEEM also provided algorithms to estimate gLVM parameters (microbial growth rates and pair-wise interactions).

## Denpendencies

BEEM was written in R (>3.3.1) and requires the following packages: 
 - foreach
 - doMC
 - lokern
 - pspline
 - monomvn

BEEM scripts can be loaded with the following command in R:
```{R}
source('path/to/this/repo/emFunctions.r')
```
## Input data

The input files for BEEM should have the same format described in the manual of [MDSINE](https://bitbucket.org/MDSINE/mdsine/). The following two files are required by BEEM:

### OTU table

This should be an abundance table file (tab-delimited text file), whose first row has the sample IDs and the and first column has the OTU IDs (or taxonomic annotations). Each row is the abundance of one OTU across all samples and each column contains the abundances of all OTUs in that sample. 

### Metadata

The metadata should be a table (tab-delimited text file) with the following columns:
```
sampleID    isIncluded    subjectID    measurementID
```
 - `sampleID`: the sample IDs matching the first row of the OTU table
 - `isIncluded`: wether the sample should be included in the analysis (1-include, 0-exclude)
 - `subjectID`: indicator for which biological replicate the sample belongs to
 - `measurementID`: time stamp for the sample

### Sample data

We provided several sample input files that were also analyzed in our manuscript.

#### Data from [Props et. al. (2016)](https://www.nature.com/articles/ismej2016117)

 - OTU count `table: isme_analysis/counts.sel.txt`
 - Metadata: `isme_analysis/metadata.sel.txt`

#### Data from [Gibbons et. al. (2017)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005364)

 - OTU count table: `time_series_analysis/{DA,DB,M3,F4}.counts.txt`
 - Metadata: `time_series_analysis/{DA,DB,M3,F4}.metadata.txt`

## Usage

### Basic Usage (R commands)

```{R}
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
 - `target_taxon`: target taxon for interaction 
 - `value`: parmater value 
 - `significance`: confidence level of the inferred interaction (only meaningful for interactions)
 
### Analyses in the manuscript

The commands for reproducing the analysis in the manuscript were presented as two jupyter notebooks: (1) [notebook for Props et. al.](https://github.com/CSB5/BEEM/blob/master/isme.ipynb) and (2) [notebook for Gibbons et. al.](https://github.com/CSB5/BEEM/blob/master/time_series_meta.ipynb).
 
