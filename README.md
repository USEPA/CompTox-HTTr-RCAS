# Integrating Transcriptomic and Targeted New Approach Methodologies into a Tiered Framework for Chemical Bioactivity Screening
Authors: Rogers JD, Bundy JL, Harrill JA, Judson RJ, Paul Friedman K, Everett LJ

All necessary code for reproducing analyses in [Rogers et. al. *in prep* (Link TBD)]()

## Overview

This analysis makes use of several data sources to prioritize chemicals for specific mechanisms-of-action (MoAs): 

1. RefChemDB: semi-automated mining of source databases for molecular target associations for over 40,000 chemicals
2. High-throughput transcriptomics (HTTr): whole transcriptome profiling of 1,201 chemicals in multi-concentration format across HepaRG and U-2 OS cell lines
3. High-throughput Screening (HTS): targeted receptor-level assay readouts of thousands of chemicals in multi-concentration format from US EPA's ToxCast program

This code makes use of these sources in several major steps: assigning reference chemicals to molecular target clusters as a putative mechanism-of-action (MoA), generating Reference Chemical-Associated Signatures (**RCAS**) as a means to profile transcriptomics data for activity related to MoA-associated reference chemicals, screening HTTr data for RCAS, and validating RCAS-based predictions with orthogonal assays from ToxCast.

Also included are scripts used to conduct all statistical analyses in the paper, as well as vignettes for generating reference chemicals, RCAS, and running the prioritization pipeline.  

## Installation

Clone this directory onto your desired path for access to all scripts:

```
git clone https://github.com/jessedrogers/NAM-integration/
```

Note that these scripts will generate two directories in the parent folder to your desired path (`../input` and `../output`), as expected by httrpathway. These directories and sub-directories will store intermediate and final transcriptomic data for use in signature concentration-response profiling, and ~8GB storage should be dedicated for intermediate files of large HTTr screens.

## Usage

See vignettes sub-directory for full examples of code usage. Each R script contains a runtime function for performing steps in the pipeline:

| Order | Script Name | Runtime Function | Intended Usage | Related Vignette |
|-------|-------------|------------------|----------------|------------------|
| 1 | `pipeline_refchem_assignment.R` | `assignRefChems()` | Assign reference chemical clusters from RefChemDB | `analysis_RCAS.Rmd` |
| 2 | `pipeline_RCAS_generation.R` | `analyzeHTrANOVA()` | Generate RCAS from gene-level HTTr data | `analysis_RCAS.Rmd` |
| 3 | `pipeline_RCAS_profiling.R` | `profileRCAS()`| Perform concentration-response profiling of RCAS from HTTr data | `analysis_RCAS.Rmd` |
| 4 | `pipeline_HTS_selection.R` | `selectHTSEndpoints()` | Select orthogonal ToxCast endpoints from InvitroDB with high specificity | `analysis_framework.Rmd` |
| 5 | `pipeline_framework.R` | `runFramework()` | Perform target-based hazard assessment of HTTr/ToxCast-profiled chemicals using tiered framework | `analysis_framework.Rmd` |

## Dependencies

All code was written and tested using R 3.6.0, and should run using later versions.

| Use Case | Package(s)
|---------|---------|
| General Data <br> Manipulation | tidyr 1.1.4 <br> dplyr 1.0.7 <br> tibble 3.1.6 <br> magrittr 2.0.1 <br> tidyselect 1.1.1 <br> broom 0.7.11 <br> openxlsx 4.2.5 <br> reshape2 1.4.4 <br> stringi 1.7.6 <br> stringr 1.4.0 |
Plotting and <br> Visualization | ggplot2 3.3.5 <br> ggrepel 0.9.1 <br> RColorBrewer 1.1-2 <br> ComplexHeatmap 2.2.0 <br> circlize 0.4.13 |
Statistical <br> Analyses | foreach 1.5.1 <br> doParallel 1.0.16 <br> philentropy 0.5.0 <br> tcplfit2 0.1.3
ToxCast Data <br> Download | tcpl 2.0.2 <br> RMySQL 0.10.22 <br> DBI 1.1.2

### Additional code dependencies

- HTTr Pipeline Code (httrpl)
    - available from [USEPA httrpl repository](https://github.com/USEPA/CompTox-httrpl)
- HTTr Concentration-Response Profiling Package (httrpathway)
    - available from [USEPA CompTox-httrpathway repository](https://github.com/USEPA/CompTox-httrpathway)
    - **Currently requires manual changes to NAMESPACE for successful installation**, recommended to clone repository into working directory

### External data requirements
- *RefChemDB*
    - Excel spreadsheet available from [Judson et. al. ALTEX 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6784312/) (Supplemental Table S12)
- *HTTr HepaRG/U-2 OS Datasets*
    - Raw and count data available on NCBI Gene Expression Omnibus: [U-2 OS (Accession TBD)](), [HepaRG (Accession TBD)]()
    - Guidance for differential expression analysis via DESeq2 is available via the httrpl package
- *ToxCast HTS Endpoints*
    - MySQL database available from [InvitroDB v3.4](https://doi.org/10.23645/epacomptox.6062623)