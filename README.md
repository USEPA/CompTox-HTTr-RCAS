# Integration of High-Throughput Transcriptomics and High-Throughput Screening Data Streams for MoA-Specific Chemical Prioritization
Authors: Rogers JD, Paul-Friedman K, Everett LJ

All necessary code for reproducing results in [Rogers et. al. *in prep* (broken link)]()

## Overview

This analysis makes use of several data sources to prioritize chemicals for specific mechanisms-of-action (MoAs): 

1. RefChemDB: semi-automated mining of multiple literature databases for molecular target associations for over 40,000 chemicals
2. High-throughput transcriptomics (HTTr): whole transcriptome sequencing of over 1,200 chemicals in multi-concentration format across HepaRG and U-2 OS cell lines
3. High-throughput Screening (HTS): targeted receptor-level assay readouts of thousands of chemicals in multi-concentration format from US EPA's ToxCast program

This code makes use of these sources in several major steps: assigning reference chemicals to molecular target clusters as a putative MoA, generating Reference Chemical-Associated Signatures (**RCAS**) as a means to profile transcriptomics data for activity related to MoA-associated reference chemicals, screening HTTr data for RCAS, and validating RCAS-based predictions with orthogonal assays from ToxCast.

Also included are scripts used to conduct all statistical analyses in the paper, as well as vignettes for generating reference chemicals, RCAS, and running the prioritization pipeline.  

## Installation

Clone this directory onto your desired path for access to all scripts:

```
git clone https://github.com/jessedrogers/NAM-integration/
```

Note that these scripts will generate two directories in the parent folder to your desired path (`../input` and `../output`), as expected by httrpathway. These directories and sub-directories will store intermediate and final transcriptomic data for use in signature concentration-response profiling, and necessary filespace should be dedicated for intermediate files of large HTTr screens (~8GB).

## Usage

See vignettes (TBD) for full examples of code usage. Each pipeline script contains a runtime function for performing steps in the function:

Step | Script | Runtime Function | Intended Function
|----|------|------------------|------------------|
| 1 | `pipeline_refchem_assignment.R` | `assignRefChems()` | Assign reference chemical clusters from RefChemDB |
| 2 | `pipeline_RCAS_generation.R` | `analyzeHTrANOVA()` | Generate RCAS from gene-level HTTr data |
| 3 | `pipeline_RCAS_profiling.R` | `profileRCAS()`| Perform concentration-response profiling of RCAS from HTTr data |
| 4 | `pipeline_HTS_selection.R` | `selectHTSEndpoints()` | Select orthogonal ToxCast endpoints from InvitroDB with high specificity |
| 5 | `pipeline_framework.R` | `runFramework()` | Perform target-based hazard assessment of HTTr/ToxCast-profiled chemicals using tiered framework | 

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
    - available from [USEPA httrpl repository](https://github.com/USEPA/httrpl_pilot)
- HTTr Concentration-Response Profiling Package (httrpathway)
    - available from [USEPA CompTox-httrpathway repository](https://github.com/USEPA/CompTox-httrpathway)
    - **Currently requires manual changes to NAMESPACE for successful installation**, recommended to clone repository into working directory

### External data requirements
- *RefChemDB*
    - Excel spreadsheet available from [Judson et. al. ALTEX 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6784312/) (Supplemental Table S12)
- *HTTr HepaRG/U-2 OS Datasets*
    - Full dataset may be provided upon request, and will be made available at a future date via [USEPA CompTox Chemicals Dashboard](https://comptox.epa.gov/dashboard)
- *ToxCast InvitroDB Endpoints*
    - MySQL database available from [InvitroDB v3.4](https://doi.org/10.23645/epacomptox.6062623)