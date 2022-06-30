# Integration of High-Throughput Transcriptomics and High-Throughput Screening Data Streams for MoA-Specific Chemical Prioritization
Authors: Rogers JD, Paul-Friedman K, Everett LJ

All necessary code for reproducing results in [Rogers et. al. *in prep* (broken link)]()

## Usage and Examples

## Software Requirements

- R v3.6.0
    - tidyr v1.1.4
    - dplyr v1.0.7
    - broom v0.7.11
    - ggplot2 v3.3.5
    - ggrepel v0.9.1
    - RColorBrewer v1.1-2
    - ComplexHeatmap v2.2.0
    - circlize v0.4.13
    - readxl v1.3.1
    - data.table v1.14.2
    - foreach v1.5.1
    - doParallel v1.0.16
    - RMySQL v0.10.22
    - DBI v1.1.2
    - tcpl v2.0.2
    - tcplfit2: available from [USEPA CompTox-ToxCast-tcplFit2 repository](https://github.com/USEPA/CompTox-ToxCast-tcplFit2)
- HTTr Pipeline Code (httrpl)
    - available from [USEPA httrpl repository](https://github.com/USEPA/httrpl_pilot)
- HTTr Concentration-Response Profiling Package (httrpathway)
    - available from [USEPA CompTox-httrpathway repository](https://github.com/USEPA/CompTox-httrpathway)
    - **Currently requires manual changes to NAMESPACE for successful package installation**

## External Data Requirements
- *RefChemDB*
    - Excel spreadsheet available from [Judson et. al. ALTEX 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6784312/) (Supplemental Table S12)
- *HTTr HepaRG/U-2 OS Datasets*
    - Full dataset may be provided upon request, and will be made available at a future date via [CompTox Chemicals Dashboard](https://comptox.epa.gov/dashboard)
- *ToxCast InvitroDB Endpoints*
    - MySQL database available from [InvitroDB v3.4](https://doi.org/10.23645/epacomptox.6062623)