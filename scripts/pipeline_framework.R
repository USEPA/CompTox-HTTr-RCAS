# NAM-integration-pilot workflow 5
# Tiered assessment of chemical screening data for selective target perturbagens
library(dplyr)
library(tidyr)
library(openxlsx)
source("scripts/pipeline_RCAS_profiling.R")
source("scripts/pipeline_HTS_selection.R")

loadRCASprofiles <- function(filename_cr) {
    #' Import of RCAS concentration-response profiling results
    #' 
    #' @param filename_cr character | location of file containing concentration-
    #'  response profiling results as output by profileRCAS(). Note that files
    #'  must be located in the directory output by httrpathway
    #'  ("../output/singature_conc_resp_summary/")
    #' @return data.frame of curve fit estimates for each RCAS
    #' @example
    #' @export
    dir_cr <- "../output/signature_conc_resp_summary/"
    filepath_cr <- paste0(dir_cr, filename_cr)
    # check for file compatibility
    if (!file.exists(filepath_cr)) {
        stop("Error: Problem with `filepath_cr`. \n `filepath_cr` not found")
    } else if (tools::file_ext(filepath_cr) != "RData") {
        stop("Error: `filepath_cr` must be an RData file")
    }
    # import + harmonize object name
    load(filepath_cr)
    cr_rcas <- get(ls()[!ls() %in% c("filename_cr", "dir_cr", "filepath_cr")])
    return(cr_rcas)
}

loadCatalogProfiles <- function(filepath_catalog) {
    #' Import of HTTr Catalog concentration-response profiling results
    #' 
    #' @param filepath_catalog character | path of file containing concentration-
    #'  response profiling results as output by httrpathway
    #' @return data.frame of curve fit estimates for all catalog signatures
    #' @example
    #' @export
    # check for file compatibility
    if (!file.exists(filepath_catalog)) {
        stop("Error: Problem with `filepath_catalog`. \n `filepath_catalog` not found")
    } else if (tools::file_ext(filepath_catalog) != "RData") {
        stop("Error: `filepath_catalog` must be an RData file")
    }
    # import + harmonize object name
    load(filepath_catalog)
    cr_catalog <- get(ls()[!ls() %in% c("filepath_catalog")])
    # log-transforma BMD values
    cr_processed <- cr_catalog %>%
        mutate(bmd_log = case_when(
            hitcall >= 0.9 & top_over_cutoff >= 1.5 ~ log10(bmd), TRUE ~ 2.5
        ))
    return(cr_processed)
}

getOverallBMD <- function(cr_catalog, rcas) {
    #' Calculate overall benchmark doses for HTTr catalog signatures
    #' 
    #' For each dtxsid, summary statistics for each profile of BMDs are computed
    #' for use in determining Tier 1 assessment outcomes. Note that chemicals
    #' with <10 bioactive signatures (hitcall >= 0.9 & top_over_cutoff >= 1.5)
    #' are delineated here calculation of medians/modes, as bioactivity among
    #' any RCAS for these chemicals would be considered "selective", and
    #' therefore modes of signature-level bioactivity are set to the default
    #' inactive dose (2.5-log10(uM)).
    #' 
    #' @param cr_catalog data.frame | curve fit estimates for all catalog
    #'  signatures as output by loadCatalogProfiles()
    #' @param rcas list | RCAS object as output by selectRCASGenes()
    #' @return tibble of potency summary statistics for each dtxsid, including
    #'  (1) median log10BMD, (2) 5th percentile log10BMD, (3) mode log10BMD,
    #'  (4) SD log10BMD, (5) mode - 1SD, (6) mode - 2SD
    #' @example
    #' @export
    # calculate median + 5th percentile log(bmd) for all signatures
    overall <- cr_catalog %>%
        mutate(bioactive = hitcall >= 0.9 & top_over_cutoff >= 1.5) %>%
        group_by(dtxsid, bioactive) %>%
        mutate(
            n_bioactive = n(),
            specific_crit = case_when(
                n_bioactive < 10 ~ "n<10",
                TRUE ~ "median_bioactive"
            )
        ) %>%
        group_by(dtxsid) %>%
        filter(
            (specific_crit == "median_bioactive" & bmd_log < 2.5) |
            (specific_crit == "n<10")
        ) %>%
        summarise(
            across(
                .cols = c(
                    sample_id, casrn, name, time,
                    fit_method, n_bioactive, specific_crit
                ),
                .fns = first
            ),
            bmd_log_med = median(bmd_log),
            bmd_log_5 = quantile(bmd_log, probs = c(0.05)),
            bmd_log_abs5 = nth(bmd_log, 5, order_by = bmd_log),
            bmd_log_mode = calcMode(bmd_log),
            bmd_log_sd = sd(bmd_log),
            threshold_1sd = bmd_log_mode - bmd_log_sd,
            threshold_2sd = bmd_log_mode - (2 * bmd_log_sd)
        ) %>%
        mutate(
            ref_class = rcas$httr.wide$ref_class[
                match(dtxsid, rcas$httr.wide$dtxsid)
            ]
        )
    return(overall)
}

loadAllHTS <- function() {}