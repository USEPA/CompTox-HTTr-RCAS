# NAM-integration-pilot workflow 5
# Tiered assessment of chemical screening data for selective target perturbagens
library(dplyr)
library(tidyr)
library(openxlsx)
source("scripts/pipeline_RCAS_profiling.R")
source("scripts/pipeline_HTS_selection.R")

runFramework <- function(
    filepath_rcas,
    filename_cr,
    filepath_catalog,
    filepath_mc5_selected,
    filepath_mc5_burst,
    filepath_mc5_ace,
    filepath_export = "data/examples/framework_results.RData"
) {
    #' Runtime function for assessing HTTr/ToxCast chemicals in tiered framework
    #' 
    #' @param filepath_rcas character | path to file containing RCAS
    #'  results as output by analyzeHTTrANOVA()
    #' @param filename_cr character | name of file containing HTTr RCAS
    #'  concentration-response profiles for chemical screen as output by
    #'  profileRCAS(). File must be located in
    #'  "../output/signature_conc_resp_summary/" directory.
    #' @param filepath_catalog character | name of file containing HTTr
    #'  concentration-response profiles for HTTr catalog signatures
    #' @param filepath_mc5_selected character | path to file containing
    #'  invitrodb concentration-response estimates (mc5 table) for orthogonal
    #'  endpoints as output by selectHTSEndpoints()
    #' @param filepath_mc5_burst character | path to file containing
    #'  non-selective point of departure estimates as output by
    #'  selectHTSEndpoints()
    #' @param filepath_mc5_ace character | path to file containing invitrodb
    #'  endpoint annotations as output by loadAnnotations()
    #' @param filepath_export character | path to save assessment outcomes
    #'  for chemical screen
    #' @return tibble of outcomes for all chemicals screened in HTTr and
    #'  ToxCast, designated as "tier1-positive", "tier1-selective",
    #'  "tier2-positive", and "tier2-selective"
    #' @examples
    #'  filepath_rcas <- "data/examples/rcas_u2os_gene.RData"
    #'  filename_cr <- "SIGNATURE_CR_u2os_res_httr_u2os_toxcast_gsea_0.05_conthits.RData"
    #'  filepath_catalog <- "/ccte/projects1/HTTr/screen_signature_cr/signature_conc_resp/SIGNATURE_CR_screen_large_u2os_toxcast_pfas_pe1_normal_v2_gsea_0.05_conthits.RData"
    #'  filepath_mc5_selected <- "data/examples/invitrodb_v3_4_selected.RData"
    #'  filepath_mc5_burst <- "data/examples/invitrodb_v3_4_burst.RData"
    #'  filepath_mc5_ace <- "data/examples/invitrodb_v3_4_assay_information.RData"
    #'  results <- runFramework(
    #'      filepath_rcas,
    #'      filename_cr,
    #'      filepath_catalog,
    #'      filepath_mc5_selected,
    #'      filepath_mc5_burst,
    #'      filepath_mc5_ace
    #'  )
    #' @export
    # load Tier 1 data:
    ## RCAS object
    message(gettextf("loading RCAS object from %s", filepath_rcas))
    rcas <- selectRCASGenes(filepath_rcas)
    ## RCAS profiling data + remove random signatures
    message(gettextf("loading RCAS cr estimates from %s", filename_cr))
    cr_rcas <- loadRCASprofiles(filename_cr)
    cr_rcas <- filter(cr_rcas, !grepl("Random", signature))
    ## Catalog profiling data
    message(gettextf(
        "loading HTTr catalog cr estimates from %s", filepath_catalog
    ))
    cr_catalog <- loadCatalogProfiles(filepath_catalog)
    ## calculate HTTr non-selective points of departure
    message("summarizing HTTr non-selective PODs")
    cr_burst <- getOverallBMD(cr_catalog, rcas)

    # load Tier 2 data:
    ## MC5 data for orthogonal endpoints
    if (file.exists(filepath_mc5_selected)) {
        message(gettextf("loading mc5 table from %s", filepath_mc5_selected))
        load(filepath_mc5_selected)
    } else {
        stop(gettextf(
            "%s not found; run selectHTSEndpoints() to generate mc5 table",
            filepath_mc5_selected
        ))
    }
    ## ace/gene/cytotox tables
    if (file.exists(filepath_mc5_ace)) {
        message(gettextf("loading ace/gene/cytotox tables from %s", filepath_mc5_ace))
        load(filepath_mc5_ace)
    } else {
        stop(gettextf(
            "%s not found; run loadAnnotations() to generate ace/gene/cytotox tables",
            filepath_mc5_ace
        ))
    }
    ## ToxCast non-selective points of departure
    if (file.exists(filepath_mc5_burst)) {
        message(gettextf(
            "loading invitrodb non-selective PODs tables from %s",
            filepath_mc5_burst
        ))
        load(filepath_mc5_burst)
    } else {
        stop(gettextf(
            "%s not found; run selectHTSEndpoints() to generate burst table",
            filepath_mc5_burst
        ))
    }
}

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