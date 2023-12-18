# NAM-integration-pilot workflow 5
# Tiered assessment of chemical screening data for selective target perturbagens
library(dplyr)
library(tidyr)
library(openxlsx)
source("scripts/pipeline_RCAS_generation.R")
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

    
    # subset Tier 1/2 data for dtxsids profiled in Tier 1 RCAS
    cr_catalog_rcas <- filter(cr_catalog, dtxsid %in% cr_rcas$dtxsid)
    cr_burst_rcas <- filter(cr_burst, dtxsid %in% cr_rcas$dtxsid)
    chems_filtered_rcas <- filter(
        chems_filtered, dsstox_substance_id %in% cr_rcas$dtxsid
    )
    mc5_burst_rcas <- filter(mc5_burst, dsstox_substance_id %in% cr_rcas$dtxsid)

    # subset Tier 2 data for endpoints matching each Tier 1 RCAS
    chems_match_rcas <- harmonizeTargetNames(cr_rcas, chems_filtered_rcas)

    # combine Tier 1/2 tables (all chemicals x all measured endpoints)
    compare_tier1 <- combineTier1(cr_rcas, cr_burst_rcas)
    compare_tier1 <- mutate(
        compare_tier1,
        # signature_label = gsub("^(.+?)_", "", gsub("_up|_dn", "", signature))
        signature_label = gsub("_up|_dn", "", signature)
    )
    compare_tier2 <- combineTier2(chems_match_rcas, mc5_burst_rcas)
    compare_full <- compare_tier1 %>%
        select(
            dtxsid, name, signature, signature_label, bmd_log, bmd_log_med,
            bmd_log_mode, threshold_1sd, tier1_positive, tier1_selective_1sd
        ) %>%
        full_join(
            ., select(
                compare_tier2,
                dsstox_substance_id, chnm, n_toxcast_meas, aeid, aenm, hitc,
                use.me, modl_acc, specific_crit, modl_acc_mode, tier2_positive,
                tier2_selective, tier2_conflict, signature
            ),
            by = c(
                "dtxsid" = "dsstox_substance_id",
                "signature_label" = "signature"
            )
        ) %>%
        rename(tier1_selective = tier1_selective_1sd)
    # condense table to chemicals x rcas
    compare_sum <- compare_full %>%
        group_by(dtxsid, signature) %>%
        summarise(
            across(
                c(
                    name, chnm, bmd_log, bmd_log_med, bmd_log_mode,
                    threshold_1sd, tier1_positive, tier1_selective,
                    n_toxcast_meas
                ),
                unique
            ),
            across(c(modl_acc, modl_acc_mode), ~ min(.x, na.rm = TRUE)),
            across(c(tier2_positive, tier2_selective, tier2_conflict), any),
            .groups = "drop"
        )
    # label tier + overall assessment outcomes
    compare_sum <- compare_sum %>%
        mutate(
            outcome_tier1 = case_when(
                tier1_positive & tier1_selective ~ "selective",
                tier1_positive & !tier1_selective ~ "positive",
                TRUE ~ "negative"
            ),
            outcome_tier2 = case_when(
                tier2_positive & tier2_selective ~ "selective",
                tier2_positive & !tier2_selective ~ "positive",
                TRUE ~ "negative"
            ),
            outcome_all = case_when(
                outcome_tier1 == "selective" & outcome_tier2 == "selective" ~ "tier2_selective",
                outcome_tier1 == "selective" & outcome_tier2 == "positive" ~ "tier2_positive",
                outcome_tier1 == "selective" & outcome_tier2 == "negative" ~ "tier1_selective",
                outcome_tier1 == "positive" ~ "tier1_positive",
                outcome_tier1 == "negative" ~ "tier1_negative"
            )
        )
    # export full/summary tables to file
    save(compare_full, compare_sum, file = filepath_export)
    return(compare_sum)
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

harmonizeTargetNames <- function(cr_rcas, chems_filtered_rcas) {
    #' Harmonize names of target classes between Tier 1/2 data
    #' 
    #' @param cr_rcas
    #' @param chems_filtered_rcas
    #' @return
    #' @example
    #' @export
    # isolate rcas names + mc5 `assay_target` names
    names_rcas <- gsub(
        "u2os_|heparg_|mcf7_",
        "",
        gsub("_up|_dn", "", unique(cr_rcas$signature))
    )
    names_mc5 <- data.frame(
        assay_target = unique(chems_filtered_rcas$assay_target)
    )
    # convert mc5 assay_targets to rcas names
    rcas_mc5 <- cleanRefLabels(names_mc5, "assay_target")
    colnames(rcas_mc5) <- "rcas_label"
    rcas_mc5 <- cbind(names_mc5, rcas_mc5)
    # filter mc5 names for those in httr
    rcas_filtered <- filter(rcas_mc5, rcas_label %in% names_rcas)
    # filter mc5 data for endpoints matching httr rcas
    chems_match_rcas <- chems_filtered_rcas %>%
        filter(assay_target %in% rcas_filtered$assay_target) %>%
        mutate(assay_rcas = rcas_filtered$rcas_label[
            match(assay_target, rcas_filtered$assay_target)
        ])
    return(chems_match_rcas)
}

combineTier1 <- function(cr_rcas, cr_burst_rcas) {
    #' Combine Tier 1 RCAS/burst readouts
    #' 
    #' @param cr_rcas
    #' @param cr_burst_rcas
    #' @return
    #' @example
    #' @export
    # combine RCAS/burst measurements + define Tier 1 criteria
    compare <- cr_rcas %>%
        mutate(
            bmd_log = case_when(
                hitcall >= 0.9 & top_over_cutoff >= 1.5 ~ log10(bmd),
                TRUE ~ 2.5
            )
        ) %>%
        group_by(dtxsid, signature) %>%
        slice_min(bmd_log, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        left_join(
            ., select(
                cr_burst_rcas,
                dtxsid, specific_crit, bmd_log_med, bmd_log_mode, threshold_1sd
            ),
            by = "dtxsid"
        ) %>%
        mutate(
            specific_crit = case_when(is.na(specific_crit) ~ "n<10", TRUE ~ specific_crit),
            bmd_log_med = case_when(is.na(bmd_log_med) ~ 2.5, TRUE ~ bmd_log_med),
            bmd_log_mode = case_when(is.na(bmd_log_mode) ~ 2.5, TRUE ~ bmd_log_mode),
            threshold_1sd = case_when(is.na(threshold_1sd) ~ 2, TRUE ~ threshold_1sd),
            tier1_positive = bmd_log < 2.5,
            tier1_selective_1sd = bmd_log < threshold_1sd
        )
    return(compare)
}

combineTier2 <- function(chems_match_rcas, mc5_burst_rcas, threshold_specific = 0.33) {
    #' Combine Tier 2 targeted/burst readouts
    #' 
    #' @param chems_match_rcas
    #' @param mc5_burst_rcas
    #' @param threshold_specific double | log10-difference between targeted
    #'  ACC and statistical mode used to denote tier2-selective chemicals
    #' @return
    #' @example
    #' @export
    # compile number of orthogonal ToxCast endopints measured per dtxsid/target
    chems_endpts_meas <- chems_match_rcas %>%
        distinct(dsstox_substance_id, aenm, assay_rcas, .keep_all = TRUE) %>%
        count(dsstox_substance_id, assay_rcas, name = "n_toxcast_meas")
    chems_match_rcas <- left_join(
        chems_match_rcas, chems_endpts_meas,
        by = c("dsstox_substance_id", "assay_rcas")
    )

    # combine targeted/burst measurements + define Tier 2 criteria
    compare_tier2 <- left_join(
        chems_match_rcas, select(mc5_burst_rcas, -spid),
        by = c("dsstox_substance_id", "casn", "chnm")
    ) %>%
        mutate(
            tier2_positive = !is.na(use.me) & use.me == 1,
            tier2_selective = case_when(
                specific_crit == "median_bioactive" ~ modl_acc <= (modl_acc_mode - threshold_specific),
                specific_crit == "n<10" ~ modl_acc <= 2.5 - threshold_specific
            ),
            tier2_selective = case_when(
                !is.na(tier2_selective) ~ tier2_selective, TRUE ~ FALSE
            ),
            signature = assay_rcas
            ) %>%
        group_by(dsstox_substance_id, signature) %>%
        mutate(
            tier2_n_specific = sum(tier2_selective),
            tier2_n_measured = n(),
            tier2_conflict = any(tier2_selective) & !all(tier2_selective)
        ) %>%
        ungroup()
    return(compare_tier2)
}