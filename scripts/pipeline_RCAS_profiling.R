# NAM-integration-pilot workflow 3
# Concentration-Response Profiling of Reference Chemical-Associated Signatures
library(dplyr)
library(tidyr)

profileRCAS <- function(
    db_host, db_name, directory_name = "./data/", cytotox = FALSE
) {
    #' Runtime function for RCAS concentration-response profiling
    #' 
    #' @param db_host character | URL of server hosting mongo database
    #' @param db_name character | valid name of study to pull data from
    #' @param directory_name character | name of top-level directory to 
    #'  construct sub-directories into. defaults to /data subdirectory
    #' @param cytotox logical | flag specifying whether to load chemical
    #'  cytotoxicity estimates for the specified study
    #' @return saved intermediate/output files: FCMAT1 [probe-level 
    #'  log2(fold-change) values], FCMAT2 [gene-level log2(fold-change) values],
    #'  SIGNATURE_CR [signature-level concentration-response curve fits], and 
    #'  associated plots
    #' @example
    #' @export
    # check for directory structure and make directory if needed
    makeDirectory(directory_name)
}

makeDirectory <- function(directory_name) {
    #' Create directory structure for RCAS profiling
    #' 
    #' Script creates series of directories as expected by httrpathway for
    #' loading input data and saving output data. See documentation of
    #' httrpathway package for further details.
    #' 
    #' @param directory_name character | name of top-level directory to 
    #'  construct sub-directories into
    #' @return directory tree for input/output data
    #' @example makeDirectory("./data/")
    #' @export
    # check if directory exists, and add "/" to end if not added
    if (!dir.exists(directory_name)) stop(gettextf("%s does not exist"))
    if (!grepl("/$", directory_name)) {
        directory_name <- paste0(directory_name, "/")
    }
    # check for input directories, and create if empty:
    # input/
    # input/fcdata
    # input/fcdata/new_versions
    # input/signatures
    dir_input <- paste0(directory_name, "input/")
    sub_input <- paste0(
        dir_input, c("fcdata/", "fcdata/new_versions", "signatures/")
    )
    for (subdir in c(dir_input, sub_input)) {
        if (!dir.exists(subdir)) {
            dir.create(subdir)
        } else {
            message(gettextf("makeDirectory: %s already exists", subdir))
        }
    }
    # repeat for output directories:
    # output/
    # output/gene_conc_resp_plots
    # output/gene_conc_resp_summary
    # output/signature_conc_resp_plots
    # output/signature_conc_resp_summary
    # output/signature_cutoff
    # output/signature_score_summary
    # output/super_target_boxplot
    dir_output <- paste0(directory_name, "output/")
    sub_output <- paste0(
        dir_output, c(
            "gene_conc_resp_plots", "gene_conc_resp_summary",
            "signature_conc_resp_plots", "signature_conc_resp_summary",
            "signature_cutoff", "signature_score_summary",
            "super_target_boxplot"
        )
    )
    for (subdir in c(dir_output, sub_output)) {
        if (!dir.exists(subdir)) {
            dir.create(subdir)
        } else {
            message(gettextf("makeDirectory: %s already exists", subdir))
        }
    }
    return()
}

selectRCASGenes <- function(filepath_rcas, prop_class = 0.7, min_variables = 10) {
    #' Perform gene selection for final RCAS determination
    #' 
    #' Gene selection is performed for each reference class to
    #' identify genes that have a significantly lower potency
    #' than a proportion of other reference classes, forming a
    #' final signature.
    #' 
    #' @param filepath_rcas character | path to RData file containing RCAS
    #'  results as output by analyzeHTTrANOVA()
    #' @param prop_class numeric | proportion of classes that a class
    #'  must have a significantly lower mean BMD than for each gene tested.
    #'  Number must be between 0-1.
    #' @param min_variables numeric | minimum number of unique genes needed to
    #'  constitute a final RCAS
    #' @return list | RCAS list object with additional data.frame named
    #'  "composite", in which final RCAS genes are listed for each class
    #' @example
    #' @export
    # check for file compatibility
    if (!file.exists(filepath_rcas)) {
        stop("Error: Problem with `filepath_rcas`. \n `filepath_rcas` not found")
    } else if (tools::file_ext(filepath_rcas) != "RData") {
        stop("Error: `filepath_rcas` must be an RData file")
    }
    # import + harmonize object name
    load(filepath_rcas)
    rcas <- get(ls()[!ls() %in% c("filepath_rcas", "random_sigs")])

    # determine number of pairwise comparisons for which a gene
    # is more potent for a reference class, calculate propotion of
    # significant/total comparisons, and keep genes with >70%
    # significant pairwise comparisons
    composite <- rcas$posthoc.estimate$posthoc.sig %>%
        mutate(n_classes = length(unique(rcas$httr.wide$ref_class))) %>%
        group_by(class_1, variable) %>%
        mutate(n_sig = n(), prop_sig = n_sig / (n_classes - 1)) %>%
        filter(prop_sig > prop_class & variable != "global") %>%
        summarise(
            across(
                diff,
                list(mean = mean, med = median, sd = sd),
                .names = "{.col}_{.fn}"
            ),
            n_sig = first(n_sig),
            prop_sig = first(prop_sig),
            mean_bmd_diff = first(diff_mean),
            mean_bmd_log = first(mean_1),
            # mean_global = first(global_1),
            .groups = "drop"
        )
    # filter for genes unique to 1 class +
    # classes with minumum number of unqiue genes
    composite_unique <- composite %>%
        filter(mean_bmd_log <= 1.5) %>%
        group_by(variable) %>%
        mutate(class_count = n()) %>%
        filter(class_count == 1) %>%
        group_by(class_1) %>%
        mutate(variable_count = n()) %>%
        ungroup() %>%
        filter(variable_count > min_variables)
    
    # add composite to rcas object
    rcas$composite <- composite_unique
    return(rcas)
}

makeCatalog <- function(rcas, random_sigs = TRUE) {
    #' Create custom signature catalog from RCAS results
    #' 
    #' Script pulls RCAS results generated from analyzeHTTrANOVA() and
    #' compiles genes into signatures for use in concentration-response
    #' profiling.
    #' 
    #' @param rcas list | RCAS list object with selected gene lists
    #'  as output by selectRCASGenes()
    #' @return saved xlsx file containing custom signature catalog. File is
    #'  stored in the created input/signatures directory
    #' @example filepath <- "~/HTTr_ANOVA_composite.RData"
    #'  makeCatalog(filepath)
    #' @export
    
}

importDESeq2 <- function(db_host, db_name) {
    #' Load gene fold-changes from MongoDB
    #' 
    #' Script pulls DESeq2-moderated fold changes from MongoDB for the specified
    #' study. Gene names are matched to each probe, and column names are
    #' formatted as expected by httrpathway::buildFCMAT1.fromDB.
}

loadCytotox <- function(db_host, db_name) {
    #' Load chemical cytotoxicity data from MongoDB
    #' 
    #' Script pulls HTPP-derived cytotoxicity flags from MongoDB for the
    #' specified study for quality control. Data are retained for chem_id
    #' instances with at least one "CELL_VIABILITY" flag across the
    #' concentration series, and the lowest concentration containing the flag
    #' is exported.
}