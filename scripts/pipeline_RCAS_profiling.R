# NAM-integration-pilot workflow 3
# Concentration-Response Profiling of Reference Chemical-Associated Signatures

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

makeCatalog <- function(filepath_rcas, random_sigs = TRUE) {
    #' Create custom signature catalog from RCAS results
    #' 
    #' Script pulls RCAS results generated from analyzeHTTrANOVA() and
    #' compiles genes into signatures for use in concentration-response
    #' profiling.
    #' 
    #' @param filepath_rcas character | path to file containing RCAS results
    #'  as output by analyzeHTTrANOVA()
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