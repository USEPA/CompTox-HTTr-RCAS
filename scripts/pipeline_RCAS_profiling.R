# NAM-integration-pilot workflow 3
# Concentration-Response Profiling of Reference Chemical-Associated Signatures

profileRCAS <- function(
    db_host, db_name, directory_name = getwd(), cytotox = FALSE
) {
    #' Runtime function for RCAS concentration-response profiling
}

makeDirectory <- function(directory_name) {
    #' Create directory structure for RCAS profiling
    #' 
    #' Script creates series of directories as expected by httrpathway for
    #' loading input data and saving output data. See documentation of
    #' httrpathway package for further details.
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