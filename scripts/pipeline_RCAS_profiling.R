# NAM-integration-pilot workflow 3
# Concentration-Response Profiling of Reference Chemical-Associated Signatures
library(dplyr)
library(tidyr)
library(openxlsx)
# library(httrlib)
library(httrpathway)
# source("../httrpathway/httrpl/Rlib/httrpl.R", chdir = TRUE)

profileRCAS <- function(db_host, db_name, cytotox = FALSE) {
    #' Runtime function for RCAS concentration-response profiling
    #' 
    #' @param db_host character | URL of server hosting mongo database
    #' @param db_name character | valid name of study to pull data from
    #' @param cytotox logical | flag specifying whether to load chemical
    #'  cytotoxicity estimates for the specified study
    #' @return saved intermediate/output files: FCMAT1 [probe-level
    #'  log2(fold-change) values], FCMAT2 [gene-level log2(fold-change) values],
    #'  SIGNATURE_CR [signature-level concentration-response curve fits], and
    #'  associated plots
    #' @example
    #' @export
    # check for directory structure and make directory if needed
    makeDirectory()
}

makeDirectory <- function() {
    #' Create directory structure for RCAS profiling
    #' 
    #' Script creates series of directories as expected by httrpathway for
    #' loading input data and saving output data. See documentation of
    #' httrpathway package for further details.
    #' 
    #' @return directory tree for input/output data
    #' @example makeDirectory("./data/")
    #' @export
    directory_name <- "../"
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
    #'  filepath_rcas <- "data/examples/rcas_u2os_gene.RData"
    #'  rcas <- selectRCASGenes(filepath_rcas)
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
        filter(variable_count >= min_variables)
    
    # add composite to rcas object
    rcas$composite <- composite_unique
    return(rcas)
}

assignGeneDirection <- function(rcas) {
    # filter httr$bmd for genes in each reference class +
    # assign direction based on sign of top +
    # choose winning direction based on count (and mean_top if tied)
    sigdir <- rcas$bmd %>%
        mutate(
            gene_class = rcas$composite$class_1[
                match(gene, rcas$composite$variable)
            ]
        ) %>%
        filter(gene_class == ref_class) %>%
        mutate(
            direction = case_when(
                sign(top) == 1 ~ "up",
                sign(top) == -1 ~ "dn",
                TRUE ~ "nondirectional"
            )
        ) %>%
        group_by(ref_class, gene, direction) %>%
        summarise(count = n(), mean_top = mean(top), .groups = "drop_last") %>%
        mutate(
            keep = case_when(
                count == max(count) ~ TRUE,
                TRUE ~ FALSE
            )
        ) %>%
        filter(keep)

    # annotate composite df with directions +
    # add to rcas object
    composite_dir <- left_join(
        rcas$composite, sigdir,
        by = c("class_1" = "ref_class", "variable" = "gene")
    )
    rcas$composite <- composite_dir
    return(rcas)
}

generateRandomSignatures <- function(rcas, nsigs = 100) {
    gene_lengths <- pull(count(rcas$composite, class_1), n)
    gene_list <- pull(rcas$composite, variable)

    # randomize signature lengths
    sig_lengths <- sample(
        seq(min(gene_lengths), max(gene_lengths), 1), nsigs, replace = TRUE
    )
    # randomly draw genes for each signature length + store as list
    sigs_all <- list()
    for (idx in seq_len(length(sig_lengths))) {
        sig <- sample(gene_list, sig_lengths[idx], replace = FALSE)
        sigs_all[[paste0("Random_", idx)]] <- sig
    }
    return(sigs_all)
}

makeCatalog <- function(rcas, catalog_name = "signatureDB_master_catalog.xlsx", random_sigs = NULL) {
    #' Create custom signature catalog from RCAS results
    #' 
    #' Script pulls RCAS results generated from analyzeHTTrANOVA() and
    #' compiles genes into signatures for use in concentration-response
    #' profiling.
    #' 
    #' @param rcas list | RCAS list object with selected gene lists
    #'  and assigned directionality as output by assignGeneDirection()
    #' @param catalog_name character | name of saved xlsx file containing
    #'  signature catalog. Must end with ".xlsx".
    #' @param random_sigs (optional) list | named vectors of random
    #'  genes as output by generateRandomSignatures()
    #' @return saved xlsx file containing custom signature catalog. File is
    #'  stored in the created input/signatures directory
    #' @example filepath <- "~/HTTr_ANOVA_composite.RData"
    #'  makeCatalog(filepath)
    #' @export
    # check if directory exists, and add "/" to end if not added
    directory_name <- "../"
    dir_catalog <- paste0(directory_name, "input/signatures/")
    # convert rcas to catalog form
    sigcatalog <- rcas$composite %>%
        unite("signature", c(class_1, direction), sep = "_", remove = FALSE) %>%
        group_by(signature) %>%
        summarise(
            ngene = n(),
            gene.list = paste(variable, collapse = "|"),
            parent = paste0(first(class_1)),
            source = "RCAS",
            subsource = "-",
            type = "directional",
            direction = first(direction),
            description = paste0(
                "RCAS for ", first(class_1), ": ", first(direction)
            ),
            target_class = "molecular target",
            super_target = gsub("_.*", "", target_class),
            .groups = "drop"
        )
    # if given random_sigs, expand lists to catalog form
    if (!is.null(random_sigs)) {
        randomcat <- data.frame(
            signature = names(random_sigs),
            ngene = sapply(random_sigs, length),
            gene.list = sapply(random_sigs, function(x) paste(x, collapse = "|")),
            parent = names(random_sigs),
            source = "Random",
            subsource = "-",
            type = "nondirectional",
            direction = "nondirectional",
            description = paste0("Signature for randomization test: ", names(random.heparg)),
            target_class = "Random",
            super_target = "Random"
        )
        sigcatalog <- rbind(sigcatalog, randomcat)
    }
    # remove signatures with <10 genes +
    # rename non-paired directional signatures to nondirectional
    sigcatalog <- sigcatalog %>%
        group_by(parent) %>%
        mutate(
            type = case_when(
                any(ngene < 10) ~ "nondirectional",
                TRUE ~ type
            ),
            direction = case_when(
                any(ngene < 10) ~ "nondirectional",
                TRUE ~ direction
            )
        ) %>%
        ungroup() %>%
        filter(ngene >= 10)

    # export catalog:
    ## signatureDB.RData: subset from signature catalog
    sigdb <- sigcatalog %>%
        select(
            signature, parent, source, subsource, type,
            direction, ngene, description, gene.list
        )
    save(sigdb, file = paste0(dir_catalog, "signatureDB.RData"))

    ## signatureDB_genelists.RData: subset from signature catalog as list
    genelist.df <- sigcatalog %>%
        select(signature, gene.list)
    genelists <- as.list(genelist.df$gene.list)
    names(genelists) <- genelist.df$signature
    genelists <- lapply(genelists, function(x) {unlist(strsplit(x, "\\|"))})

    save(genelists, file = paste0(dir_catalog, "signatureDB_genelists.RData"))

    ## {catalog_name}.xlsx: subset from signature catalog
    sigcatalog <- sigcatalog %>%
        select(
            signature, parent, source, subsource, type,
            direction, ngene, description, super_target,
            target_class
        )
    write.xlsx(sigcatalog, paste0(dir_catalog, catalog_name))
    return()
}

importDESeq2 <- function(db_host, db_name, fc1_name) {
    #' Load gene fold-changes from MongoDB
    #' 
    #' Script pulls DESeq2-moderated fold changes from MongoDB for the specified
    #' study. Gene names are matched to each probe, and column names are
    #' formatted as expected by httrpathway::buildFCMAT1.fromDB.
    #' 
    #' @param db_host character | alias of host MongoDB server
    #' @param db_name character | name of study to pull data from. See
    #'  httr-db Wiki for options, all study names and credentials must
    #'  be stored in a .mongopw file on home directory for access.
    #' @param fc1_name character | name of saved RData file containing
    #'  fold-change matrix. Must end with ".RData".
    #' @return saved RData file containing gene fold-change values for
    #'  all chemicals in selected study. File is stored in the created
    #'  input/fcdata directory.
    #' @example
    #'  db_host <- "ccte-mongodb-res.epa.gov"
    #'  db_name <- "res_httr_u2os_toxcast"      # Must have this study and credentials in .mongopw
    #'  file_fc1 <- "httr_u2os_fc1.RData"
    #'  importDESeq2(db_host, db_name, file_fc1)
    #' @export
    # get chem_ids for full study
    directory_name <- "../"
    db.con <- httrlib::openMongo(
        host = db_host, db = db_name, collection = "httr_chem"
    )
    httr.chem <- db.con$find()

    # get probe_ids matching any rcas gene symbol
    load(paste0(directory_name, "input/signatures/signatureDB_genelists.RData"))
    genes_all <- unname(unlist(genelists))
    probes <- httrlib::getProbeManifest(db_host = db_host, db_name = db_name)
    probes_rcas <- probes[probes$gene_symbol %in% genes_all, ]
    
    # retrieve DESeq2 results
    FCMAT1 <- httrlib::getDEGs(
        db_host = db_host,
        db_name = db_name,
        flatten = TRUE,
        chem_id = httr.chem$chem_id[1],
        # probe_id = probes_rcas$probe_name[1],
        anl_name = "meanncnt0_5-plateteffect_1-shrinkage_normal"
    )
    # add missing columns via chemical manifest
    FCMAT1$dtxsid <- httr.chem$dtxsid[match(FCMAT1$chem_id, httr.chem$chem_id)]
    FCMAT1$chem_name <- httr.chem$chem_name[
        match(FCMAT1$chem_id, httr.chem$chem_id)
    ]

    # add gene-level annotations via probe manifest
    FCMAT1$gene_symbol <- probes$gene_symbol[
        match(FCMAT1$probe_id, probes$probe_name)
    ]

    # rename other columns to match expected FCMAT2 names
    colnames.old <- c("trt_grp_id", "log2FoldChange", "lfcSE")
    colnames.new <- c("sample_key", "l2fc", "se")
    colnames(FCMAT1)[which(colnames(FCMAT1) %in% colnames.old)] <- colnames.new

    # save to expected dir for DESeq2 output
    save(FCMAT1, file = paste0(
        directory_name, "input/fcdata/new_versions/", fc1_name
    ))
    return()
}

makeFCMAT2 <- function(db_name, fc1_name) {
    #' Reshape DESeq2-moderated fold-changes for signature profiling
    #' 
    #' Script processes DESeq2-moderated fold changes from MongoDB to
    #'  a form usable by httrpathway::driver(). Serves as a wrapper
    #'  for httrpathway::buildFCMAT1.fromDB() and
    #'  httrpathway::buildFCMAT2.fromDB().
    #' 
    #' @param db_name character | name of study to pull data from. See
    #'  httr-db Wiki for options, all study names and credentials must
    #'  be stored in a .mongopw file on home directory for access.
    #' @param fc1_name character | name of saved RData file containing
    #'  fold-change matrix as output by importDESeq2(). Must end with
    #'  ".RData".
    #' @return saved RData files containing gene fold-change values for
    #'  all chemicals in selected study. Files are stored in the created
    #'  input/fcdata directory.
    #' @example
    #'  directory_name <- "./data/"
    #'  file_fc1 <- "httr_u2os_fc1.RData"
    #'  makeFCMAT2(directory_name, db_host, db_name, file_fc1)
    #' @export
    # CHEM_DICT + FCMAT1 objects
    directory_name <- "../"
    buildFCMAT1.fromDB(
        dataset = db_name,
        dir = paste0(directory_name, "input/fcdata/new_versions/"),
        infile = fc1_name
    )
    # FCMAT2 + PV + SE + FCoverSE + CHEM_DICT objects
    buildFCMAT2.fromDB(
        dataset = db_name,
        dir = paste0(directory_name, "input/fcdata/"),
        method = "gene"
    )
}

importCytotox <- function(db_host, db_name) {
    #' Load chemical cytotoxicity data from MongoDB
    #' 
    #' Script pulls HTPP-derived cytotoxicity flags from MongoDB for the
    #' specified study for quality control. Data are retained for chem_id
    #' instances with at least one "CELL_VIABILITY" flag across the
    #' concentration series, and the lowest concentration containing the flag
    #' is exported.
}