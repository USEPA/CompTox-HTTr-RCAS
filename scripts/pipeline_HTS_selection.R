# NAM-integration-pilot workflow 5
# Tiered assessment of chemical screening data for selective target perturbagens
library(dplyr)
library(tidyr)
library(tcpl)
library(DBI)
library(data.table)
source("scripts/pipeline_refchem_assignment.R")

selectHTSEndpoints <- function(
    filepath_rcas_cr,
    db_host = "ccte-mysql-res.epa.gov",
    db_name = "prod_internal_invitrodb_v3_4",
    usernm = "_dataminer",
    passwd = "pass",
    filepath_save_mc5 = "data/examples/invitrodb_v3_4_filtered_mc5.RData",
    filepath_save_ace = "data/examples/invitordb_v3_4assay_information.RData",
    filepath_load_refchem = NULL
) {
    #' Conduct selection of selective endpoints from invitrodb
    #' 
    #' @param filepath_rcas_cr character | path to concentration-response
    #'  estimates for RCAS as output by profileRCAS()
    #' @param db_host character | name of host url
    #' @param db_name character | name of invitrodb instance to query
    #' @param usernm character | username for host url
    #' @param passwd character | password for host url
    #' @param filepath_save_mc5 character | path to save invitrodb
    #'  concentration-response estimates table
    #' @param filepath_save_ace character | path to save invitrodb
    #'  assay-component-endpoint annotation tables
    #' @param filepath_load_refchem character | path containing refchemdb
    #'  cluster assignments as output by assignRefChems(). If null, cluster
    #'  assignments will be generated and into the `data/examples` directory.
    #' @return 
    #' @example
    #' @export
    # load MC5 table + add flag-based filters + load ace tables
    if (file.exists(filepath_save_mc5)) {
        message("loading mc5 table from %s", filepath_save_mc5)
        load(filepath_save_mc5)
    } else {
        message("%s not found; loading mc5 table from invitrodb")
        mc5 <- loadMC5(db_host, db_name, usernm, passwd)
        mc5 <- filterMC5(mc5, filepath_save_mc5)
    }
    if (!file.exists(filepath_save_ace)) {
        message("%s not found; loading ace/gene/cytotox tables from invitrodb")
        loadAnnotations(filepath_save_ace, db_host, db_name, usernm, passwd)
    }
    load(filepath_save_ace)

    # load refchem assignments or generate if null
    if (!is.null(filepath_load_refchem)) {
        load(filepath_load_refchem)
    } else {
        message("no refchem filepath specified; generating assignments from RefChemDB")
        refchem_assign <- assignRefChems()
    }

    # load RCAS
}

loadMC5 <- function(
    db_host = "ccte-mysql-res.epa.gov",
    db_name = "prod_internal_invitrodb_v3_4",
    usernm = "_dataminer",
    passwd = "pass"
) {
    #' Load multi-concentration potency estimates from invitrodb
    #' 
    #' @param db_host character | name of host url
    #' @param db_name character | name of invitrodb instance to query
    #' @param usernm character | username for host url
    #' @param passwd character | password for host url
    #' @return data.table of potency estimates for all spid-aeid combinations 
    #'  (MC5 table)
    #' @example
    #' @export
    # initialize connection
    tcplConf(
        user = usernm,
        pass = passwd,
        db = db_name,
        drvr = "MySQL",
        host = db_host
    )
    con <- dbConnect(
        drv = RMySQL::MySQL(),
        user = usernm,
        password = passwd,
        host = db_host,
        database = db_name
    )
    # import mc5 table + subset duplicate chids
    mc5 <- tcplPrepOtpt(tcplLoadData(lvl = 5, type = "mc"))
    mc5 <- tcplSubsetChid(mc5)
    # import mc6 QC flags
    mc6 <- tcplPrepOtpt(tcplLoadData(lvl = 6, fld = "m4id", val = mc5$m4id, type = "mc"))
    # parse QC flags + annotate mc5
    setDT(mc6)
    mc6_mthds <- mc6[, .(mc6_mthd_id = paste(mc6_mthd_id, collapse = ",")), by = m4id]
    mc6_flags <- mc6[, .(flag = paste(flag, collapse = ";")), by = m4id]
    mc5$mc6_flags <- mc6_mthds$mc6_mthd_id[match(mc5$m4id, mc6_mthds$m4id)]
    mc5[, flag.length := ifelse(!is.na(mc6_flags), count.fields(textConnection(mc6_flags), sep = ","), NA)]
    # parse assay name
    mc5[, asnm := tstrsplit(aenm, "_", fixed = TRUE, keep = c(1))]
    return(mc5)
}

filterMC5 <- function(
    mc5, filepath_save = "data/invitrodb_v3_4_filtered_mc5.RData"
) {
    #' Perform flag-based filtering of invitrodb curve-fits
    #' 
    #' @param filepath_save character | filepath to save RData file
    #' @param mc5 data.table | MC5 table as output by loadMC5()
    #' @return RData file containing MC5 data.table with additional column
    #'  "use.me", denoting bioactive and QC-pass rows (1) versus inactive/
    #'  QC-fail rows (0)
    #' @example
    #' @export
    mc5[hitc == 1 & flag.length < 3, use.me := 1]
    mc5[hitc == 1 & is.na(flag.length), use.me := 1]
    mc5[hitc == 1 & flag.length >= 3, use.me := 0]
    mc5[fitc %in% c(36, 45), use.me := 0]
    mc5[asnm == "BSK" & hitc == 1, use.me := 1]
    mc5[hitc == -1, use.me := 0] # make hitc interpretable as a positive sum
    mc5[use.me == 0, modl_ga := as.numeric(NA)]
    mc5[use.me == 0, hitc := 0]
    mc5[hitc == 0, modl_ga := as.numeric(NA)]
    mc5[hitc == 0, modl_acc := as.numeric(NA)]
    # get curve fits in uM
    mc5[, ac50_uM := ifelse(!is.na(modl_ga), 10^modl_ga, NA)]
    mc5[, acc_uM := ifelse(!is.na(modl_acc), 10^modl_acc, NA)]
    # save to specified path
    save(mc5, file = filepath_save)
    return(mc5)
}

loadAnnotations <- function(
    filepath_save = "data/invitrodb_v3_4_assay_information.RData",
    db_host = "ccte-mysql-res.epa.gov",
    db_name = "prod_internal_invitrodb_v3_4",
    usernm = "_dataminer",
    passwd = "pass"
) {
    #' Load assay-component-endpoint annotations and gene symbols from invitrodb
    #' 
    #' @param filepath_save character | filepath to save RData file
    #' @param db_host character | name of host url
    #' @param db_name character | name of invitrodb instance to query
    #' @param usernm character | username for host url
    #' @param passwd character | password for host url
    #' @return RData file containing data.tables of annotations for all aeids
    #'  (ace), gene symbols (gene), and cytotoxicity estimates for dtxsids
    #'  (cytotox)
    #' @example
    #' @export
    # initialize connection
    tcplConf(
        user = usernm,
        pass = passwd,
        db = db_name,
        drvr = "MySQL",
        host = db_host
    )
    con <- dbConnect(
        drv = RMySQL::MySQL(),
        user = usernm,
        password = passwd,
        host = db_host,
        database = db_name
    )
    # query invitrodb for assay component endpoint (ace) annotations
    ace <- dbGetQuery(
        con,
        "select assay_component_endpoint.*
        FROM prod_internal_invitrodb_v3_4.assay_component_endpoint;"
    ) %>%
        data.table()

    gene <- dbGetQuery(
        con,
        "SELECT * FROM prod_internal_invitrodb_v3_4.intended_target
        INNER JOIN prod_internal_invitrodb_v3_4.gene
        ON intended_target.target_id=gene.gene_id;"
    ) %>%
        data.table()

    ace <- merge.data.table(
        ace,
        gene[, c("aeid", "official_symbol", "gene_name")],
        by = "aeid",
        all.x = TRUE
    )

    cytotox <- dbGetQuery(
        con,
        "SELECT * FROM prod_internal_invitrodb_v3_4.cytotox
        INNER JOIN prod_internal_invitrodb_v3_4.chemical
        ON cytotox.chid=chemical.chid;"
    ) %>%
        data.table()
    # remove duplicate column `chid`
    cytotox[, 1] <- NULL
    # set default cytotox value to 2.5 (lower bound)
    cytotox$cytotox_lower_bound_log[
        cytotox$cytotox_lower_bound_log > 2.5
    ] <- 2.5
    # save to specified directory
    save(ace, gene, cytotox, file = filepath_save)
    return()
}