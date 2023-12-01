# NAM-integration-pilot workflow 5
# Tiered assessment of chemical screening data for selective target perturbagens
library(dplyr)
library(tidyr)
library(tcpl)
library(DBI)
library(data.table)
source("scripts/pipeline_refchem_assignment.R")

selectHTSEndpoints <- function(
    # filepath_rcas_cr,
    cluster_targets,
    db_host = "ccte-mysql-res.epa.gov",
    db_name = "prod_internal_invitrodb_v3_4",
    usernm = "_dataminer",
    passwd = "pass",
    filepath_save_mc5 = "data/examples/invitrodb_v3_4_filtered_mc5.RData",
    filepath_save_ace = "data/examples/invitrodb_v3_4_assay_information.RData",
    filepath_save_sc2 = "data/examples/invitrodb_v3_4_sc2.RData",
    filepath_load_refchem = NULL
) {
    #' Conduct selection of selective endpoints from invitrodb
    #' 
    #' @param filepath_rcas_cr character | path to concentration-response
    #'  estimates for RCAS as output by profileRCAS()
    #' @param cluster_targets character | cluster annotations to select
    #'  endpoints in relation to. Must be one or more of the entries in the
    #'  `cluster_target` column of the refchemdb cluster assignments table as
    #'  output by assignRefChems().
    #' @param db_host character | name of host url
    #' @param db_name character | name of invitrodb instance to query
    #' @param usernm character | username for host url
    #' @param passwd character | password for host url
    #' @param filepath_save_mc5 character | path to save invitrodb
    #'  concentration-response estimates table
    #' @param filepath_save_ace character | path to save invitrodb
    #'  assay-component-endpoint annotation tables
    #' @param filepath_save_sc2 character | path to save invitrodb
    #'  single concentration summary table
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

    # filter mc5 for refchems associated with any cluster_target +
    # annotate with aeid gene/refchemdb target
    refchem_targets <- filter(
        refchem_assign, cluster_target %in% cluster_targets
    )
    mc5_full <- mc5 %>%
        filter(dsstox_substance_id %in% refchem_targets$dsstox_substance_id) %>%
        mutate(
            ref_class = refchem_targets$cluster_target[
                match(dsstox_substance_id, refchem_targets$dsstox_substance_id)
            ],
            ac50_uM = case_when(!is.na(modl_ga) ~ 10^modl_ga),
            acc_uM = case_when(!is.na(modl_acc) ~ 10^modl_acc),
            gene_symbol = gene$gene_symbol[match(aeid, gene$aeid)]
        ) %>%
        group_by(aeid) %>%
        mutate(chems_measured = n()) %>%
        ungroup()

    # find aeids for gene symbols related to refchemdb clusters
    mc5_measured <- data.frame()
    for (class_name in cluster_targets) {
        mc5_target <- filterForTarget(mc5_full, ace, class_name)
        mc5_measured <- rbind(
            mc5_measured, mutate(mc5_target, assay_target = class_name)
        )
    }
    mc5_measured <- mutate(
        mc5_measured,
        hitc = case_when(is.na(hitc) ~ 0L, TRUE ~ as.integer(hitc)),
        in_class = ref_class == assay_target
    )

    # evaluate aeid confidence: classification accuracy
    conf_activity <- evalConfActivity(mc5_measured)

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

filterForTarget <- function(mc5_full, ace, class_name) {
    #' Match invitrodb endpoints with refchemdb cluster annotations
    #' 
    #' @param mc5_full tibble | table of invitrodb concentration-response
    #'  estimates, subset for any chemical in refchemdb cluster assignments.
    #' @param ace data.table | table of invitrodb assay annotations as output by
    #'  loadAnnotations().
    #' @param class_name character | name of cluster annotation to match to.
    #'  Must match a single entry from refchemdb cluster assignment table as
    #'  output by assignRefChems().
    #' @return tibble of invitrodb concentration-response estimates, subset for
    #'  endpoints matching the class_name target and direction of regulation.
    #' @example
    #' @export
    # parse all genes/directions from class_name
    targets_all <- unlist(strsplit(class_name, "\\|"))
    genes_all <- unname(
        sapply(targets_all, function(x) unlist(strsplit(x, "_")))[1, ]
    )
    directions_all <- unname(
        sapply(targets_all, function(x) unlist(strsplit(x, "_")))[2, ]
    )
    # target <- toupper(unlist(strsplit(class_name, "_"))[1])

    # expand terms for other targets in clusters (or related targets)
    if (any(genes_all == "ESR1")) {
        direction_expand <- directions_all[which(genes_all == "ESR1")]
        genes_all <- c(genes_all, "ESR2")
        directions_all <- c(directions_all, direction_expand)
    } else if (any(genes_all == "RARA")) {
        direction_expand <- directions_all[which(genes_all == "RARA")]
        genes_all <- c(genes_all, "RARB", "RARG", "RXRA", "RXRB", "RXRG")
        directions_all <- c(directions_all, rep(direction_expand, 5))
    } else if (any(genes_all == "PPARA")) {
        direction_expand <- directions_all[which(genes_all == "PPARA")]
        genes_all <- c(genes_all, "PPARG", "PPARD")
        directions_all <- c(directions_all, rep(direction_expand, 2))
    } else if (any(genes_all == "AHR")) {
        direction_expand <- directions_all[which(genes_all == "AHR")]
        genes_all <- c(genes_all, "CYP1A1", "CYP1A2")
        directions_all <- c(directions_all, rep(direction_expand, 2))
    }

    # match refchem directions with aenm directions
    directions_aenm <- directions_all
    directions_aenm[directions_aenm == "Positive"] <- "_Agonist|_up"
    directions_aenm[directions_aenm == "Negative"] <- "_Antagonist|_dn"

    # apply additional direction terms for selected targets (for nondirectional aeids)
    directions_aenm[genes_all == "KCNH2"] <- paste0(directions_aenm[genes_all == "KCNH2"], "|_IC_")
    directions_aenm[genes_all == "NR3C1"] <- paste0(directions_aenm[genes_all == "NR3C1"], "|NR_hGR")
    directions_aenm[genes_all == "CYP1A1"] <- paste0(directions_aenm[genes_all == "CYP1A1"], "|_ADME_")
    directions_aenm[genes_all == "CYP1B1"] <- paste0(directions_aenm[genes_all == "CYP1B1"], "|_ADME_")
    
    # search for matching aeids per target
    ace_all <- data.frame()
    for (idx in seq_len(length(genes_all))) {
        gene <- genes_all[idx]
        direction <- directions_aenm[idx]
        ace_ref <- filter(
            ace,
            official_symbol == gene &
            grepl(direction, assay_component_endpoint_name)
        )
        ace_all <- rbind(ace_all, ace_ref)
    }
    # filter mc5 data for all matching aeids + remove additonal assays:
    ## CEETOX: measuring steroidogenesis
    ## ATG_TRANS_dn: noisier/unreliable curve fits across assay mode
    ## NVS_ADME: measuring enzymatic activity
    mc5_refs <- mc5_full %>%
        filter(
            aeid %in% ace_all$aeid &
            !grepl("CEETOX", aenm) &
            !grepl("ATG_(.+)_TRANS_dn", aenm) &
            !grepl("NVS_ADME", aenm)
        )
    return(mc5_refs)
}

evalConfActivity <- function(mc5_measured) {
    #' Evaluate classification accuracy of invitrodb endpoints versus refchemdb
    #' 
    #' @param mc5_measured tibble | table invitrodb concentration-response
    #'  estimates, subset for as output by filterForTarget().
    #' @return tibble of 
    #' @example
    #' @export
    #' 
    #' 
    # annotate with bioactivity classification
    mc5_measured <- mc5_measured %>%
        mutate(bioactive = (hitc == 1 & use.me == 1))

    # remove duplicate chemical-endpoint pairs (for sc data, keeping highest resp) +
    # count bioactivity classifications for in-class/out-of-class refs
    refs_summary <- mc5_measured %>%
        arrange(dsstox_substance_id, aenm, ref_class, hitc) %>%
        distinct(dsstox_substance_id, aenm, ref_class, .keep_all = TRUE) %>%
        group_by(aenm, assay_target, in_class, bioactive) %>%
        summarise(count = n(), .groups = "drop") %>%
        complete(aenm, assay_target, in_class, bioactive, fill = list(count = 0)) %>%
        mutate(
            class_label = case_when(
                in_class ~ "in-class references",
                TRUE ~ "out-of-class references"
            )
        )
    return(refs_summary)
}