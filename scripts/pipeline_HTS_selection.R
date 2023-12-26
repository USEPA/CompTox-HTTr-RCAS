# NAM-integration-pilot workflow 5
# Tiered assessment of chemical screening data for selective target perturbagens
library(dplyr)
library(tidyr)
library(tcpl)
library(DBI)
library(data.table)
source("scripts/pipeline_refchem_assignment.R")

selectHTSEndpoints <- function(
    cluster_targets,
    filepath_export = "data/examples/invitrodb_v3_4_selected.RData",
    db_host = "ccte-mysql-res.epa.gov",
    db_name = "prod_internal_invitrodb_v3_4",
    usernm = "_dataminer",
    passwd = "pass",
    filepath_save_mc5 = "data/examples/invitrodb_v3_4_filtered_mc5.RData",
    filepath_save_ace = "data/examples/invitrodb_v3_4_assay_information.RData",
    filepath_save_sc2 = "data/examples/invitrodb_v3_4_sc2.RData",
    filepath_save_burst = "data/examples/invitrodb_v3_4_burst.RData",
    filepath_load_refchem = NULL
) {
    #' Conduct selection of selective endpoints from invitrodb
    #' 
    #' @param cluster_targets character | cluster annotations to select
    #'  endpoints in relation to. Must be one or more of the entries in the
    #'  `cluster_target` column of the refchemdb cluster assignments table as
    #'  output by assignRefChems().
    #' @param filepath_export character | path to save final table of
    #'  invitrodb bioactivity estimates for selected endpoints
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
    # load MC5 table + add flag-based filters
    if (file.exists(filepath_save_mc5)) {
        message(gettextf("loading mc5 table from %s", filepath_save_mc5))
        load(filepath_save_mc5)
    } else {
        message(gettextf(
            "%s not found; loading mc5 table from invitrodb", filepath_save_mc5
        ))
        mc5 <- loadMC5(db_host, db_name, usernm, passwd)
        mc5 <- filterMC5(mc5, filepath_save_mc5)
    }
    # load ace/gene/cytotox tables
    if (!file.exists(filepath_save_ace)) {
        message(gettextf(
            "%s not found; loading ace/gene/cytotox tables from invitrodb",
            filepath_save_ace
        ))
        loadAnnotations(filepath_save_ace, db_host, db_name, usernm, passwd)
    }
    load(filepath_save_ace)
    # load SC2 table
    if (file.exists(filepath_save_sc2)) {
        message(gettextf("loading sc2 table from %s", filepath_save_sc2))
        load(filepath_save_sc2)
    } else {
        message(gettextf(
            "%s not found; loading sc2 table from invitrodb",
            filepath_save_sc2
        ))
        sc2 <- loadSC2(filepath_save_sc2, db_host, db_name, usernm, passwd)
    }

    # load refchem assignments or generate if null
    if (!is.null(filepath_load_refchem)) {
        message(gettextf(
            "loading refchem assignments from %s", filepath_load_refchem
        ))
        load(filepath_load_refchem)
    } else {
        message("no refchem filepath specified; generating assignments from RefChemDB")
        refchem_assign <- assignRefChems()
    }

    # filter mc5/sc2 for refchems associated with any cluster_target +
    # annotate with aeid gene/refchemdb target
    refchem_targets <- filter(
        refchem_assign, cluster_target %in% cluster_targets
    )
    mc5_full <- mc5 %>%
        mutate(
            ref_class = refchem_assign$cluster_target[
                match(dsstox_substance_id, refchem_assign$dsstox_substance_id)
            ],
            ac50_uM = case_when(!is.na(modl_ga) ~ 10^modl_ga),
            acc_uM = case_when(!is.na(modl_acc) ~ 10^modl_acc),
            gene_symbol = gene$gene_symbol[match(aeid, gene$aeid)]
        ) %>%
        group_by(aeid) %>%
        mutate(chems_measured = n()) %>%
        ungroup()

    sc2_full <- sc2 %>%
        mutate(
            ref_class = refchem_assign$cluster_target[
                match(dsstox_substance_id, refchem_assign$dsstox_substance_id)
            ],
            gene_symbol = gene$gene_symbol[match(aeid, gene$aeid)]
        ) %>%
        group_by(aeid) %>%
        mutate(chems_measured = n()) %>%
        ungroup()
    
    # calculate non-selective burst values + save to file
    mc5_burst <- calcBurstMC5(mc5_full)
    save(mc5_burst, file = filepath_save_burst)
    message(gettextf("mc5_burst saved to %s", filepath_save_burst))

    # find aeids for gene symbols related to refchemdb clusters
    mc5_measured <- data.frame()
    sc2_measured <- data.frame()
    for (class_name in cluster_targets) {
        mc5_target <- filterForTarget(mc5_full, ace, class_name)
        sc2_target <- filterForTarget(sc2_full, ace, class_name)
        mc5_measured <- rbind(
            mc5_measured, mutate(mc5_target, assay_target = class_name)
        )
        sc2_measured <- rbind(
            sc2_measured, mutate(sc2_target, assay_target = class_name)
        )
    }
    # combine mc5/sc2 tables + determine in-class/out-of-class labels
    chems_combined <- full_join(
        mc5_measured,
        sc2_measured,
        by = c(
            "chid", "casn", "chnm", "dsstox_substance_id", "ref_class",
            "aeid", "aenm", "assay_target"
        ),
        suffix = c("", "_sc2")
    ) %>%
        mutate(
            hitc = case_when(is.na(hitc) ~ 0L, TRUE ~ as.integer(hitc)),
            hitc_sc2 = case_when(is.na(hitc_sc2) ~ 0L, TRUE ~ as.integer(hitc_sc2)),
            hitc_any = hitc == 1 | hitc_sc2 == 1,
            in_class = ref_class == assay_target
        )
    chems_combined_ref <- filter(
        chems_combined,
        dsstox_substance_id %in% refchem_targets$dsstox_substance_id
    )

    # evaluate aeid confidence:
    # classification accuracy + potency of active chemicals
    conf_activity <- evalConfActivity(chems_combined_ref)
    conf_potency <- evalConfPotency(chems_combined_ref)

    # determine aeids that don't pass criteria for bioactivity + potency
    endpoints_remove <- selectEndpoints(conf_activity, conf_potency)

    # remove non-passing endpoints from mc5/sc2 data + export to file
    chems_filtered <- filter(chems_combined, !aenm %in% endpoints_remove)
    save(chems_filtered, file = filepath_export)
    message(gettextf("passing endpoint data saved to %s", filepath_export))

    return(chems_filtered)
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

loadSC2 <- function(
    filepath_save = "data/examples/invitrodb_v3_4_sc2.RData",
    db_host = "ccte-mysql-res.epa.gov",
    db_name = "prod_internal_invitrodb_v3_4",
    usernm = "_dataminer",
    passwd = "pass"
) {
    #' Load single-concentration bioactivity estimates from invitrodb
    #' 
    #' @param filepath_save character | path to save RData file
    #' @param db_host character | name of host url
    #' @param db_name character | name of invitrodb instance to query
    #' @param usernm character | username for host url
    #' @param passwd character | password for host url
    #' @return data.table of bioactivity estimates for all spid-aeid
    #'  combinations (SC2 table)
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
    # import sc2 data.table
    sc2 <- tcplPrepOtpt(tcplLoadData(lvl = 2, type = "sc"))
    # annotate concentrations from sc_agg
    agg <- tcplPrepOtpt(
        tcplLoadData(
            lvl = "agg",
            type = "sc",
            fld = "s2id",
            val = list(unique(sc2$s2id))
        )
    )
    sc2$logc <- agg$logc[match(sc2$s2id, agg$s2id)]
    # save table
    save(sc2, file = filepath_save)
    return(sc2)
}

filterForTarget <- function(mc5_refchems, ace, class_name) {
    #' Match invitrodb endpoints with refchemdb cluster annotations
    #' 
    #' @param mc5_refchems tibble | table of invitrodb concentration-response
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
    mc5_refs <- mc5_refchems %>%
        filter(
            aeid %in% ace_all$aeid &
            !grepl("CEETOX", aenm) &
            !grepl("ATG_(.+)_TRANS_dn", aenm) &
            !grepl("NVS_ADME", aenm) &
            !grepl("ERF", aenm)
        )
    return(mc5_refs)
}

evalConfActivity <- function(chems_combined) {
    #' Evaluate classification accuracy of invitrodb endpoints versus refchemdb
    #' 
    #' @param chems_combined tibble | table of invitrodb mc5/sc2 measurements as
    #'  output during selectHTSEndpoints()
    #' @return tibble summarizing numbers of chemicals with positive/negative
    #'  endpoint outcomes (`bioactive`) for in-class/ out-of-class annotated
    #'  chemicals (`in_class`)
    #' @example
    #' @export
    # annotate with bioactivity classification based on mc/sc criteria
    chems_combined <- chems_combined %>%
        mutate(bioactive = (hitc == 1 & use.me == 1) | hitc_sc2 == 1)

    # remove duplicate chemical-endpoint pairs (for sc data, keeping highest resp) +
    # count bioactivity classifications for in-class/out-of-class refs
    refs_summary <- chems_combined %>%
        arrange(dsstox_substance_id, aenm, ref_class, desc(max_med_sc2)) %>%
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

evalConfPotency <- function(chems_combined) {
    #' Evaluate classification accuracy of invitrodb endpoints versus refchemdb
    #' 
    #' @param chems_combined tibble | table of invitrodb mc5/sc2 measurements as
    #'  output during selectHTSEndpoints()
    #' @return tibble summarizing numbers of chemicals with positive/negative
    #'  endpoint outcomes (`bioactive`) for in-class/ out-of-class annotated
    #'  chemicals (`in_class`)
    #' @example
    #' @export
    # annotate with bioactivity classification based on mc/sc criteria
    chems_combined <- chems_combined %>%
        mutate(bioactive = (hitc == 1 & use.me == 1) | hitc_sc2 == 1)
    # determine aeids to test
    # (3+ positive and negative annotated chemicals bioactive in aeid)
    aenms_test <- chems_combined %>%
        arrange(dsstox_substance_id, aenm, ref_class, desc(max_med_sc2)) %>%
        distinct(dsstox_substance_id, aenm, ref_class, .keep_all = TRUE) %>%
        filter(bioactive) %>%
        count(aenm, in_class) %>%
        group_by(aenm) %>%
        mutate(
            n_classes = n(),
            test = n_classes == 2 & !any(n < 3)
        ) %>%
        filter(test)

    # run wilcoxon rank-sum tests for log10(ACC) values for each aeid
    ranksum <- chems_combined %>%
        arrange(dsstox_substance_id, aenm, ref_class, desc(max_med_sc2)) %>%
        distinct(dsstox_substance_id, aenm, ref_class, .keep_all = TRUE) %>%
        filter(bioactive & aenm %in% unique(aenms_test$aenm)) %>%
        group_by(aenm) %>%
        nest() %>%
        mutate(
            wilcox = purrr::map(data, function(x) wilcox.test(
                formula = modl_acc ~ in_class, data = x
            )),
            wilcox_df = purrr::map(wilcox, broom::tidy)
        ) %>%
        unnest(wilcox_df) %>%
        mutate(fdr = p.adjust(p.value, method = "fdr"), keep = fdr <= 0.05)

    return(ranksum)
}

selectEndpoints <- function(conf_activity, conf_potency) {
    #' Evaluate classification accuracy of invitrodb endpoints versus refchemdb
    #' 
    #' @param conf_activity tibble | table of classification summary results as
    #'  output by evalConfActivity()
    #' @param conf_potency tibble | table of potency statistical results as
    #'  output by evalConfPotency()
    #' @return vector of aenm values to remove
    #' @example
    #' @export
    # determine whether endpoints pass criteria for
    # bioactivity (count) and potency (ranksum)
    conf_summary <- conf_activity %>%
        group_by(aenm, assay_target, in_class) %>%
        mutate(
            count_total = sum(count),
            count_prop = count / count_total,
            ranksum_keep = conf_potency$keep[match(aenm, conf_potency$aenm)],
            keep = case_when(
                !is.na(ranksum_keep) ~ ranksum_keep,
                is.na(ranksum_keep) & in_class ~ count_prop > 0.5,
                is.na(ranksum_keep) & !in_class ~ count_prop > 0.75,
            )
        ) %>%
        filter(in_class == bioactive) %>%
        group_by(aenm, assay_target) %>%
        mutate(keep_all = all(keep))

    endpoints_remove <- conf_summary %>%
        filter(!is.na(count_prop) & !keep_all) %>%
        pull(aenm) %>%
        unique()
    return(endpoints_remove)
}

calcMode <- function(data) {
    #' Helper function to calculate 1st statistical mode
    #' 
    #' @param data numeric | vector of values to calculate mode from
    #' @return numeric estimate of 1st statistical mode
    #' @export
    if (length(data) < 10) {
        return(2)
    } else {
        # estimate probability density function
    dens <- density(data)
    md <- dens$x[which.max(dens$y)]
    }
    if (length(md) == 1) {
        return(md)
    } else {
        return(min(md))
    }
}

calcBurstMC5 <- function(mc5_full) {
    #' Estimate non-selective points of departure from ToxCast endpoints
    #' 
    #' For each dtxsid, summary statistics for each profile of ACCs are computed
    #' for use in determining Tier 2 assessment outcomes. Note that chemicals
    #' with <10 bioactive endpoints (hitcall == 1 & <3 flags)
    #' are delineated here calculation of medians/modes, as bioactivity among
    #' any endpoint for these chemicals would be considered "selective", and
    #' therefore modes of signature-level bioactivity are set to a default
    #' inactive dose (2-log10(uM)).
    #' 
    #' @param mc5_full data.table | table of table of invitrodb concentration
    #'  -response estimates as output during selectHTSEndpoints().
    #' @return data.table of summary points of departure for each spid,
    #'  including (1) median log10(ACC), (2) 5th percentile log10(ACC), (3)
    #'  absolute 5th log10(ACC), and (4) statistical mode of log10(ACC).
    #' @example
    #' @return
    # compile burst estimates for full screen
    mc5_burst <- mc5_full %>%
        filter(use.me == 1) %>%
        group_by(dsstox_substance_id) %>%
        mutate(
            n_bioactive = n(),
            specific_crit = case_when(
                        n_bioactive < 10 ~ "n<10",
                        TRUE ~ "median_bioactive"
            )
        ) %>%
        summarise(
            across(
                .cols = c(
                    spid, casn, chnm, n_bioactive, specific_crit
                ),
                .fns = first
            ),
            modl_acc_med = median(modl_acc, na.rm = TRUE),
            modl_acc_5 = quantile(modl_acc, probs = c(0.05), na.rm = TRUE),
            modl_acc_abs5 = nth(modl_acc, 5, order_by = modl_acc),
            modl_acc_mode = calcMode(modl_acc),
            modl_acc_sd = sd(modl_acc, na.rm = TRUE)
        )
    return(mc5_burst)
}