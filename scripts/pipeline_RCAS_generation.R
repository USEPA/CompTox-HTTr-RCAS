# NAM-integration-pilot workflow 2
# Generation of Reference Chemical-Associated Chemicals by univariate testing

library(tidyr)
library(dplyr)


analyzeHTTrANOVA <- function(
  filepath.data, filepath.ref, col.httr.type,
  col.ref.class = "cluster",
  col.ref.dtxsid = "dsstox_substance_id",
  col.ref.name = "name", ref.class.clean = TRUE,
  ref.class.curate = TRUE,
  coff.hitc = 0.9, coff.toc = 1.5,
  col.httr.metric = "bmd_log",
  metric.fill = 3, class.count = 3,
  p.adjust.method = "fdr",
  coff.aov.p = 0.05, coff.posthoc.p = 0.05
) {
  #' runtime function for ANOVA gene- or signature-level analysis
  #' @param filepath.data character | path of .RData file containing HTTr
  #'   concentration-response data
  #' @param filepath.ref character | path of .RData file containing reference
  #'   chemicals + classes
  #' @param col.httr.type character | name of column designating type of
  #'   variable. Must be one of c("gene", "signature")
  #' @param col.ref.class character | name of column designating reference
  #'   classes
  #' @param col.ref.dtxsid character | name of column designating chemical
  #'   DTXSIDs
  #' @param col.ref.name character | name of column designating chemical names
  #' @param ref.class.clean logical | if TRUE, performs curated cleaning of
  #'   class target names (i.e. PPARa -> PPAR)
  #' @param ref.class.curate logical | if TRUE, performs further curation of
  #'   selected chemical-class annotations (i.e. Budesonide -> NR3C1 Agonist)
  #' @param coff.hitc double | value of lower bound for `hitcall` variable
  #'   (default = 0.9)
  #' @param coff.toc double | value of lower bound for 'top_over_cutoff'
  #'   variable (default = 1.5)
  #' @param col.httr..metric character | name of column designating independent
  #'   variable to test. Column must be of type `double`.
  #' @param metric.fill double | default value to assign to inactive genes/
  #'   signatures. Value should be chosen based on the metric used.
  #' @param class.count int | minimum number of chemicals per reference class
  #'   used to determine classes to keep
  #' @param p.adjust.method character | name of method to pass to p.adjust()
  #'   for ANOVA p-value multiple comparison correction
  #' @param coff.aov.p double | value of cutoff to select significant
  #'   genes/signatures for posthoc analysis
  #' @param coff.posthoc.p double | value of cutoff to select significant
  #'   class-pairs after posthoc analysis
  #' @return list of objects: [ref], [httr.wide], [aov.return], [posthoc.return]
  #' @example httr.anova <- analyzeHTTrANOVA(filepath.data, filepath.ref, "gene", "cluster_target")
  #' @export
  # load httr/ref data
  httr <- loadData(filepath.data)
  ref <- loadRefList(
    filepath.ref,
    col.dtxsid = col.ref.dtxsid,
    col.name = col.ref.name,
    col.class = col.ref.class,
    class.clean = ref.class.clean,
    class.curate = ref.class.curate
  )

  # filter httr data for refs/active genes
  httr.filtered <- filterData(
    httr,
    ref,
    col.type = col.httr.type,
    coff.hitc = coff.hitc, coff.toc = coff.toc
  )

  # widen httr data + remove underrepresented ref classes
  httr.wide <- widenData(
    httr.filtered,
    ref,
    col.type = col.httr.type,
    col.metric = col.httr.metric,
    metric.fill = metric.fill,
    class.count = class.count
  )

  # conduct ANOVAs for ref classes across all genes/signatures
  aov.return <- calcANOVA(httr.wide, p.adjust.method = p.adjust.method)

  # conduct Tukey's HSD posthoc for ref classes across significant
  # genes/signatures
  posthoc.return <- calcPosthoc(
    httr.wide,
    aov.return$aov.p.adjust,
    aov.p.coff = coff.aov.p,
    posthoc.p.coff = coff.posthoc.p
  )

  bmd <- filter(
    httr,
    gene %in% posthoc.return$posthoc.vars &
    dtxsid %in% httr.wide$dtxsid
  ) %>%
    mutate(
      ref_class = httr.wide$ref_class[match(dtxsid, httr.wide$dtxsid)]
    )

  # return object list
  return(list(
    ref = ref,
    httr.wide = httr.wide,
    aov.estimate = aov.return,
    posthoc.estimate = posthoc.return,
    bmd = bmd
  ))
}


multiFactorHTTrAnova <- function(filepath.data, filepath.ref, col.httr.type,
                                 col.ref.class = "cluster",
                                 col.ref.dtxsid = "dsstox_substance_id",
                                 col.ref.name = "name", ref.class.clean = TRUE,
                                 ref.class.curate = FALSE,
                                 coff.hitc = 0.9, coff.toc = 1.5,
                                 col.httr.metric = "bmd_log",
                                 metric.fill = 3, class.count = 3,
                                 p.adjust.method = "fdr",
                                 coff.aov.p = 0.05) {
# load httr/ref data
  httr <- loadData(filepath.data)
  ref <- loadRefList(filepath.ref,
    col.dtxsid = col.ref.dtxsid,
    col.name = col.ref.name,
    col.class = col.ref.class,
    class.clean = ref.class.clean,
    class.curate = ref.class.curate
  )

  # filter httr data for refs/active genes
  httr.filtered <- filterData(httr, ref,
    col.type = col.httr.type,
    coff.hitc = coff.hitc, coff.toc = coff.toc
  )

  # conduct ANOVAs for ref classes across all genes/signatures
  aov.return <- calcANOVAFactorial(httr.filtered, ref,
                                   col.type = col.httr.type,
                                   col.metric = col.httr.metric,
                                   metric.fill = metric.fill,
                                   class.count = class.count,
                                   p.adjust.method = p.adjust.method)
  # unpack httr.wide object (for alignment with `analyzeHTTrANOVA` output)
  httr.wide <- aov.return$httr.wide
  aov.return$httr.wide <- NULL
  
  # return object list
  return(list(
    ref = ref,
    httr.wide = httr.wide,
    aov.estimate = aov.return
  ))
}


loadData <- function(filepath.data) {
  #' utility function for loading HTTr concentration-response curve fit data
  #' @param filepath.data character | path of .RData file containing HTTr
  #'   concentration-response data
  #' @return tibble of concentration-response data witth additional log(BMD)
  #'   column `bmd_log`
  #' @example httr <- loadData(filepath)
  #' @export
  # check for file compatibility
  if (!file.exists(filepath.data)) {
    stop("Error: Problem with `filepath.data`. \n `filepath.data` not found")
  } else if (tools::file_ext(filepath.data) != "RData") {
    stop("Error: `filepath.data` must be an RData file")
  }
  # import + harmonize object name + add bmd_log
  load(filepath.data)
  httr <- get(ls()[ls() != "filepath.data"])
  httr <- httr %>%
    as_tibble() %>%
    mutate(bmd_log = case_when(!is.na(bmd) ~ log10(bmd)))
  return(httr)
}


cleanRefLabels <- function(ref, col.class) {
  #' utility function to clean up reference class labels (manually)
  #' @param ref tibble | reference chemicals with class assignments as passed
  #'   from loadRefList()
  #' @return tibble of ref with cleaned up labels
  #' @export
  # get list of references only + perform cleanup
  ref.list <- ref[[col.class]]
  ref.list <- gsub("ESR[1-2]_", "ER_", ref.list)
  ref.list <- gsub("PTGS[1-2]_", "COX_", ref.list)
  ref.list <- gsub("PPAR[A-G]_", "PPAR_", ref.list)
  ref.list <- gsub("RAR[A-G]_", "RAR_", ref.list)
  ref.list <- gsub("RXR[A-G]_", "RXR_", ref.list)
  ref.list <- gsub("KCNH2_", "hERG_", ref.list)
  ref.list[grepl("[A-B]CHE_Negative", ref.list)] <- "AChE_Inhibitor"
  ref.list[grepl("CA[3-9]|CA1[1-9]&Negative", ref.list)] <- "PanCA_Inhibitor"
  ref.list[grepl("CA[1-2]_Negative", ref.list)] <- "CA1_Inhibitor"
  ref.list[grepl("CA[1-2]_Positive", ref.list)] <- "CA1_Agonist"
  ref.list[grepl("CHRN[A-B][1-7]_Positive", ref.list)] <- "nAChR_Agonist"
  ref.list[grepl("CHRM[1-5]_Negative", ref.list)] <- "mAChR_Inhibitor"
  ref.list[grepl("CYP19A1_Negative", ref.list)] <- "CYP19A1_Inhibitor"
  ref.list[grepl("HDAC[1-3]_Negative", ref.list)] <- "HDAC_Inhibitor"
  ref.list[grepl("HRH[1-4]_Positive", ref.list)] <- "HRH_Agonist"
  ref.list[grepl("HRH[1-4]_Negative", ref.list)] <- "HRH_Inhibitor"
  ref.list[grepl("PDE[1-4][A-D]_Negative", ref.list)] <- "PDE_Inhibitor"
  ref.list[grepl("MAO[A-B]_Negative", ref.list)] <- "MAO_Inhibitor"
  ref.list[grepl("TUB[A-B].*_Negative", ref.list)] <- "TUB_Inhibitor"
  ref.list[grepl("NR3C[1-2]_Negative", ref.list)] <- "NR3C1_Inhibitor"
  ref.list[grepl("NR1I2_Positive", ref.list)] <- "PXR_Agonist"
  ref.list[grepl("NR1I3_Positive", ref.list)] <- "CAR_Agonist"
  ref.list[grepl("NR1H4_Positive", ref.list)] <- "FXR_Agonist"
  ref.list[grepl("RAR|RXR&Positive&Negative", ref.list)] <- "RAR_RXR_Agonist"
  ref.list <- gsub("_Positive", "_Agonist", ref.list)
  ref.list <- gsub("_Negative", "_Inhibitor", ref.list)

  # remove duplicates (within entries) after cleaning
  removeDupNames <- function(entry) {
    entry.unique <- unique(unlist(strsplit(entry, "\\|")))
    if (length(entry.unique) > 1) {
      entry.return <- paste(entry.unique, collapse = "_")
    } else {
      entry.return <- entry.unique
    }
    return(entry.return)
  }
  
  ref.list <- sapply(ref.list, removeDupNames)

  # integrate clean list into full df
  ref.clean <- ref
  ref.clean[[col.class]] <- unname(ref.list)
  return(ref.clean)
}

curateRefLabels <- function(ref, col.class, curated.list = NULL) {
  #' utility function to manually adjust classes in cases not accounted for by
  #' RefChemDB
  #' @param ref tibble | reference chemicals with class assignments as passed
  #'   from loadRefList()
  #' @param curated.list list of 2 | list of vectors detailing (1) chemical
  #'   names, and (2) reference classes to assign
  #' @return tibble of ref with cleaned up labels
  #' @export
  # check curated.list for same length of component vectors
  # use default list of curated chemical-class annotations if not specified
  if (is.null(curated.list)) {
    curated.list <- list(
      name = c(
        "Budesonide", "Cortisone acetate", "Desoximetasone", "Paramethasone"
      ),
      class = c(
        "NR3C1 Agonist", "NR3C1 Agonist", "NR3C1 Agonist", "NR3C1 Agonist"
      )
    )
  }
  # replace current class annotation with curated list
  ref.curated <- ref
  for (term in seq_len(length(curated.list$name))) {
    ref.curated[[col.class]][ref.curated$name == curated.list$name[term]] <- curated.list$class[term]
  }
  return(ref.curated)
}

loadRefList <- function(filepath.ref,
                        col.dtxsid = "dsstox_substance_id",
                        col.name = "name",
                        col.class = "cluster",
                        class.clean = FALSE,
                        class.curate = FALSE) {
  #' utility function for loading reference chemical sets/clusters
  #' @param filepath.ref character | path of .RData file containing reference
  #'   chemicals + classes
  #' @param col.dtxsid character | name of column designating chemical DTXSIDs
  #' @param col.name character | name of column designating chemical names
  #' @param col.class character | name of column designating reference classes
  #' @param class.clean logical | flag to implement cleanup of class labels
  #' @param class.curate logical | flag to implement curation of select chemical
  #'   class annotations
  #' @return tibble of reference chemicals with class assignments
  #' @example ref <- loadRefList(filepath)
  #' @export
  # check for file compatibility
  if (!file.exists(filepath.ref)) {
    stop("Error: Problem with `filepath.ref`.\n `filepath.ref` not found")
  } else if (tools::file_ext(filepath.ref) != "RData") {
    stop("Error: `filepath.ref` must be an RData file")
  }
  # import + harmonize object/column names
  load(filepath.ref)
  ref <- get(ls()[!ls() %in% c(
    "filepath.ref", "col.dtxsid", "col.name",
    "col.class", "class.clean", "class.curate"
  )])
  if (!"tbl_df" %in% class(ref)) {
    ref <- tibble::as_tibble(ref)
  } else if (class.clean & class.curate) {
    ref <- cleanRefLabels(ref, col.class)
    ref <- curateRefLabels(ref, col.class)
  } else if (class.clean & !class.curate) {
    ref <- cleanRefLabels(ref, col.class)
  }
  ref.export <- dplyr::select(
    ref,
    tidyselect::all_of(c(col.dtxsid, col.name, col.class))
  )
  colnames(ref.export) <- c("dtxsid", "name", "ref_class")
  return(ref.export)
}


filterData <- function(httr, ref, col.type,
                       ref.count = 3, coff.hitc = 0.9, coff.toc = 1.5) {
  #' utility function for applying filtering strategy to HTTr data, based on
  #' reference list + active gene representation
  #' @param httr tibble | HTTr gene/signature concentration-response curve fits
  #'   as passed from loadData()
  #' @param ref tibble | reference chemical/class data as passed from
  #'   loadRefList()
  #' @param col.type character | name of column designating type of variable.
  #'   Must be one of c("gene", "signature")
  #' @param ref.count int | minimum number of reference chemicals in a class
  #'   needed to keep that class
  #' @param coff.hitc double | value of lower bound for `hitcall` variable
  #'   (default = 0.9)
  #' @param coff.toc double | value of lower bound for 'top_over_cutoff'
  #'   variable (default = 1.5)
  #' @return tibble of HTTr concentration-response data for reference chemicals
  #'   and active genes/signatures only
  #' @example httr.filtered <- filterData(httr, ref, "gene")
  #' @export
  # validate col.type argument
  if (!col.type %in% c("gene", "signature")) {
    stop("Error: Problem with `col.type`.\n `col.type` must be one of c(\"gene\", \"signature\")")
  } else if (!col.type %in% colnames(httr)) {
    stop("Error: Problem with `httr`.\n `httr` must contain a column matching `col.type`")
  }

  # find references in HTTr data (via dtxsid) + count class representation
  ref.in_httr <- ref %>%
    filter(dtxsid %in% httr$dtxsid) %>%
    group_by(ref_class) %>%
    mutate(ref_class_count = n()) %>%
    ungroup()

  # remove chemicals that either has no assigned class (is.na) + classes with
  # <3 chemicals
  ref.filtered <- ref.in_httr %>%
    filter(!is.na(ref_class) & ref_class_count >= ref.count)

  # subset httr for resulting reference substances
  if (col.type == "signature") {
    httr.filtered <- httr %>%
      filter(dtxsid %in% unique(ref.filtered$dtxsid) &
        !grepl("Random_", signature)) %>%
      distinct(dtxsid, signature, .keep_all = TRUE) %>%
      filter(hitcall > 0.9 & top_over_cutoff > 1.5)
  } else {
    httr.filtered <- httr %>%
      filter(dtxsid %in% unique(ref.filtered$dtxsid)) %>%
      distinct(dtxsid, gene, .keep_all = TRUE) %>%
      filter(hitcall > coff.hitc & top_over_cutoff > coff.toc)
  }
  return(httr.filtered)
}


widenData <- function(httr.filtered, ref, col.type,
                      col.metric = "bmd_log",
                      metric.fill = 3,
                      class.count = 3) {
  #' utility function for shaping HTTr data to wide form for ANOVA
  #' @param httr.filtered tibble | reference-only active concentration-response
  #'   curve fits as passed from filterData()
  #' @param ref tibble | reference chemical/class data as passed from
  #'   loadRefList()
  #' @param col.type character | name of column designating type of variable.
  #'   Must be one of c("gene", "signature")
  #' @param col.metric character | name of column designating independent
  #'   variable to test. Column must be of type `double`.
  #' @param metric.fill double | default value to assign to inactive genes/
  #'   signatures. Value should be chosen based on the metric used.
  #' @param class.count int | minimum number of chemicals per reference class
  #'   used to determine classes to keep
  #' @return wide-form tibble of HTTr concentration-response metrics,
  #'   structured as chemicals x genes/signatures
  #' @example httr.wide <- widenData(httr.filtered, "gene")
  #' @export
  # validate col.type argument
  if (!col.type %in% c("gene", "signature")) {
    stop("Error: Problem with `col.type`.\n `col.type` must be one of c(\"gene\", \"signature\")")
  } else if (!col.type %in% colnames(httr.filtered)) {
    stop("Error: Problem with `httr.filtered`.\n `httr` must contain a column matching `col.type`")
  }
  # validate col.metric argument
  if (!col.metric %in% colnames(httr.filtered)) {
    stop("Error: Problem with `col.metric`.\n `col.metric` must be a column in `httr.filtered`")
  } else if (!is.numeric(httr.filtered[[col.metric]])) {
    stop("Error: Problem with `col.metric`.\n `col.metric` must be a numeric column")
  }

  # widen httr data + remove dtxsids for classes with <class.count chemicals
  # after widening
  httr.wide <- httr.filtered %>%
    pivot_wider(
      id_cols = c(name, dtxsid),
      names_from = as.name(col.type),
      values_from = as.name(col.metric)
    ) %>%
    replace(is.na(.), metric.fill) %>%
    arrange(dtxsid) %>%
    mutate(ref_class = ref$ref_class[match(dtxsid, ref$dtxsid)]) %>%
    group_by(ref_class) %>%
    mutate(ref_class_count = n()) %>%
    filter(ref_class_count >= class.count) %>%
    select(-ref_class_count) %>%
    ungroup()

  # replace special characters from column names
  colnames(httr.wide) <- gsub("/", "_", colnames(httr.wide))

  return(httr.wide)
}


calcANOVA <- function(httr.wide, p.adjust.method = "fdr") {
  #' utility function for calculating one-way ANOVA test statistics for
  #' reference class variable across HTTr genes/signatures
  #' @param httr.wide tibble | wide-form, reference-only active concentration-
  #'   response curve fits as passed from widenData()
  #' @param p.adjust.method character | name of method to pass to p.adjust()
  #'   for ANOVA p-value multiple comparison correction
  #' @return list of vectors: [aov.vars] variables tested via ANOVA, [aov.p]
  #'   raw p-values of ANOVA F-tests, and [aov.p.adjust] multiple
  #'   comparison-adjusted p-values of ANOVA F-tests
  #' @example aov.return <- calcANOVA(httr.wide)
  #' @export
  # get list of variables for iteration + coerce names for use in formulas
  require(dplyr)
  aov.vars <- colnames(dplyr::select(httr.wide, where(is.numeric)))
  aov.vars <- aov.vars[!grepl("-", aov.vars)]
  aov.vars[grepl(" ", aov.vars)] <- paste0(
    "`", aov.vars[grepl(" ", aov.vars)], "`"
  )

  # perform ANOVA for each variable + store variables/p-values
  aov.p <- NULL
  aov.p.vars <- NULL
  for (variable in aov.vars) {
    aov.formula <- as.formula(paste(variable, "~", "ref_class"))
    aov.estimate <- summary(aov(aov.formula, data = httr.wide))
    aov.p <- c(aov.p, aov.estimate[[1]][["Pr(>F)"]][1])
    aov.p.vars <- c(aov.p.vars, variable)
  }
  names(aov.p) <- aov.p.vars

  # adjust p-value vector via specified method
  aov.p.adjust <- p.adjust(aov.p, method = p.adjust.method)
  return(list(
    aov.vars = aov.vars,
    aov.p = aov.p,
    aov.p.adjust = aov.p.adjust
  ))
}

calcANOVAFactorial <- function(httr.filtered, ref, col.type,
                               col.metric = "bmd_log",
                               metric.fill = 3,
                               class.count = 3,
                               p.adjust.method = "fdr") {
  #' utility function for calculating test statistics for reference classes
  #' across HTTr genes/signatures, using a multi-factor design
  #' @param httr.filtered tibble | reference-only active concentration-response
  #'   curve fits as passed from filterData()
  #' @param ref tibble | reference chemical/class data as passed from
  #'   loadRefList()
  #' @param col.type character | name of column designating type of variable.
  #'   Must be one of c("gene", "signature")
  #' @param col.metric character | name of column designating independent
  #'   variable to test. Column must be of type `double`.
  #' @param metric.fill double | default value to assign to inactive genes/
  #'   signatures. Value should be chosen based on the metric used.
  #' @param class.count int | minimum number of chemicals per reference class
  #'   used to determine classes to keep
  #' @param p.adjust.method character | name of method to pass to p.adjust()
  #'   for ANOVA p-value multiple comparison correction
  #' @return list of objects: [aov.vars] variables tested via ANOVA, [aov.p]
  #'   raw p-values of ANOVA F-tests, and [aov.p.adjust] multiple
  #'   comparison-adjusted p-values of ANOVA F-tests
  #' @example aov.return <- calcANOVA(httr.wide)
  #' @export
  require(dplyr)
  # find references in HTTr data (via dtxsid) + count class representation
  ref.in_httr <- ref %>%
    filter(dtxsid %in% httr.filtered$dtxsid) %>%
    group_by(ref_class) %>%
    mutate(ref_class_count = n()) %>%
    ungroup()

  # remove chemicals that either has no assigned class (is.na) + classes with
  # <3 chemicals
  ref.filtered <- ref.in_httr %>%
    filter(!is.na(ref_class) & ref_class_count >= class.count)

  # get design matrix from ref.filtered + remove classes with few chemicals
  # (< class.count)
  ref.design <- ref.filtered %>%
    mutate(measured = 1) %>%
    pivot_wider(
      id_cols = c(dtxsid, name),
      names_from = ref_class,
      values_from = measured,
      values_fill = 0
    )

  # widen httr data + concatenate design matrix
  httr.wide <- httr.filtered %>%
    pivot_wider(
      id_cols = c(name, dtxsid),
      names_from = as.name(col.type),
      values_from = as.name(col.metric)
    ) %>%
    replace(is.na(.), metric.fill) %>%
    arrange(dtxsid) %>%
    left_join(., ref.design, by = c("dtxsid", "name"))

  # get list of variables for iteration + coerce names for use in formulas
  aov.vars <- colnames(select(httr.wide, where(is.numeric)))
  aov.vars <- aov.vars[!grepl("Agonist|Inhibitor", aov.vars)]
  aov.vars <- aov.vars[!grepl("-", aov.vars)]
  aov.vars[grepl(" ", aov.vars)] <- paste0(
    "`", aov.vars[grepl(" ", aov.vars)], "`"
  )
  aov.factors <- colnames(httr.wide)[
    grepl("Agonist|Inhibitor", colnames(httr.wide))
  ]
  # perform ANOVA for each variable + store variables/p-values
  aov.p <- list()
  aov.p.vars <- NULL
  for (variable in aov.vars) {
    aov.formula <- as.formula(
      paste(variable, "~", paste(aov.factors, collapse = " + "))
    )
    aov.estimate <- summary(aov(aov.formula, data = httr.wide))
    aov.p[[variable]] <- aov.estimate[[1]][["Pr(>F)"]]
    names(aov.p[[variable]]) <- trimws(rownames(aov.estimate[[1]]))
    aov.p.vars <- c(aov.p.vars, variable)
  }
  aov.df <- as.data.frame(aov.p)

  # adjust p-value vector via specified method + coerce back to df
  aov.df.adjust <- p.adjust(
    as.vector(as.matrix(aov.df)), method = p.adjust.method
  )
  aov.df.adjust <- as.data.frame(matrix(
    aov.df.adjust, nrow = nrow(aov.df), ncol = ncol(aov.df)
  ))
  rownames(aov.df.adjust) <- rownames(aov.df)
  colnames(aov.df.adjust) <- gsub("`", "", aov.p.vars)
  return(list(
    httr.wide = httr.wide,
    aov.vars = aov.p.vars,
    aov.p = aov.df,
    aov.p.adjust = aov.df.adjust
  ))
}

calcPosthoc <- function(httr.wide, aov.p.adjust,
                        aov.p.coff = 0.05,
                        posthoc.p.coff = 0.05) {
  #' utility function for conducting multiple comparisons for genes/signatures
  #' found significant via ANOVA. Tests are conducted via Tukey's HSD
  #' @param httr.wide tibble | wide-form, reference-only active concentration-
  #'   response curve fits as passed from widenData()
  #' @param aov.p.adjust named vector | vector of multiple comparison-adjusted
  #'   p-values (names from colnames(httr.wide)) passed from calcANOVA()
  #' @param aov.p.coff double | value of cutoff to select significant
  #'   genes/signatures for posthoc analysis
  #' @param posthoc.p.coff double | value of cutoff to select significant
  #'   class-pairs after posthoc analysis
  #' @return list of objects: [posthoc.vars] vector of variables tested,
  #'   [posthoc.all] tibble of all class comparisons for all variables,
  #'   [posthoc.sig] tibble of significant class comparisons for all variables
  #' @example posthoc.return <- calcPosthoc(httr.wide, aov.return$aov.p.adjust)
  #' @export
  # get list of variables with aov.p.adjust <= aov.p.coff
  posthoc.vars <- names(aov.p.adjust[which(aov.p.adjust <= aov.p.coff)])

  # perform multiple comparison (TukeyHSD) for each variable + concatenate dfs
  posthoc.all <- NULL
  for (variable in posthoc.vars) {
    aov.formula <- as.formula(paste(variable, "~", "ref_class"))
    aov.estimate <- aov(aov.formula, data = httr.wide)
    posthoc.estimate <- TukeyHSD(aov.estimate)
    posthoc.df <- posthoc.estimate$ref_class %>%
      as_tibble(rownames = "levels") %>%
      separate(levels, into = c("class_1", "class_2"), sep = "-") %>%
      mutate(variable = variable)
    # add reference class means to posthoc df
    variable_adj <- gsub("`", "", variable)
    ref.means <- httr.wide[, c(variable_adj, "ref_class")] %>%
      group_by(ref_class) %>%
      summarise(across(all_of(variable_adj), mean))
    posthoc.df$mean_1 <- ref.means[[variable_adj]][match(posthoc.df$class_1, ref.means$ref_class)]
    posthoc.df$mean_2 <- ref.means[[variable_adj]][match(posthoc.df$class_2, ref.means$ref_class)]
    posthoc.all <- rbind(posthoc.all, posthoc.df)
  }

  # filter for significant comparisons + unique comparisons (via diff direction)
  # negative diff: keep direction + class assignments
  posthoc.sig.neg <- posthoc.all %>%
    filter(`p adj` <= posthoc.p.coff & diff < 0) %>%
    mutate(variable = gsub("`", "", variable)) %>%
    unite(class_both, class_1, class_2, sep = "|", remove = FALSE)
  # positive diff: reverse direction + class assignments (for comparison with
  # negatives)
  posthoc.sig.pos <- posthoc.all %>%
    filter(`p adj` <= posthoc.p.coff & diff > 0) %>%
    mutate(variable = gsub("`", "", variable)) %>%
    unite(class_both, class_1, class_2, sep = "|", remove = FALSE) %>%
    relocate(
      class_1 = class_2,
      class_2 = class_1,
      lwr = upr,
      upr = lwr,
      mean_1 = mean_2,
      mean_2 = mean_1
    ) %>%
    mutate(across(c(lwr, upr, diff), ~ -1 * .x))

  posthoc.sig <- rbind(posthoc.sig.neg, posthoc.sig.pos) %>%
    distinct(class_1, class_2, variable, .keep_all = TRUE)

  return(list(
    posthoc.vars = posthoc.vars,
    posthoc.all = posthoc.all,
    posthoc.sig = posthoc.sig
  ))
}