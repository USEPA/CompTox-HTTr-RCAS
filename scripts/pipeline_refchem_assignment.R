# NAM-integration-pilot workflow 1
# Selection of MoA-associated reference chemicals by literature association

loadRefChemDB <- function(
    homedir = getwd(),
    filename = "NIHMS1537541-supplement-Supplement1.xlsx",
    sheet = "S12 Data",
    fileurl = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6784312/bin/NIHMS1537541-supplement-Supplement1.xlsx"
) {
    filedir <- paste0(homedir, "/data")
    if (!filename %in% list.files(filedir)) {
        refchemdb <- openxlsx::read.xlsx(fileurl, sheet = sheet)
    } else {
        refchemdb <- openxlsx::read.xlsx(
            paste0(filedir, filename),
            sheet = sheet
        )
    }
    return(refchemdb)
}

filterSupport <- function(refchemdb, support) {
    require(dplyr)
    require(tidyr)

    # filter refchemdb for threshold support level + distinct mode
    # (i.e. not unspecified)
    filtered <- refchemdb %>%
        filter(support >= support & mode %in% c("Negative", "Positive")) %>%
        distinct(dsstox_substance_id, target, mode, .keep_all = TRUE)

    # for multiple modes annotated for same chemical/target,
    # keep mode with the highest support +
    # concatenate target/mode into single annotation
    filtered <- filtered
        group_by(dsstox_substance_id, target) %>%
        mutate(
            mode_count = n(),
            mode_ratio = case_when(
                mode_count == 2 ~ max(support) / min(support), TRUE ~ 0
            ),
            mode_keep = case_when(
                mode_count == 1 ~ TRUE,
                mode_count == 2 & support == max(support) ~ TRUE,
                TRUE ~ FALSE
            )
        ) %>%
        ungroup() %>%
        unite(target_mode, c(target, mode), remove = FALSE)

    return(filtered)
}

calcJaccard <- function(filtered) {
    # convert filtered refchemdb to wide form (target_mode x chemical), 
    # filling in missing pairs as inactive
    widened <- filtered %>%
        filter(mode_keep == TRUE) %>%
        mutate(active == 1) %>%
        pivot_wider(
            id_cols = target_mode,
            names_from = dsstox_substance_id,
            values_from = active,
            values_fill = 0
        )

    # calculate jaccard distances across chemical pairs
    jaccard <- widened %>%
    select(where(is.numeric)) %>%
    as.matrix() %>%
    philentropy::distance(., method = "jaccard")

    jaccard_obj <- list(widened = widened, jaccard = jaccard)
    return(jaccard_obj)
}

clusterTargets <- function(jaccard_obj, method = "average", dend_height = 0.8) {
    require(tibble)
    require(dplyr)
    widened <- jaccard_obj$widened
    jaccard <- jaccard_obj$jaccard

    # perform hierarchical clustering on target_mode annotations
    hc <- hclust(jaccard, method = method)

    # cut dendrogram to specified height + get clusters
    clusters <- hc %>%
        cutree(h = dend_height) %>%
        as_tibble(rownames = "row") %>%
        mutate(target_mode = widened$target_mode)

    return(clusters)
}

assignChemsToClusters <- function(filtered, clusters, cluster_size = 5) {
    require(tidyr)
    require(dplyr)

    # link target_mode annotations to clusters
    chems_assigned <- filtered %>%
        filter(mode_keep == TRUE) %>%
        mutate(
            cluster = clusters$value[
                match(target_mode, clusters$target_mode)
            ]
        )

    # assign cluster to each chemical:
    ## top cluster by total number of matching target_mode terms
    ## ties: top cluster by sum(support) of remaining clusters
    ## remaining ties: removed from list
    chems_assigned <- chems_assigned %>%
        group_by(dsstox_substance_id, cluster) %>%
        summarise(
            cluster_count = n(),
            support_sum = sum(support),
            .groups = "drop_last"
        ) %>%
        slice_max(cluster_count) %>%    # choose top cluster
        slice_max(support_sum) %>%      # break ties
        mutate(clusters_n = n()) %>%    # remove remaining ties
        filter(clusters_n == 1)

    # remove clusters with <5 chemicals represented
    chems_assigned_rep <- filtered %>%
        filter(dsstox_substance_id %in% chems_assigned$dsstox_substance_id) %>%
        mutate(
            cluster = chems_assigned$cluster[
                match(dsstox_substance_id, chems_assigned$dsstox_substance_id)
            ]
        ) %>%
        distinct(dsstox_substance_id, .keep_all = TRUE) %>%
        group_by(cluster) %>%
        summarise(size = n()) %>%
        mutate(size_total = sum(size)) %>%
        filter(size_total >= cluster_size)


}