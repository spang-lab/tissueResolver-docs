source("R/simulation/helpers.R")
# fits bulks using tissueResolver
fit_bulks_tr <- function(bulks, sccounts) {
    scmat <- sccounts %>%
        dplyr::select(starts_with("ENSG")) %>%
        as.matrix()
    rownames(scmat) <- sccounts %>% pull(cell_id)
    bmat <- bulks %>%
        dplyr::select(-bulk_id) %>%
        as.matrix()
    rownames(bmat) <- bulks %>% pull(bulk_id)
    # call here the parallel version of fit tissue
    # for efficient computation of bulks in parallel on multiple cores
    tissuemodel <- tissueResolver::fit_tissue(
        t(bmat),
        t(scmat),
        maxit = 1e4,
        bootstrap = FALSE, # no bootstrapping
        ncores = 50
    )
    return(tissuemodel)
}

# uses tr's tissuemodel to infer cell type specific proportions
benchmark_tr <- function(bulks, sc_library, mapping, test_samples, irun) {
    # Deconvolute only test data
    bulks <- bulks %>% filter(bulk_id %in% test_samples)
    # determine cell proportions
    tissuemodel <- fit_bulks_tr(bulks, sc_library)
    mapping_tr <- sc_library %>% dplyr::select(cell_id, celltype)
    weights <- tissueResolver::cell_proportions(
        tissuemodel = tissuemodel,
        mapping = mapping_tr,
        by = "celltype",
        weight = "weight",
        bulk_id = "bulk_id",
        cell_id = "cell_id"
    )
    cell_ids <- sc_library %>% pull(cell_id)
    scmat <- sc_library %>%
        dplyr::select(starts_with("ENSG")) %>%
        as.matrix()
    rownames(scmat) <- cell_ids
    smat <- t(scmat)
    csre <- tissueResolver::specific_expression_regulation(tissuemodel,
                                    smat,
                                    mapping_tr,
                                    by="celltype",
                                    weight="weight",
                                    bulk_id = "bulk_id",
                                    cell_id="cell_id",
                                    compute_total = TRUE
                                  )
    cse <- csre %>% dplyr::select(bulk_id, celltype, gene, expression) %>%
                dplyr::rename(cse = expression) %>%
                filter(celltype!='total_explained')

    cell_proportions <- weights %>% 
        group_by(bulk_id) %>%
        mutate(bulk_sum = sum(weight)) %>%
        ungroup() %>% 
        mutate(prp = weight/bulk_sum) %>%
        dplyr::select(bulk_id, celltype, prp)

    # just reshape to fit the other algos
    cse <- cse %>%
            pivot_wider(names_from = "celltype", values_from = "cse")
    return(list(cse = cse, cell_proportions = cell_proportions))
}