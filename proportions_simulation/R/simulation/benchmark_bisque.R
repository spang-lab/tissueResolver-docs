prepare_bmind_input <- function(bulk_counts, sc_library) {
    bulks <- bulk_counts %>%
        dplyr::select(-bulk_id) %>%
        # bmind can not deal with genes identical to 0
        dplyr::select(where(~ n_distinct(.) > 1)) %>%
        as.matrix()
    rownames(bulks) <- bulk_counts %>% pull(bulk_id)
    bulks <- t(bulks)

    available_celltypes <- sc_library %>%
        pull(celltype) %>%
        unique()

    sc <- sc_library %>%
        dplyr::select(starts_with("ENSG")) %>%
        as.matrix()
    rownames(sc) <- sc_library %>% pull(cell_id)
    sc <- t(sc)

    sc_meta <- sc_library %>%
        dplyr::select(cell_id, celltype) %>%
        dplyr::rename(subject = cell_id, cell_type = celltype) %>%
        data.frame()
    row.names(sc_meta) <- sc_meta$subject

    return(list(bulk_counts = bulks, sc = sc, sc_meta = sc_meta))
}

benchmark_bisque <- function(sc_library, bulk_counts, test_samples) {
    bulk_counts <- bulk_counts %>% filter(bulk_id %in% test_samples)
    input_bisque <- prepare_bmind_input(bulk_counts, sc_library)
    # infer proportions with Bisque from sc data
    bisque_proportions <- est_frac_sc(
        bulk = input_bisque$bulk_counts,
        frac_method = "Bisque",
        sc_count = input_bisque$sc,
        sc_meta = input_bisque$sc_meta
    )
    proportions <- bisque_proportions %>%
        as_tibble(rownames = "bulk_id") %>%
        pivot_longer(-bulk_id, names_to = "celltype", values_to = "prp") %>%
        mutate(bulk_id = as.integer(bulk_id))
    return(proportions)
}