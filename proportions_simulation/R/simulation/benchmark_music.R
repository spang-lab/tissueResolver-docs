# give the bulk expression in the appropriate format
compute_bulk_counts_music <- function(bulk_counts, test_samples) {
    # we deconvolute only the test samples
    bulk_counts <- bulk_counts %>% filter(bulk_id %in% test_samples)
    bulks_music <- bulk_counts %>%
        dplyr::select(-bulk_id) %>%
        as.matrix()
    rownames(bulks_music) <- bulk_counts %>% pull(bulk_id)
    bulks_music <- t(bulks_music)

    return(bulks_music)
}

benchmark_music <- function(sc_library, sclong, bulk_counts, test_samples) {
    # music expects the sc library as a SingleCellExperiment
    # dataset, so prepare that here
    sc_library <- sc_library %>% inner_join(sclong %>% dplyr::select(cell_id, Patient), by = "cell_id")
    celltype <- sc_library %>% pull(celltype)
    cell_id <- sc_library %>% pull(cell_id)
    patient <- sc_library %>% pull(Patient)

    sc_counts <- sc_library %>%
        dplyr::select(starts_with("ENSG")) %>%
        as.matrix()
    gene_names <- colnames(sc_counts)
    rownames(sc_counts) <- cell_id
    sc_counts <- t(sc_counts)

    steen_sce <- SingleCellExperiment(list(counts = sc_counts),
        rowData = DataFrame(gene.name = gene_names),
        colData = DataFrame(sampleID = patient, cellType = celltype)
    )

    bulk_counts <- compute_bulk_counts_music(bulk_counts, test_samples)

    music_output <- music_prop(
        bulk.mtx = bulk_counts, sc.sce = steen_sce,
        clusters = "cellType", samples = "sampleID"
    )

    # perform the music algorithm
    proportions <- music_output$Est.prop.weighted %>% as_tibble() %>%
        cbind(tibble(bulk_id = rownames(music_output$Est.prop.weighted))) %>%
        pivot_longer(!bulk_id, names_to = "celltype", values_to = "prp")

    return(proportions)
}