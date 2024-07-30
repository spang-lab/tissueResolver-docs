# helpers to read and write cibersort input/output

compute_reference_cibersort <- function(sc_library) {
    reference <- tibble()
    reference_tibble <- sc_library %>%
        pivot_longer(cols = starts_with("ENSG"), names_to = "gene", values_to = "count") %>%
        group_by(celltype, gene) %>%
        mutate(cs_expr = mean(count)) %>%
        ungroup() %>%
        dplyr::select(celltype, gene, cs_expr) %>%
        unique()
    reference_tibble <- reference_tibble %>%
        pivot_wider(names_from = gene, values_from = cs_expr)
    reference <- reference_tibble %>%
        dplyr::select(-celltype) %>%
        as.matrix()
    cell_names <- reference_tibble %>% pull(celltype)
    rownames(reference) <- cell_names
    reference <- t(reference)
    gene_symbols <- rownames(reference)
    reference_gene_names <- cbind("gene_name" = gene_symbols, reference)
    reference_matrix <- rbind(c("name", cell_names), reference_gene_names)
    return(reference_matrix)
}

compute_refsample_cibersort <- function(sc_library) {
    celltypes <- sc_library %>% pull(celltype)
    refsample <- sc_library %>%
        dplyr::select(starts_with("ENSG")) %>%
        as.matrix()
    # rownames(refsample) <- celltypes
    refsample <- t(refsample)
    gene_symbols <- rownames(refsample)
    refsample_gene_names <- cbind("gene_name" = gene_symbols, refsample)
    refsample_matrix <- rbind(c("name", celltypes), refsample_gene_names)
    return(refsample_matrix)
}

compute_bulk_counts_cibersort <- function(bulk_counts, test_samples) {
    # we deconvolute only the test samples
    bulk_counts <- bulk_counts %>% filter(bulk_id %in% test_samples)
    bulk_matrix <- bulk_counts %>%
        dplyr::select(-bulk_id) %>%
        as.matrix()
    sample_names <- bulk_counts %>% pull(bulk_id)
    bulk_matrix <- t(bulk_matrix)
    gene_symbols <- rownames(bulk_matrix)
    bulk_counts_gene_names <- cbind("gene_name" = gene_symbols, bulk_matrix)
    bulk_counts_matrix <- rbind(c("name", sample_names), bulk_counts_gene_names)
    return(list(bulk_counts_matrix = bulk_counts_matrix))
}



compute_input_cibersort <- function(sc_library, bulk_counts, test_samples) {
    reference_matrix <- compute_reference_cibersort(sc_library)
    refsample_matrix <- compute_refsample_cibersort(sc_library)
    bulk_counts_list <- compute_bulk_counts_cibersort(bulk_counts, test_samples)
    return(list(
        reference = reference_matrix,
        refsample = refsample_matrix,
        bulk_counts = bulk_counts_list
    ))
}

write_cibersort_input_txt <- function(reference, refsample, bulk_counts, output_dir, irun) {

    write.table(reference,
        file = file.path(output_dir, paste0("sigmatrix_run_", irun, ".txt")),
        quote = FALSE,
        sep = "\t ",
        col.names = FALSE,
        row.names = FALSE
    )
    write.table(refsample,
        file = file.path(output_dir, paste0("refsample_run_", irun, ".txt")),
        quote = FALSE,
        sep = "\t ",
        col.names = FALSE,
        row.names = FALSE
    )
    write.table(bulk_counts,
        file = file.path(output_dir, paste0("bulks_run_", irun, ".txt")),
        quote = FALSE,
        sep = "\t",
        col.names = FALSE,
        row.names = FALSE)
}

read_cibersort_proportions <- function(output_dir, cibersortfile, algoname, irun) {
    cibersort_prp <- read.table(
        file.path(output_dir, paste0(cibersortfile, "_", algoname, "_run_", irun, ".txt")),
        sep = "\t",
        header = TRUE
    )
    cibersort_proportions <- cibersort_prp %>%
        as_tibble() %>%
        dplyr::select(-P.value, -Correlation, -RMSE) %>%
        dplyr::rename(bulk_id = Mixture) %>%
        pivot_longer(!bulk_id, names_to = "celltype", values_to = "prp") %>%
        add_column(run = irun, algo = paste0("cibersort_", algoname)) %>%
        mutate(bulk_id = as.character(bulk_id))
    return(cibersort_proportions)
}

add_cibersort_proportions <- function(output_dir, proportions_all, cibersortfile, nruns) {
    for (irun in 1:nruns) {
        cibersort_proportions <- read_cibersort_proportions(output_dir, cibersortfile, "sc", irun)
        proportions_all <- rbind(proportions_all, cibersort_proportions)
    }
    proportions_all <- proportions_all %>% arrange(run)
    return(proportions_all)
}

