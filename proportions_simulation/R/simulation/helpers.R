# computes the proportions of a cell type within each bulk sample,
# by comparing the celltype-wise total gene count to the bulk's total count
compute_proportion <- function(bulk_counts, cse) {
    bulk_total <- compute_bulk_total(bulk_counts)
    proportion <- cse %>%
        group_by(bulk_id, celltype) %>%
        mutate(cse_total = sum(cse)) %>%
        ungroup() %>%
        dplyr::select(-cse) %>%
        dplyr::distinct() %>%
        inner_join(bulk_total, by = "bulk_id") %>%
        mutate(prp = cse_total / bulk_sum) %>%
        # TODO: positive select
        dplyr::select(-gene, -cse_total, -bulk_sum) %>%
        dplyr::distinct()

    return(proportion)
}

# computes the cell type specific expression,
# by weighting the reference according to the weights
# resulting from the deconvolution
compute_cse <- function(weights, reference) {
    cse <- weights %>%
        # celltype is our notion of cell-type
        inner_join(reference, by = "celltype") %>%
        pivot_longer(cols = starts_with("ENSG"), names_to = "gene", values_to = "expression") %>%
        mutate(expression = expression * weight) %>%
        group_by(bulk_id, gene, celltype) %>%
        mutate(cse = sum(expression)) %>%
        ungroup() %>%
        # TODO: maybe positive select better
        dplyr::select(-expression, -weight) %>%
        dplyr::distinct()
    return(cse)
}

# the total bulk count for a single bulk
# is the sum over all gene counts within that bulk
compute_bulk_total <- function(bulk_counts) {
    bulk_total <- bulk_counts %>%
        pivot_longer(cols = starts_with("ENSG"), names_to = "gene", values_to = "expression") %>%
        group_by(bulk_id) %>%
        mutate(bulk_sum = sum(expression)) %>%
        ungroup() %>%
        dplyr::select(bulk_id, bulk_sum) %>%
        dplyr::distinct()

    return(bulk_total)
}

# computes the cell-type specific expression for the
# simulated bulks by using the true cell weights stored
# in the composition
# This will yield our ground truth
compute_bulk_cse <- function(composition, sclong) {
    scwide <- sclong %>% pivot_wider(names_from = gene, values_from = expression, values_fill = 0)
    celltype_comp <- composition %>%
        inner_join(scwide, by = "cell_id") %>%
        pivot_longer(cols = starts_with("ENSG"), names_to = "gene", values_to = "expression") %>%
        # the cse equals each weighted sc profile cumulated over the respective cell type
        mutate(expression = expression * weight) %>%
        group_by(bulk_id, celltype, gene) %>%
        mutate(cse = sum(expression)) %>%
        ungroup()
    cse <- celltype_comp %>%
        dplyr::select(bulk_id, gene, cse, celltype) %>%
        dplyr::distinct()
    return(cse)
}

# compute the actual proportions, by weighting the reference with c
# and determine the proportions between the sum of cse counts
# and the sum of bulk counts
convert_weight_to_proportion <- function(weight, reference, bulk_counts) {
    ref_wide <- reference %>%
        as_tibble(rownames = "gene") %>%
        pivot_longer(-gene, names_to = "celltype", values_to = "expression") %>%
        pivot_wider(names_from = gene, values_from = expression)

    reshaped_weights <- tibble(celltype = rownames(weight)) %>%
        cbind(as_tibble(weight)) %>%
        pivot_longer(!celltype, names_to = "bulk_id", values_to = "weight") %>%
        relocate(bulk_id)

    reshaped_bulk_counts <- tibble(bulk_id = colnames(bulk_counts)) %>%
        cbind(as_tibble(t(bulk_counts))) %>% as_tibble()

    cse <- compute_cse(reshaped_weights, ref_wide)
    proportions <- compute_proportion(reshaped_bulk_counts, cse)

    cse <- cse %>%
        pivot_wider(names_from = "celltype", values_from = "cse")
    return(list(proportions = proportions, cse = cse))
}