source("R/simulation/helpers.R")

# SIMULATE GROUND TRUTH

# function to simulate nbulks bulks from sc data.
# Every bulk is simulated as a weighted sum over nperbulk single cell profiles.
# If fuzzy is not FALSE and fuzzyness != 0.0,
# the weights are drawn from a normal distribution around 1 with std "fuzzyness".
simulate_bulks <- function(scdata, nbulks, npatients, patients, fuzzy = NA, fuzzyness = 0.3) {
    if (fuzzyness < 0.0) {
        stop("negative value for fuzzyness.")
    }
    if (is.na(fuzzy)) {
        fuzzy <- (fuzzyness != 0.0)
    }
    # we uniformly sample profiles from the scdata frame
    # and cumulated weighted profiles to form each bulk
    if (!"cell_id" %in% names(scdata)) {
        stop("no cell_id column in scdata frame")
    }
    bulkcompo <- tibble()
    for (ibulk in 1:nbulks) {
        cell_ids <- sample_cells_patientwise(scdata, patients)
        if (fuzzy) {
            weights <- abs(rnorm(length(cell_ids), 1.0, fuzzyness))
        } else {
            weights <- rep(c(1.0), length(cell_ids))
        }
        bulkcompo <- bulkcompo %>% rbind(tibble(cell_id = cell_ids, bulk_id = ibulk, weight = weights))
    }
    bulks <- bulkcompo %>%
        inner_join(scdata %>% dplyr::select(cell_id, starts_with("ENSG")), by = "cell_id") %>%
        pivot_longer(starts_with("ENSG"), names_to = "gene", values_to = "expression") %>%
        mutate(expression = expression * weight) %>%
        dplyr::select(-weight) %>% # multiply by weight
        group_by(bulk_id, gene) %>%
        summarize(expression = sum(expression)) %>% # sum over each simulated bulk
        spread(gene, expression) %>% # bring to wide format
        ungroup()

    return(list(bulks = bulks, composition = bulkcompo))
}

# We sample a uniform percentage for each celltype in order to have
# each celltype covered in our simulated bulks
# However, the proportion we sample is distorted in order to arrive
# at qualitatively different bulks
sample_cells_patientwise <- function(scdata, patients) {

    cell_distribution_patients <- scdata %>% group_by(Patient) %>% dplyr::count(celltype)
    p <- sample(patients, 1, replace = FALSE)
    cell_distribution <- cell_distribution_patients %>% filter(Patient == p)
    cell_types <- cell_distribution %>% pull(celltype)
    all_cells <- scdata %>% dplyr::select(celltype, cell_id)
    cell_ids <- c()
    for (ct in cell_types) {
        cells_of_type <- all_cells %>%
            filter(celltype == ct) %>%
            pull(cell_id)
        n_cells <- round(
            # as we have npatients patients we take here the appropriate amount of
            # single cells distorted by a standard deviation of "half a patient's" cells
            pmin(
                abs(rnorm(
                    1, 1, 1
                )),
                1.5
            ) * cell_distribution[cell_distribution$celltype == ct, ]$n
        )
        cell_ids <- append(cell_ids, sample(cells_of_type, n_cells, replace = TRUE))
    }
    return(cell_ids)
}

# compute the true proportions given as the fraction
# of total cse (sum over all genes) compared to the
# total bulk expression
true_quantities <- function(bulks, sclong) {
    bulk_counts <- bulks$bulks
    composition <- bulks$composition
    cse <- compute_bulk_cse(composition, sclong)
    cse_output <- cse %>%
            pivot_wider(names_from = "celltype", values_from = "cse") %>%
            mutate(bulk_id = as.character(bulk_id))
    cse_total <- cse %>%
        group_by(bulk_id, celltype) %>%
        # compute total cse
        mutate(cse_sum = sum(cse)) %>%
        ungroup() %>%
        dplyr::select(bulk_id, celltype, cse_sum) %>%
        dplyr::distinct()

    bulk_total <- compute_bulk_total(bulk_counts)

    # compute fraction of cse over total bulk
    proportion <- cse_total %>%
        inner_join(bulk_total, by = "bulk_id") %>%
        mutate(prp = cse_sum / bulk_sum)

    # now we need to include also celltypes that
    # were not available in the composition libraray
    # these then obtain weight 0
    bulk_ids <- proportion %>%
        pull(bulk_id) %>%
        unique()
    all_celltypes <- sclong %>%
        dplyr::select(celltype) %>%
        unique()
    proportion_full <- tibble()
    for (b_id in bulk_ids) {
        prop <- proportion %>%
            filter(bulk_id == b_id) %>%
            full_join(all_celltypes, by = "celltype") %>%
            replace_na(list(bulk_id = b_id, prp = 0))
        proportion_full <- rbind(proportion_full, prop)
    }
    proportion_full <- proportion_full %>%
        dplyr::select(-cse_sum, -bulk_sum)

    return(list(proportions = proportion_full, cse = cse_output))
}