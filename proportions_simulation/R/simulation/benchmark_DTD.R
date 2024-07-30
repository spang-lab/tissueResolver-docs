source("R/simulation/helpers.R")
source("R/simulation/true_quantities.R")

# the profiles in the reference for DTD are given as cell-type mean over sc_library
compute_reference_DTD <- function(sc_library) {
    reference <- tibble()
    clusters <- sc_library %>%
        pull(celltype) %>%
        unique()
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
    rownames(reference) <- reference_tibble %>% pull(celltype)
    reference <- t(reference)
    # reference <- sweep(reference, 2, 0.001 * colSums(reference), FUN = "/")
    return(reference)
}

# give the bulk expression in the appropriate format
compute_bulk_counts_DTD <- function(bulk_counts) {
    bulks_DTD <- bulk_counts %>%
        dplyr::select(-bulk_id) %>%
        as.matrix()
    rownames(bulks_DTD) <- bulk_counts %>% pull(bulk_id)
    bulks_DTD <- t(bulks_DTD)

    return(bulks_DTD)
}

# give the ground truth proportions in the adequate format
compute_bulk_pheno_DTD <- function(bulks, sclong, true_prp, sc_library) {
    # we only use the clusters available in the sc_library
    # as other clusters are not part of the deconvolution
    available_celltypes <- sc_library %>%
        pull(celltype) %>% unique()
    bulk_pheno_mat_true <- bulk_pheno_to_matrix(true_prp, available_celltypes)
    return(list(true = bulk_pheno_mat_true))

}

# this just converts the proportions tibble to a matrix
# serving as valuable input for DTD
bulk_pheno_to_matrix <- function(proportions, available_celltypes) {
    bulk_pheno <- proportions %>%
        filter(celltype %in% available_celltypes) %>%
        pivot_wider(names_from = bulk_id, values_from = prp) %>%
        replace(is.na(.), 0)
    bulk_pheno_mat <- bulk_pheno %>%
        dplyr::select(-celltype) %>%
        as.matrix()
    rownames(bulk_pheno_mat) <- bulk_pheno %>% pull(celltype)
    return(bulk_pheno_mat)
}

preprocess_DTD_input <- function(bulks, sclong, true_prp, sc_library) {
    # prepare DTD input
    message("Computing reference for DTD")
    reference <- compute_reference_DTD(sc_library)
    message("Computing bulk pheno for DTD")
    bulk_pheno <- compute_bulk_pheno_DTD(bulks, sclong, true_prp, sc_library)
    message("Computing bulk counts for DTD")
    bulk_counts <- compute_bulk_counts_DTD(bulks$bulks)
    input.DTD <- list(
        "reference" = reference,
        "bulk.counts" = bulk_counts,
        "bulk.pheno" = bulk_pheno
    )
    return(input.DTD)
}

generate_input_DTD <- function(prep_DTD, train_samples, test_samples) {
    # filter data by predefined train and test samples
    training_data_true <- list(
        "mixtures" = prep_DTD$bulk.counts[, train_samples],
        "quantities" = prep_DTD$bulk.pheno$true[, train_samples]
    )
    test_data_true <- list(
        "mixtures" = prep_DTD$bulk.counts[, test_samples],
        "quantities" = prep_DTD$bulk.pheno$true[, test_samples]
    )
    return(list(reference = prep_DTD$reference,
            train_samples = train_samples,
            test_samples = test_samples,
            training_data_true = training_data_true,
            test_data_true = test_data_true))
}

benchmark_dtd <- function(reference, training_data, test_data) {
    # traing DTD
    start.tweak <- rep(1, nrow(reference))
    names(start.tweak) <- rownames(reference)

    reference <- DTD::normalize_to_count(reference)
    training_data$mixtures <- DTD::normalize_to_count(training_data$mixtures)
    test_data$mixtures <- DTD::normalize_to_count(test_data$mixtures)

    # the correct ordering of genes is essential in order for the reference
    # to be compatible with the bulk mixtures
    gene_order <- rownames(reference)
    training_data <- list(
        "mixtures" = training_data$mixtures[gene_order, ],
        "quantities" = training_data$quantities
    )

    test_data <- list(
        "mixtures" = test_data$mixture[gene_order, ],
        "quantities" = test_data$quantities
    )

    message("DTD: Train DTD")
    DTD.DTD.x <- train_deconvolution_model(
        tweak = start.tweak,
        X.matrix = reference,
        train.data.list = training_data,
        test.data.list = NULL,
        estimate.c.type = "non_negative",
        use.implementation = "cxx",
        lambda.seq = "none",
        cv.verbose = TRUE,
        verbose = TRUE,
        NORM.FUN = "norm1"
    )


    # calculate c on test data
    message("DTD: Estimate c")
    esti.c.DTD.test <- estimate_c(
        X.matrix = reference,
        new.data = test_data$mixtures,
        DTD.model = DTD.DTD.x$best.model$Tweak,
        estimate.c.type = "non_negative"
    )
    esti.c.DTD.test[esti.c.DTD.test < 0] <- 0

    proportions_cse <- convert_weight_to_proportion(esti.c.DTD.test, reference, test_data$mixtures)
    proportions <- proportions_cse$proportions
    cse <- proportions_cse$cse

    return(list(proportions = proportions, cse = cse, model = DTD.DTD.x))
}