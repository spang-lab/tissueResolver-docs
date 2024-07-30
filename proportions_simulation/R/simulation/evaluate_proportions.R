# Here we evaluate the proportions output by the algorithms
# We use the metrics from the MuSic package (see their vignette)

# this function transforms the proportions dataframe output
# in our benchmarking to a matrix that can be input to the MuSic
# evaluation function
generate_proportion_matrices <- function(proportions) {
    algos <- proportions %>%
        pull(algo) %>%
        unique()
    proportions <- proportions %>% mutate(prp = ifelse(prp < 0, 0, prp))
    proportion_matrices <- list()
    for (alg in algos) {
        alg_proportions <- proportions %>%
            filter(algo == alg) %>%
            dplyr::select(bulk_id, celltype, prp) %>%
            pivot_wider(names_from = celltype, values_from = prp)
        bulk_ids <- alg_proportions %>% pull(bulk_id)
        alg_prp_mat <- alg_proportions %>%
            dplyr::select(-bulk_id) %>%
            as.matrix()
        rownames(alg_prp_mat) <- bulk_ids
        proportion_matrices <- append(proportion_matrices, list(alg_prp_mat))
    }

    names(proportion_matrices) <- algos
    return(proportion_matrices)

}

# This function computes quality scores for each algorithm and returns
# a proportions heatmap for each algorithm
perform_evaluation <- function(proportion_matrices, output_dir, simulation) {

    true_proportions <- proportion_matrices$true
    estimated_proportions <- proportion_matrices[names(proportion_matrices) != "true"]

    eval_statistics <- Eval_multi(
        prop.real = true_proportions,
        prop.est = estimated_proportions,
        method.name = names(estimated_proportions)
    )

    # compute also celltype-wise correlation
    eval_statistics_celltype <- tibble()
    celltypes <- colnames(true_proportions)
    for (alg in names(estimated_proportions)) {
        for(ct in celltypes) {
            co <- cor(true_proportions[, ct], estimated_proportions[[alg]][, ct])
            correlation_celltype <- tibble(algo = alg, celltype = ct, R = co)
            eval_statistics_celltype <- rbind(eval_statistics_celltype, correlation_celltype)
        }
    }
    message("Final evaluation statistics are: ")
    message(eval_statistics)

    algos <- names(proportion_matrices)

    proportion_matrices_sampled <- list()
    for (algo in algos){
        prop_algo <- proportion_matrices[[algo]]
        proportion_matrices_sampled <- append(proportion_matrices_sampled, list(prop_algo))
    }
    names(proportion_matrices_sampled) <- algos

    true_proportions <- proportion_matrices_sampled$true
    estimated_proportions <- proportion_matrices_sampled[names(proportion_matrices) != "true"]

    plot <- Prop_comp_multi(
        prop.real = true_proportions,
        prop.est = estimated_proportions,
        eval = FALSE,
        method.name = names(estimated_proportions)
    )
    return(list(eval_statistics = eval_statistics, eval_statistics_celltype = eval_statistics_celltype, plot = plot))
}