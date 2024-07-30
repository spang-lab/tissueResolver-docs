# Simulation benchmarking the performance of tissueResolver against MuSiC, Cibersort, BayesPrism, DTD and Bisque
# regarding correlation of cell type proportions across cell types

# Make sure to have run main_except_cibersort.R to have all other algorithms proportions stored in proportions_all
# Here we now add the proportions computed from CIBERSORT in simulation_run_cibersort_local.sh

# Here we generate the figures for cell type abundance in our publication using the MuSiC quality scores

# PARAMETERS

# paths
base_dir <- "data"
simulation <- "steen"
output_dir <- file.path(base_dir, simulation)

# repetition of runs
nruns <- 5

# base name of the .txt storing CIBERSORT fractions
cibersortfile <- "CIBERSORTx_Adjusted"


# DEPENDENCIES
library(tidyverse)
library(SingleCellExperiment)
library(MuSiC)
library(ggplot2)
library(gridExtra)

set.seed(42)

# load proportions computed for other algos
proportions_all <- readRDS(file.path(output_dir, "proportions_all.rds"))

# add CIBERSORT proportions (computed externally via DOCKER container)
source("R/simulation/benchmark_cibersort.R")
proportions_all <- add_cibersort_proportions(output_dir, proportions_all, cibersortfile, nruns)
saveRDS(proportions_all, file.path(output_dir, "proportions_all.rds"))

# now that proportions_all is complete: COMPUTE QUALITY SCORES
source("R/simulation/evaluate_proportions.R")
eval_statistics_all <- tibble()
eval_statistics_celltype_all <- tibble()
p_stat <- list()
p_ev <- list()
algos <- proportions_all %>% pull(algo) %>% unique()
if ("DTD" %in% algos) {
    bulk_indices <- proportions_all %>% filter(algo == "DTD") %>% pull(bulk_id) %>% unique()
} else {
    bulk_indices <- proportions_all %>%
        filter("DTD") %>%
        pull(bulk_id) %>%
        unique()
}
for (irun in 1:nruns) {
    test_bulks <- sample(bulk_indices, 20, replace = FALSE)
    message("Computing quality metrics for comparison")
    proportions_run <- proportions_all %>% filter(run == irun) %>% filter(bulk_id %in% test_bulks)
    proportion_matrices <- generate_proportion_matrices(proportions_run)
    # compute the scores
    ev <- perform_evaluation(proportion_matrices, output_dir, simulation)
    message(paste0("Computed quality metrics for run", irun))

    p_ev[[irun]] <- ev$plot + 
        ggtitle(paste0("Cell type fractions: Run ", irun)) + 
        ylab("Sample") + theme(axis.text.y = element_blank())
    eval_statistics <- ev$eval_statistics %>%
        as_tibble(rownames = "algo") %>%
        add_column(run = irun)
    eval_statistics_celltype <- ev$eval_statistics_celltype %>%
        add_column(run = irun)
    saveRDS(eval_statistics, file.path(output_dir, paste0("eval_statistics_run_", irun, ".rds")))
    saveRDS(ev$eval_statistics_celltype, file.path(output_dir, paste0("eval_statistics_celltype_run_", irun, ".rds")))
    p_stat[[irun]] <- tableGrob(eval_statistics)
    eval_statistics_all <- rbind(
        eval_statistics_all,
        eval_statistics)
    eval_statistics_celltype_all <- rbind(
        eval_statistics_celltype_all,
        eval_statistics_celltype
    )
}

ggsave(file.path(output_dir,
    paste0("heatmap_", simulation, "_all.pdf")),
    grid.arrange(grobs = p_ev, ncol = 1),
    width = 10, height = 20)
ggsave(file.path(output_dir,
    paste0("eval_statistics_", simulation, "_all.pdf")),
    grid.arrange(grobs = p_stat, ncol=1),
    width = 10, height = 15)

message("Saving all eval_statistics")
saveRDS(eval_statistics_all %>% arrange(algo), file.path(output_dir, "eval_statistics_all.rds"))

# SUMARRY STATISTICS
eval_statistics_all <- readRDS(file.path(output_dir, "eval_statistics_all.rds"))
eval_statistics_summary <- eval_statistics_all %>%
    group_by(algo) %>%
    mutate(mean_RMSD = mean(RMSD),
        mean_mAD = mean(mAD),
        mean_R = mean(R),
        sd_RMSD = sd(RMSD),
        sd_mAD = sd(mAD),
        sd_R = sd(R)) %>%
    ungroup() %>%
    dplyr::select(algo, starts_with("mean_"), starts_with("sd_")) %>%
    dplyr::distinct() %>%
    arrange(desc(mean_R))
saveRDS(eval_statistics_summary, file.path(output_dir, "eval_statistics_summary.rds"))
ggsave(
    file.path(
        output_dir,
        paste0("eval_statistics_summary.pdf")
    ),
    tableGrob(eval_statistics_summary),
    width = 10, height = 15
)

eval_summary_mean <- eval_statistics_summary %>%
    dplyr::select(algo, starts_with("mean_")) %>%
    pivot_longer(-algo, names_to = "metric", values_to = "value") %>%
    mutate(metric = str_extract(metric, "[^_]+$"))

eval_summary_sd <- eval_statistics_summary %>%
    dplyr::select(algo, starts_with("sd_")) %>%
    pivot_longer(-algo, names_to = "metric", values_to = "sd") %>%
    mutate(metric = str_extract(metric, "[^_]+$"))

eval_statistics_plot <- inner_join(eval_summary_mean, eval_summary_sd, by = c("algo", "metric"))


eval_statistics_celltype_summary <- eval_statistics_celltype_all %>%
    group_by(algo, celltype) %>%
    mutate(
        mean_R = mean(R),
        sd_R = sd(R)
    ) %>%
    ungroup() %>%
    dplyr::select(algo, celltype, starts_with("mean_"), starts_with("sd_")) %>%
    dplyr::distinct()
saveRDS(eval_statistics_celltype_summary, file.path(output_dir, "eval_statistics_summary.rds"))

# BAR PlOTS
metrics <- c("R", "RMSD", "mAD")
metric_plots <- list()
for (met in metrics) {
    bar_summary <- ggplot(eval_statistics_plot %>% filter(metric == met), aes(x = algo, y = value, fill = algo)) +
        geom_bar(
            stat = "identity", color = "black",
            position = position_dodge()
        ) +
        geom_errorbar(aes(ymin = value - sd, ymax = value + sd),
            width = .2,
            position = position_dodge(.9)
        ) + 
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        ggtitle(met)

    if (met == "R") {
        bar_summary <- bar_summary + coord_cartesian(ylim = c(0.6, 1))
    }
    metric_plots[[met]] <- bar_summary
}
ggsave(
    file.path(
        output_dir,
        paste0("bar_summary.pdf")
    ),
    grid.arrange(grobs = metric_plots, nrow = 1),
    width = 10, height = 5
)

bar_summary <- ggplot(eval_statistics_celltype_summary, aes(x = celltype, y = mean_R, fill = algo)) +
        geom_bar(
            stat = "identity", color = "black",
            position = position_dodge()
        ) +
        geom_errorbar(aes(ymin = mean_R - sd_R, ymax = mean_R + sd_R),
            width = .2,
            position = position_dodge(.9)
        ) +
        ggtitle("Celltype-wise correlation") +
        coord_cartesian(ylim = c(0.5, 1))
ggsave(
    file.path(
        output_dir,
        paste0("bar_summary_ct.pdf")
    ),
    bar_summary,
    width = 10, height = 5
)