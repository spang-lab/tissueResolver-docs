# Simulation benchmarking the performance of tissueResolver against MuSiC, Cibersort, BayesPrism, DTD and Bisque
# regarding correlation of cell type proportions across cell types

# In this part we perform all algos except CIBERSORT, becuase it needs special
# tooling as it is not accessible as R script but only via a DOCKER container

# Make sure to run main_generate_input.R to have all data input files available

# PARAMETERS

# paths
base_dir <- "data"
simulation <- "steen"
output_dir <- file.path(base_dir, simulation)

# number of total simulated bulks
nbulks <- 100
# repetition of runs
nruns <- 5


# algorithms to benchmark:
# algos <- c("tr", "DTD", "music", "bp_subtypes", "Bisque")
algos <- c("DTD")

# DEPENDENCIES
library(tidyverse)
library(ggplot2)


if ("DTD" %in% algos) {
    library(DTD)
}
if ("music" %in% algos) {
    library(SingleCellExperiment)
    library(MuSiC)
}
if ("tr" %in% algos) {
    # Note that in tissueResolver we make use of
    # 50 cores for parallel computing specify the respective parts to your system as needed
    library(tissueResolver)
}
if ("bp_subtypes" %in% algos) {
    library(Seurat)
    # Note that in bayesPrism we make use of
    # 50 cores for parallel computing specify the respective parts to your system as needed
    library(BayesPrism)
}
if ("Bisque" %in% algos) {
    library(MIND)
}

# REPRODUCABILITY
set.seed(42)

# INFER PROPORTION FOR ALL ALGOS
proportions_all <- tibble()
# cse_all <- tibble()
for (irun in 1:nruns) {
    message(paste0("run ", irun, " / ", nruns))
    sc_library <- readRDS(file.path(output_dir, paste0("sc_library_run_", irun, ".rds")))
    sclong <- readRDS(file.path(output_dir, paste0("sclong_run_", irun, ".rds")))
    bulks_sim <- readRDS(file.path(output_dir, paste0("bulks_sim_run_", irun, ".rds")))

    proportions_true <- readRDS(file.path(output_dir, paste0("proportions_true_run_", irun, ".rds")))
    proportions_all <- rbind(
        proportions_all,
        proportions_true
    )

    if ("bp_subtypes" %in% algos) {
        source("R/simulation/benchmark_bayesprism.R")
        message("Computing subclustering of sc data for BayesPrism")
        sub_clustering <- compute_subclustering(sc_library)
        saveRDS(sub_clustering, file.path(output_dir, paste0("clustering_bp_run_", irun, ".rds")))
        message("Done computing subclustering")
    }

    # infer proportions of DTD
    if ("DTD" %in% algos) {
        source("R/simulation/benchmark_DTD.R")
        true_prp <- readRDS(file.path(output_dir, paste0("proportions_true_run_", irun, ".rds"))) %>% dplyr::select(-run, -algo)
        prep_DTD <- preprocess_DTD_input(bulks_sim, sclong, true_prp, sc_library)
        input_DTD <- generate_input_DTD(prep_DTD, 1:round(nbulks / 2), (round(nbulks / 2) + 1):nbulks)
        saveRDS(input_DTD, file.path(output_dir, paste0("input_DTD_run_", irun, ".rds")))
        output_DTD <- benchmark_dtd(
            input_DTD$reference,
            input_DTD$training_data_true,
            input_DTD$test_data_true
        )
        DTD_p <- output_DTD$proportions %>% add_column(run = irun, algo = "DTD")
        saveRDS(DTD_p, file.path(output_dir, paste0("proportions_DTD_run_", irun, ".rds")))
        proportions_all <- rbind(
            proportions_all,
            DTD_p
        )
    }

    # we test all 100 bulk samples (but later compare only among the train samples of DTD, i.e. 51: 2*nbulks)
    test_samples <- 1:nbulks

    # COMPETING ALGORITHMS on TEST SET
    if ("Bisque" %in% algos) {
        message("Inferring bisque proportions")
        source("R/simulation/benchmark_bisque.R")
        bisque_p <- benchmark_bisque(sc_library, bulks_sim$bulks, test_samples) %>%
            add_column(run = irun, algo = "Bisque")
        message("Done with Bisque proportions")
        saveRDS(bisque_p, file.path(output_dir, paste0("proportions_bisque_run_", irun, ".rds")))
        proportions_all <- rbind(
            proportions_all,
            bisque_p
        )
    }

    if ("music" %in% algos) {
        message("Inferring MuSiC proportions")
        source("R/simulation/benchmark_music.R")
        music_p <- benchmark_music(sc_library, sclong, bulks_sim$bulks, test_samples) %>%
            add_column(run = irun, algo = "music")
        message("Done with music proportions")
        saveRDS(music_p, file.path(output_dir, paste0("proportions_music_run_", irun, ".rds")))
        proportions_all <- rbind(
            proportions_all,
            music_p
        )
    }

    if ("bp_subtypes" %in% algos) {
        message("Inferring BayesPrism proportions")
        bp.res <- fit_bulks_bayesPrism(bulks_sim$bulks, sc_library, sub_clustering, test_samples)
        message("Done with BayesPrism")

        bp_proportions <- get_bayesPrism_proportions(bp.res) %>%
            add_column(run = irun, algo = "bp_subtypes")
        saveRDS(bp_proportions, file.path(output_dir, paste0("proportions_bayesPrism_run_", irun, ".rds")))
        proportions_all <- rbind(
            proportions_all,
            bp_proportions
        )
    }

    if ("tr" %in% algos) {
        source("R/simulation/benchmark_tr.R")
        message("Inferring proportions with tissueResolver")
        tr_output <- benchmark_tr(bulks_sim$bulks, sc_library, mapping, test_samples, irun)
        message("Done with tissueResolver")

        tr_proportions <- tr_output$cell_proportions %>% add_column(run = irun, algo = "tr")
        saveRDS(tr_proportions, file.path(output_dir, paste0("proportions_tissueResolver_run_", irun, ".rds")))
        proportions_all <- rbind(
            proportions_all,
            tr_proportions)
    }
}

message("Saving all proportions")
saveRDS(proportions_all, file.path(output_dir, "proportions_all.rds"))