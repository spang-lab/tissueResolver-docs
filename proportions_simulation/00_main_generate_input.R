# Simulation benchmarking the performance of tissueResolver against MuSiC, Cibersort, BayesPrism, DTD and Bisque
# regarding correlation of cell type proportions across cell types

# We use real sc data (Steen DLBCL https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182436)
# to simulate bulks by summing over sc profiles.
# The sc data also serves as reference in the deconvolution tools.
# This simulation then computes the inferred proportions and correlates them to the ground truth
# proportions.

# DEPENDENCIES
library(tidyverse)

# paths
base_dir <- "data"
simulation <- "steen"
output_dir <- file.path(base_dir, simulation)

# input files
scfile <- file.path(base_dir, "sc.rds")
nbulks <- 100

# fraction of non-modified to modified cells of the selected cell type
trainfrac <- 0.5

# repetition of runs
nruns <- 5

# this is the fraction of patients we use in order to generate
# artifical bulks from their single cells used for TRAINIG DTD
# in order to have a distinct train/test split on the patient level
bulk_train_patient_fraction <- 0.6
# this is the relative amount wrt nbulks of artificial bulks
# we use to TRAIN DTD on
bulk_train_fraction <- 0.8

# REPRODUCABILITY
set.seed(42)

# SPLIT FOR BULK AND SC
# the number of patients in the sc library is important for the bulk construction
# in simulate_bulks

# GENERATE INPUT FOR ALGORITHMS
proportions_all <- tibble()

for (irun in 1:nruns) {
    message(paste0("run ", irun, " / ", nruns))
    source("R/simulation/preprocess_steen.R")
    message("Preparing data")
    prep <- preprocess_steen(scfile, trainfrac)
    sclong <- prep$sclong
    npatients <- prep$npatients

    # the source single cell df for our simulation
    sclong <- filter_patient_common_celltypes(sclong)
    saveRDS(sclong, file.path(output_dir, paste0("sclong_run_", irun, ".rds")))

    # the part of the sc that simulates the bulks
    sc_bulksim <- sclong %>%
        filter(testtrain == "test") %>%
        dplyr::select(cell_id, celltype, expression, gene, Patient) %>%
        pivot_wider(names_from = "gene", values_from = "expression")
    saveRDS(sc_bulksim, file.path(output_dir, paste0("sc_bulks_run_", irun, ".rds")))

    # the part of the sc that constitutes the sc library
    sc_library <- sclong %>%
        filter(testtrain == "train") %>%
        dplyr::select(cell_id, celltype, expression, gene) %>%
        pivot_wider(names_from = "gene", values_from = "expression")
    saveRDS(sc_library, file.path(output_dir, paste0("sc_library_run_", irun, ".rds")))

    # simulate bulks
    source("R/simulation/true_quantities.R")
    message("Simulating bulks")

    # in order to separate the training samples for DTD on the patient level
    bulk_patients <- sc_bulksim %>% pull(Patient) %>% unique()
    train_pos <- sample(1:length(bulk_patients),
            round(bulk_train_patient_fraction * length(bulk_patients)),
            replace = FALSE)
    train_patients <- bulk_patients[train_pos]
    test_patients <- bulk_patients[-train_pos]
    amount_train_bulks <- round(bulk_train_fraction * nbulks)
    amount_test_bulks <- nbulks - amount_train_bulks
    # generate now train/test artificial bulks
    train_bulks_sim <- simulate_bulks(sc_bulksim,
                amount_train_bulks,
                npatients,
                train_patients,
                fuzzyness = 0.3)
    test_bulks_sim <- simulate_bulks(sc_bulksim,
                amount_test_bulks,
                npatients,
                test_patients,
                fuzzyness = 0.3)
    # Combine train and test bulks to a single tibble with a training and test block
    test_bulks_sim$bulks <- test_bulks_sim$bulks %>% mutate(bulk_id = bulk_id + amount_train_bulks)
    test_bulks_sim$composition <- test_bulks_sim$composition %>% mutate(bulk_id = bulk_id + amount_train_bulks)
    split <- c()
    split$train_samples <- 1:amount_train_bulks
    split$test_samples <- (amount_train_bulks + 1):nbulks

    bulks_sim <- c()
    bulks_sim$bulks <- rbind(train_bulks_sim$bulks, test_bulks_sim$bulks)
    bulks_sim$composition <- rbind(train_bulks_sim$composition, test_bulks_sim$composition)
    message("Done simulating bulks")
    saveRDS(bulks_sim, file.path(output_dir, paste0("bulks_sim_run_", irun, ".rds")))

    # preprocess input for cibersort
    source("R/simulation/benchmark_cibersort.R")
    test_samples <- bulks_sim$bulks %>%
        pull(bulk_id) %>%
        unique()
    input_cibersort <- compute_input_cibersort(sc_library, bulks_sim$bulks, test_samples)
    write_cibersort_input_txt(
        input_cibersort$reference,
        input_cibersort$refsample,
        input_cibersort$bulk_counts$bulk_counts_matrix,
        output_dir,
        irun
    )

    # infer the true proportions
    message("Inferring true proportions")
    true <- true_quantities(bulks_sim, sclong)
    true_prp <- true$proportions
    message("Done with true proportions")
    tp_save <- true_prp %>% add_column(run = irun, algo = "true")
    saveRDS(tp_save, file.path(output_dir, paste0("proportions_true_run_", irun, ".rds")))
    saveRDS(true$cse, file.path(output_dir, paste0("cse_true_run_", irun, ".rds")))
    proportions_all <- rbind(
            proportions_all,
            tp_save)
}

# this tibble stores the true proportion
message("Saving all true proportions")
saveRDS(proportions_all, file.path(output_dir, "proportions_all.rds"))