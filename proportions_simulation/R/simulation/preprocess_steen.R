# prepares steen data for simulation
preprocess_steen <- function(scfile, trainfrac) {
    # Use sc file to prepare for simulating bulks and consituting the reference matrix
    if (!file.exists(scfile)) {
        stop(paste0("Single cell file is not available in ",
            scfile,
            ". Please download from https://doi.org/10.5281/zenodo.10568550."))
    }
    sc <- readRDS(scfile)
    # use Steen et al. dataset
    sc$sc.pheno <- sc$sc.pheno %>% filter(origin == "Steen")
    sc$sc.counts <- sc$sc.counts[, sc$sc.pheno %>% pull(colnames)]
    scpheno <- as_tibble(sc$sc.pheno) %>%
        dplyr::select(celltype, colnames, Patient) %>%
        dplyr::rename(cell_id = colnames)
    sccounts <- tibble(cell_id = colnames(sc$sc.counts)) %>% cbind(as_tibble(t(sc$sc.counts)))

    # split by patients into a dataset for bulk creation and a single cell library:
    patients <- scpheno %>%
        pull(Patient) %>%
        unique()
    trainpatients <- patients[sample(trainfrac * length(patients), replace = FALSE)]
    npatients_sc <- length(trainpatients)
    scpheno <- scpheno %>%
        mutate(testtrain = ifelse(Patient %in% trainpatients, "train", "test")) %>%
        mutate(testtrain = as_factor(testtrain))
    sclong <- sccounts %>%
        pivot_longer(starts_with("ENSG"), names_to = "gene", values_to = "expression") %>%
        inner_join(scpheno)
    return(list(sclong = sclong, npatients = npatients_sc))
}

# filter only those cells contained in all patients
filter_patient_common_celltypes <- function(sclong) {
    patients <- sclong %>%
        pull(Patient) %>%
        unique()
    celltypes_isec <- c()
    count <- 0
    for (patient in patients) {
        celltypes_patient <- sclong %>%
            filter(Patient == patient) %>%
            pull(celltype) %>%
            unique()
        if(count > 0){
            celltypes_isec <- intersect(celltypes_patient, celltypes_isec)
        } else {
            celltypes_isec <- celltypes_patient
        }
        count <- count + 1
    }
    sclong <- sclong %>% filter(celltype %in% celltypes_isec)
    return(sclong)
}

