# Preparing input for csDGE simulations
set.seed(42)

# PARAMETERS
# a suffix appended to all output
suffix <- "B"
# suffix <- "cd8"
# single cell dataset (is expected to follow a certain annotation format, e.g. for cell labels)
scfile <- "data/sc.rds"
# celltype that is split into the "modified" and "unmodified" groups
celltype_to_split <- "B"
# celltype_to_split <- "T_CD8"
# fraction of cells (of that type) that should be in the modified class
changefrac <- 0.5
# number of genes that should be modified in the respective group (a list of percentages)
pct_genes_to_change <- c(2, 5, 10, 17, 20, 30, 40, 50, 60, 70, 80, 90, 95)
# log2 fold changes: selected genes in the selected cell type are multiplied by this factor
log2foldchanges <- c(0.5, 0.8, 1.0, 1.5, 2)
# gene selection:
# minimal expression value in the single cell dataset for a gene to be considered "expressed"
mod_genes_min_expression <- 1.0
# every gene with less than this fraction of zeros in the selected cell type will be suitable for modification,
# i.e. only genes expressed in a larger cell population will be changed to avoid singling out cells and that the genes act as marker genes.
# e.g.: if the gene is not expressed in less than 40 % of all cells in all T cells, then it is considered for modification
mod_genes_max_zero_freq <- 0.4
# if set to true, demand that the above condition are also met in all other cell types
# this will mean that the genes are not specific for the selected cell type, again avoiding that marker genes are picked.
apply_gene_filtering_to_all_cells <- TRUE
# fraction of non-modified to modified cells of the selected cell type
trainfrac <- 0.5

# statistics parameters
# number of bulks in each of the two (modified and unmodified) classes
nbulks <- 50
# number of (randomly selected) single cells to sum over to form one bulk profile
nperbulk <- 500
# repetition of runs, to compute errors and means.
nruns <- 1

# ONLY FOR BAYESPRISM WITH SUBCLUSTERING
cluster_resolution <- 0.5

library(tidyverse)

# function to simulate nbulks bulks from single cell data.
# every bulk is simulated as a weighted sum over nperbulk single cell profiles.
# if fuzzy is not FALSE and fuzzyness != 0.0, the weights are drawn from a normal distribution around 1 with width "fuzzyness".
simulate_bulks <- function(scdata, nbulks, nperbulk, fuzzy = NA, fuzzyness = 0.3) {
  if (fuzzyness < 0.0) {
    stop("negative value for fuzzyness.")
  }
  if (is.na(fuzzy)) {
    fuzzy <- (fuzzyness != 0.0)
  }
  # we uniformly sample from the scdata frame
  if (!"id" %in% names(scdata)) {
    stop("no id column in scdata frame")
  }
  bulkcompo <- tibble()
  # all other columns are assumed to be genes.
  all_ids <- scdata %>% pull(id)
  for (ibulk in 1:nbulks) {
    ids <- sample(all_ids, nperbulk, replace = TRUE)
    if (fuzzy) {
      weights <- abs(rnorm(nperbulk, 1.0, fuzzyness))
    } else {
      weights <- rep(c(1.0), nperbulk)
    }
    bulkcompo <- bulkcompo %>% rbind(tibble(id = ids, bulkid = ibulk, weight = weights))
  }
  bulks <- bulkcompo %>%
        inner_join(scdata) %>%
        gather(gene, expression, -id, -bulkid, -weight) %>%
    mutate(expression = expression * weight) %>%
        dplyr::select(-weight) %>% # multiply by weight
    group_by(bulkid, gene) %>%
        summarize(expression = sum(expression) / nperbulk) %>% # sum over each simulated bulk
    spread(gene, expression) %>% # bring to wide format
    ungroup()
  return(list(bulks = bulks, composition = bulkcompo))
}


# READ SC FILE
if (!file.exists(scfile)) {
 stop(paste0("Single cell file is not available in ", scfile, ". Please download from https://doi.org/10.5281/zenodo.10568550 or let the vignette run first."))
}
sc <- readRDS(scfile)
# use Steen et al. dataset
sc$sc.pheno <- sc$sc.pheno %>% filter(origin == "Steen")
sc$sc.counts <- sc$sc.counts[, sc$sc.pheno %>% pull(colnames)]
scpheno <- as_tibble(sc$sc.pheno) %>%
    dplyr::select(celltype, colnames, Patient) %>%
    dplyr::rename(id = colnames)
sccounts <- tibble(id = colnames(sc$sc.counts)) %>% cbind(as_tibble(t(sc$sc.counts)))

# split by patients into a dataset for bulk creation and a single cell library:
patients <- scpheno %>%
    pull(Patient) %>%
    unique()
trainpatients <- patients[sample(trainfrac * length(patients), replace = FALSE)]
scpheno <- scpheno %>%
    mutate(testtrain = ifelse(Patient %in% trainpatients, "train", "test")) %>%
    mutate(testtrain = as_factor(testtrain))

# now select cells to change:
all_these_cells <- scpheno %>%
    filter(celltype == celltype_to_split) %>%
    pull(id)
selected_cells <- sample(all_these_cells, changefrac * length(all_these_cells), replace = FALSE)
# rename "mutated" or changed celltypes with suffix "m"
scpheno <- scpheno %>% mutate(is_selected = id %in% selected_cells)

# find genes suitable to change:
sclong <- sccounts %>%
    pivot_longer(starts_with("ENSG"), names_to = "gene", values_to = "expression") %>%
    inner_join(scpheno)
saveRDS(sclong, paste0("data/sclong_", suffix,".rds"))

n_cells_splittype <- scpheno %>%
    filter(celltype == celltype_to_split) %>%
    nrow()
n_cells <- scpheno %>% nrow()
suitable_genes <- sclong %>%
  filter(celltype == celltype_to_split) %>%
    group_by(gene) %>%
    mutate(zerofreq = sum(expression < mod_genes_min_expression) / n_cells_splittype) %>%
    filter(zerofreq < mod_genes_max_zero_freq) %>%
    pull(gene) %>%
    unique()
if (apply_gene_filtering_to_all_cells) {
    suitable_genes_in_all_cells <- sclong %>%
        group_by(gene) %>%
        mutate(zerofreq = sum(expression < mod_genes_min_expression) / n_cells) %>%
        filter(zerofreq < mod_genes_max_zero_freq) %>%
        pull(gene) %>%
        unique()

  suitable_genes <- intersect(suitable_genes, suitable_genes_in_all_cells)
  print("selected genes restricting on ALL cell types")
} else {
  print(paste("selected genes restricting on", celltype_to_split))
}

print(paste("selected ngenes: ", length(suitable_genes)))


# HELPERS TO GENERATE THE CIBERSORT INPUT TXT FILES
# (needed as input for DOCKER container)
compute_refsample_cibersort <- function(train_counts) {
  celltypes <- train_counts %>% pull(celltype)
  refsample <- train_counts %>%
    dplyr::select(starts_with("ENSG")) %>%
    as.matrix()
  refsample <- t(refsample)
  gene_symbols <- rownames(refsample)
  refsample_gene_names <- cbind("gene_name" = gene_symbols, refsample)
  refsample_matrix <- rbind(c("name", celltypes), refsample_gene_names)
  return(refsample_matrix)
}

compute_bulk_counts_cibersort <- function(bulk_counts) {
  bulk_matrix <- bulk_counts %>%
    dplyr::select(-bulkid) %>%
    as.matrix()
  sample_names <- bulk_counts %>% pull(bulkid)
  bulk_matrix <- t(bulk_matrix)
  gene_symbols <- rownames(bulk_matrix)
  bulk_counts_gene_names <- cbind("gene_name" = gene_symbols, bulk_matrix)
  bulk_counts_matrix <- rbind(c("name", sample_names), bulk_counts_gene_names)
  return(bulk_counts_matrix)
}


compute_input_cibersort <- function(train_counts, bulk_counts) {
  refsample_matrix <- compute_refsample_cibersort(train_counts)
  bulk_counts <- compute_bulk_counts_cibersort(bulk_counts)
  return(list(
    refsample = refsample_matrix,
    bulk_counts = bulk_counts
  ))
}

write_cibersort_input_txt <- function(refsample, bulk_counts, logfc, pctchange, irun) {
  write.table(refsample,
    file = file.path("data", paste0("refsample_fc_", logfc, "_pct_", pctchange,"_run_", irun,"_", suffix, ".txt")),
    quote = FALSE,
    sep = "\t ",
    col.names = FALSE,
    row.names = FALSE
  )
  write.table(bulk_counts,
    file = file.path("data", paste0("bulks_fc_", logfc, "_pct_", pctchange,"_run_", irun,"_", suffix, ".txt")),
    quote = FALSE,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE
  )
}

prepare_input_cibersort <- function(train_counts, bulks_mut, bulks_nonmut, logfc, pctchange, irun) {
  allbulks <- rbind(
    bulks_mut$bulks %>% mutate(bulkid = paste0(bulkid, "_mut")),
    bulks_nonmut$bulks %>% mutate(bulkid = paste0(bulkid, "_nonmut"))
  )
  input_cibersort <- compute_input_cibersort(train_counts, allbulks)
  write_cibersort_input_txt(
    input_cibersort$refsample,
    input_cibersort$bulk_counts,
    logfc,
    pctchange,
    irun
  )
}


# ACTUAL DATA GENERATION
for (pctchange in pct_genes_to_change) {
  ngenes <- round(pctchange * length(suitable_genes) / 100)
  # gene list storing actually modified genes
  genelist <- sample(suitable_genes, ngenes, replace = FALSE)
  saveRDS(genelist, paste0("data/genelist_", pctchange, "_", suffix, ".rds"))
  for (logfc in log2foldchanges) {
      # this mutates only the selected genes in the selected cell population:
      mutsc <- sclong %>% rows_update(sclong %>% filter(gene %in% genelist) %>% # only selected genes
                                       filter(is_selected) %>% # only selected cells
                                       mutate(expression = (2^logfc) * expression), by = c("id", "gene"))
      train_counts <- mutsc %>%
            filter(testtrain == "train") %>%
            dplyr::select(id, celltype, expression, gene) %>%
            pivot_wider(names_from = "gene", values_from = "expression")
      saveRDS(train_counts, paste0("data/train_counts_", pctchange, "_", logfc, "_", suffix, ".rds"))

      # from the celltype_to_split take only the modified ones
      test_and_mod_counts <- mutsc %>%
            filter(testtrain == "test") %>%
            filter((celltype != celltype_to_split) | is_selected) %>%
            dplyr::select(id, celltype, expression, gene) %>%
            pivot_wider(names_from = "gene", values_from = "expression")
      # from the celltype_to_split take only the non-modified ones
      test_and_nonmod_counts <- mutsc %>%
            filter(testtrain == "test") %>%
            filter((celltype != celltype_to_split) | !is_selected) %>%
            dplyr::select(id, celltype, expression, gene) %>%
            pivot_wider(names_from = "gene", values_from = "expression")

    for (irun in 1:nruns) {
      print(paste0("run ", irun, " / ", nruns, " logfc = ", logfc, ", ngenes = ", ngenes, " (", pctchange, " % of all)"))

      # now, simulate bulks:
      bulks_mut <- simulate_bulks(test_and_mod_counts %>% dplyr::select(id, starts_with("ENSG")), nbulks, nperbulk, fuzzyness = 0.3)
      saveRDS(bulks_mut, paste0("data/bulks_mut_", pctchange, "_", logfc, "_run_", irun, "_", suffix, ".rds"))
      bulks_nonmut <- simulate_bulks(test_and_nonmod_counts %>% dplyr::select(id, starts_with("ENSG")), nbulks, nperbulk, fuzzyness = 0.3)
      saveRDS(bulks_nonmut, paste0("data/bulks_nonmut_", pctchange, "_", logfc, "_run_", irun, "_", suffix, ".rds"))

      if(pctchange == 17 || logfc == 0.8) {
        prepare_input_cibersort(train_counts, bulks_mut, bulks_nonmut, logfc, pctchange, irun)
      }
    }
  }
}
