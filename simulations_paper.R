# Simulations benchmarking the performance of tissueResolver and BayesPrism

set.seed(42)

# Note that both in the bayesPrism and in the tissueResolver algorithm we make use of 50 cores for parallel computing
# specify the respective parts to your system as needed

# PARAMETERS
# a suffix appended to all output
suffix <- "cd8"
# algorithms to benchmark:
algos <- c("tr", "bp_nosubtypes", "bp_subtypes")
# single cell dataset (is expected to follow a certain annotation format, e.g. for cell labels)
scfile <- "data/sc.rds"
# celltype that is split into the "modified" and "unmodified" groups
celltype_to_split <- "T_CD8"
# fraction of cells (of that type) that should be in the modified class
changefrac <- 0.5
# number of genes that should be modified in the respective group (a list of percentages)
pct_genes_to_change <- c(2,5,10,17,20,30,40,50,60,70,80,90,95)
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
nruns <- 5

# ONLY FOR BAYESPRISM WITH SUBCLUSTERING 
cluster_resolution <- 0.5

library(tidyverse)
library(edgeR)
library(Seurat)
library(precrec)
library(ggplot2)
library(parallel)
if ("tr" %in% algos) {
  library(tissueResolver)
}

if ("bp_nosubtypes" %in% algos | "bp_subtypes" %in% algos) {
  library(BayesPrism)
}

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
	stop(paste0("Single cell file is not available in ", scfile, ". Please download from https://doi.org/10.5281/zenodo.10139154 or let the vignette run first."))
}
sc <- readRDS(scfile)
# use Steen et al. dataset
sc$sc.pheno <- sc$sc.pheno %>% filter(origin == "Steen")
sc$sc.counts <- sc$sc.counts[, sc$sc.pheno %>% pull(colnames)]
scpheno <- as_tibble(sc$sc.pheno) %>%
    dplyr::select(celltype, colnames, Patient) %>%
    rename(id = colnames)
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

n_cells_splittype <- scpheno %>%
    filter(celltype == celltype_to_split) %>%
    nrow()
n_cells <- scpheno %>% nrow()
sclong %>%
    filter(celltype == celltype_to_split) %>%
    group_by(gene) %>%
    mutate(zerofreq = sum(expression == 0) / n_cells_splittype)
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

fit_bulks_tr <- function(bulks, sccounts) {
    scmat <- sccounts %>%
        dplyr::select(starts_with("ENSG")) %>%
        as.matrix()
    rownames(scmat) <- sccounts %>% pull(id)
    bmat <- bulks$bulks %>%
        dplyr::select(-bulkid) %>%
        as.matrix()
    rownames(bmat) <- bulks$bulks %>% pull(bulkid)
    # call here the parallel version of fit tissue
    # for efficient computation of bulks in parallel on multiple cores
    tissuemodel <- tissueResolver::fit_tissue_par(
      t(bmat),
      t(scmat),
      maxit = 1e4,
      bootstrap = FALSE, # no bootstrapping
      ncores = 50
      )
    return(tissueResolver::specific_expression_regulation(tissuemodel, t(scmat), sccounts %>% dplyr::select(id, celltype) %>% rename(cell_id = id)))
}

compute_subclustering <- function(train_counts) {
	clustering <- train_counts %>% pull(id)
	names(clustering) <- clustering
	for (ctype in levels(factor(train_counts %>% pull(celltype)))) {
		exprmat <- train_counts %>%
            filter(celltype == ctype) %>%
            dplyr::select(starts_with("ENSG")) %>%
            as.matrix()
		rownames(exprmat) <- train_counts %>%
            filter(celltype == ctype) %>%
            pull(id)
		seurobj <- CreateSeuratObject(counts = t(exprmat), min.cells = 3, min.features = 5)
		if (dim(seurobj)[2] < 100) {
			cellids <- colnames(seurobj)
			clustering[cellids] <- ctype
			next
		}
		seurobj <- ScaleData(seurobj, features = rownames(seurobj))
		seurobj <- RunPCA(seurobj, features = rownames(seurobj))

		message("find neigbours")
		seurobj <- FindNeighbors(seurobj)
		message("find clusters")
		seurobj <- FindClusters(seurobj, resolution = cluster_resolution)

		clusters <- seurobj@meta.data %>%
            rownames_to_column() %>%
            as_tibble()
		cellids <- clusters %>% pull(rowname)
		clusternames <- paste0(ctype, "_", clusters %>% pull(seurat_clusters))
		clustering[cellids] <- clusternames
	}
	return(clustering)
}

fit_bulks_bayesPrism <- function(bulks, sccounts, clustering) {
	bk.dat <- bulks %>%
        dplyr::select(-bulkid, -mut) %>%
        as.matrix()
	rownames(bk.dat) <- bulks %>% pull(bulkid)
	sc.dat <- sccounts %>%
        dplyr::select(starts_with("ENSG")) %>%
        as.matrix()
	rownames(sc.dat) <- sccounts %>% pull(id)
		
	cell.type.labels <- sccounts %>% pull(celltype)
	names(cell.type.labels) <- sccounts %>% pull(id)

	if (is.null(clustering)) {
		cell.state.labels <- cell.type.labels
	} else {
		cell.state.labels <- clustering[sccounts %>% pull(id)]
	}

	myPrism <- new.prism(reference = sc.dat, mixture = bk.dat, input.type = "count.matrix", cell.type.labels = cell.type.labels, cell.state.labels = cell.state.labels, key = NULL, outlier.cut = 0.01, outlier.fraction = 0.1)
	bp.res <- run.prism(prism = myPrism, n.cores = 50)
	return(bp.res)
}

benchmark_tr <- function(train_counts, bulks_mut, bulks_nonmut, logfc, ngenes, irun, genelist) {
      cse_all <- tibble()
      # fit the bulks:
      cse_mut <- fit_bulks_tr(bulks_mut, train_counts)
      cse_nonmut <- fit_bulks_tr(bulks_nonmut, train_counts)

      # bind them together and put into edger:
      cse <- rbind(cse_mut %>% add_column(mut = TRUE), cse_nonmut %>% add_column(mut = FALSE))
      for (etype in c("expression", "regulation")) {
        # compute roc curves only for celltypes that are present in both mut and nonmut classes:
        for (ctype in intersect(cse_mut %>% pull(celltype) %>% unique(), cse_nonmut %>% pull(celltype) %>% unique())) {
          thiscse <- cse %>%
                filter(celltype == ctype) %>%
                dplyr::select(bulk_id, mut, gene, all_of(etype)) %>%
                pivot_wider(names_from = "gene", values_from = etype, values_fill = 0.0) %>%
                dplyr::select(mut, starts_with("ENSG"))
          sampleclasses <- thiscse %>%
                pull(mut) %>%
                as.factor()
	        csemat <- thiscse %>%
                dplyr::select(-mut) %>%
                as.matrix()
          if (nrow(csemat) > 2) { # dispersion cannot be fitted from only two bulks.
            tryCatch({
                edgerobj <- DGEList(counts = t(csemat), group = sampleclasses)
                mm <- model.matrix(~group, data = edgerobj$samples)
                # edgerobj <- calcNormFactors(edgerobj)
                edgerobj <- estimateDisp(edgerobj, mm, prior.df = 0)
                fit <- glmQLFit(edgerobj, mm)
                fit <- glmQLFTest(fit, coef = 2) # group 2 vs 1, mut vs nonmut
                computed_genelist <- fit$table %>%
                  mutate(gene = rownames(fit$table)) %>%
                  as_tibble()
                computed_genelist <- computed_genelist %>%
                  mutate(score = 1 - PValue) %>% # add score col
                  mutate(mutated = gene %in% genelist)
                precrec_obj <- evalmod(scores = computed_genelist$score, labels = computed_genelist$mutated)
                ggsave(paste0("data/plots/roc_cse_tr_", etype, "_sim", suffix, "_", ngenes, "genes_", logfc, "lfc_", ctype, "_run_", irun, ".pdf"), autoplot(precrec_obj))
                aucprc <- attr(precrec_obj, "aucs")
                cse_all <- rbind(cse_all, tibble(aucprc) %>% add_column(celltype = ctype, ereg = etype) %>% add_column(algo = "tr"))
              }, error = function(e) {
                  print(paste0("The fit of celltype", ctype, "in tissueResolver possesses at least one sample with zero library size and thus this celltype is excluded from downstream analysis."))
              })
          }
        }
      }
      return(cse_all)
}

benchmark_bayesPrism <- function(train_counts, bulks_mut, bulks_nonmut, clustering, logfc, ngenes, name, irun, genelist) {
      cse_all <- tibble()
      allbulks <- rbind(bulks_mut$bulks %>% add_column(mut = TRUE) %>% mutate(bulkid = paste0(bulkid, "_mut")), bulks_nonmut$bulks %>% add_column(mut = FALSE) %>% mutate(bulkid = paste0(bulkid, "_nonmut")))
      bp.res <- fit_bulks_bayesPrism(allbulks, train_counts, clustering)

      # compute roc curves only for celltypes that are present in both mut and nonmut classes:
      for (ctype in levels(factor(train_counts %>% pull(celltype)))) {
	      message(paste0("computing specific expression for ", ctype))
        z <- get.exp(bp = bp.res, state.or.type = "type", cell.name = ctype)
	      tryCatch({
          edgerobj <- DGEList(counts = t(z), group = allbulks$mut)
          mm <- model.matrix(~group, data = edgerobj$samples)
	      	edgerobj <- estimateDisp(edgerobj, mm, prior.df = 0)
	      	fit <- glmQLFit(edgerobj, mm)
	      	fit <- glmQLFTest(fit, coef = 2)
          computed_genelist <- fit$table %>%
            mutate(gene = rownames(fit$table)) %>%
            as_tibble()
          computed_genelist <- computed_genelist %>%
            mutate(score = 1 - PValue) %>% # add score col
            mutate(mutated = gene %in% genelist)
          precrec_obj <- evalmod(scores = computed_genelist$score, labels = computed_genelist$mutated)
          ggsave(paste0("data/plots/roc_cse_bp_", name, "_sim", suffix, "_", ngenes, "genes_", logfc, "lfc_", ctype, "_run_", irun, ".pdf"), autoplot(precrec_obj))
          aucprc <- attr(precrec_obj, "aucs")
          cse_all <- rbind(cse_all, tibble(aucprc) %>% add_column(celltype = ctype, ereg = "expression") %>% add_column(algo = paste0("bp_", name)))
	      }, error = function(e) {
            print(paste0("The fit of celltype ", ctype, " in bayesPrism possesses at least one sample with zero library size and thus this celltype is excluded from downstream analysis."))
        })
      }
      return(cse_all)
}


# ACTUAL SIMULATION
# varying number of genes and log fold changes
# in every iteration 
cse_all <- tibble()

for (pctchange in pct_genes_to_change) {
  ngenes <- round(pctchange * length(suitable_genes) / 100)
  genelist <- sample(suitable_genes, ngenes, replace = FALSE)
  for (logfc in log2foldchanges) {
            # this mutates only the selected genes in the selected cell population:
      mutsc <- sclong %>% rows_update(sclong %>% filter(gene %in% genelist) %>% # only selected genes
                                       filter(is_selected) %>% # only selected cells
                                       mutate(expression = (2^logfc) * expression), by = c("id", "gene"))
      train_counts <- mutsc %>%
            filter(testtrain == "train") %>%
            dplyr::select(id, celltype, expression, gene) %>%
            pivot_wider(names_from = "gene", values_from = "expression")
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

      if ("bp_subtypes" %in% algos) {
              message("computing clustering")
              clustering <- compute_subclustering(train_counts)
      }
      for (irun in 1:nruns) {
              print(paste0("run ", irun, " / ", nruns, " logfc = ", logfc, ", ngenes = ", ngenes, " (", pctchange, " % of all)"))

        # now, simulate bulks:
        bulks_mut <- simulate_bulks(test_and_mod_counts %>% dplyr::select(id, starts_with("ENSG")), nbulks, nperbulk, fuzzyness = 0.3)
        bulks_nonmut <- simulate_bulks(test_and_nonmod_counts %>% dplyr::select(id, starts_with("ENSG")), nbulks, nperbulk, fuzzyness = 0.3)

        if ("bp_nosubtypes" %in% algos) {
        cse_all <- rbind(cse_all, benchmark_bayesPrism(train_counts, bulks_mut, bulks_nonmut, NULL, logfc, ngenes, "nosub", irun, genelist) %>% add_column(run = irun, logfc = logfc, ngenes = ngenes, pct_genes = ngenes / length(suitable_genes) * 100))
        }
            if ("bp_subtypes" %in% algos) {
              message("fitting")
        cse_all <- rbind(cse_all, benchmark_bayesPrism(train_counts, bulks_mut, bulks_nonmut, clustering, logfc, ngenes, "sub", irun, genelist) %>% add_column(run = irun, logfc = logfc, ngenes = ngenes, pct_genes = ngenes / length(suitable_genes) * 100))
        }
        if ("tr" %in% algos) {
          cse_all <- rbind(
                      cse_all,
                            benchmark_tr(train_counts, bulks_mut, bulks_nonmut, logfc, ngenes, irun, genelist) %>% 
                            add_column(run = irun, logfc = logfc, ngenes = ngenes, pct_genes = ngenes / length(suitable_genes) * 100)
                  )
        }
      }
      saveRDS(cse_all, file = paste0("data/plots/sim_roc_", suffix, "_up_to_", logfc, "_", pctchange, ".rds"))
  }
}
saveRDS(cse_all, file = paste0("data/plots/sim_roc", suffix, ".rds"))