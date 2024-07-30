# Simulations benchmarking the performance of tissueResolver and BayesPrism
# Make sure to run simulations_paper_generate.R in order to generate all input data needed here

# Note that for tissueResolver, bayesPrism, bMind and ISLET algorithm we make use of 50 cores for parallel computing
# specify the respective parts to your system as needed

# PARAMETERS
# a suffix appended to all output
suffix <- "B"
# suffix <- "cd8"

# algorithms to benchmark:
# Note that if cibersort_hires is contained here, the output .txt files of the cibersort/hires DOCKER container
# need to be pre computed via simulation_run_cibersort_tr.sh
algos <- c("tr", "bp_nosubtypes", "bp_subtypes", "bmind", "islet", "islet_cibersort", "cibersort_hires")
# single cell dataset (is expected to follow a certain annotation format, e.g. for cell labels)
scfile <- "data/sc.rds"
# celltype that is split into the "modified" and "unmodified" groups
celltype_to_split <- "B"
# celltype_to_split <- "T_CD8"
# number of genes that should be modified in the respective group (a list of percentages)
pct_genes_to_change <- c(2, 5, 10, 17, 20, 30, 40, 50, 60, 70, 80, 90, 95)
# log2 fold changes: selected genes in the selected cell type are multiplied by this factor
log2foldchanges <- c(0.5, 0.8, 1.0, 1.5, 2)

# statistics parameters
# number of bulks in each of the two (modified and unmodified) classes, i.e., total amount of bulks is 2*nbulks
nbulks <- 50
# repetition of runs, to optionally compute errors and means.
nruns <- 1

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

if ("bmind" %in% algos) {
  library(MIND)
}

if ("islet" %in% algos || "islet_cibersort" %in% algos) {
  library(ISLET)
  library(SummarizedExperiment)
}

# TISSUERESOLVER

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
    tissuemodel <- tissueResolver::fit_tissue(
      t(bmat),
      t(scmat),
      maxit = 1e4,
      bootstrap = FALSE, # no bootstrapping
      ncores = 50)
    return(tissueResolver::specific_expression_regulation(tissuemodel, t(scmat), sccounts %>% dplyr::select(id, celltype) %>% dplyr::rename(cell_id = id)))
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
        tryCatch(
          {
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
          },
          error = function(e) {
            print(paste0("The fit of celltype", ctype, "in tissueResolver possesses at least one sample with zero library size and thus this celltype is excluded from downstream analysis."))
          }
        )
      }
    }
  }
  return(cse_all)
}

# BAYESPRISM

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
    # Seurat v5 needs this preprocessing for PCA, due to somewhat new required layers
    # in SeuratObjects, see: https://satijalab.org/seurat/articles/seurat5_integration
    seurobj <- NormalizeData(seurobj)
    seurobj <- FindVariableFeatures(seurobj)
    seurobj <- ScaleData(seurobj)
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

benchmark_bayesPrism <- function(train_counts, bulks_mut, bulks_nonmut, clustering, logfc, ngenes, name, irun, genelist) {
  cse_all <- tibble()
  allbulks <- rbind(bulks_mut$bulks %>% add_column(mut = TRUE) %>% mutate(bulkid = paste0(bulkid, "_mut")), bulks_nonmut$bulks %>% add_column(mut = FALSE) %>% mutate(bulkid = paste0(bulkid, "_nonmut")))
  bp.res <- fit_bulks_bayesPrism(allbulks, train_counts, clustering)

  # compute roc curves only for celltypes that are present in both mut and nonmut classes:
  for (ctype in levels(factor(train_counts %>% pull(celltype)))) {
    message(paste0("computing specific expression for ", ctype))
    z <- get.exp(bp = bp.res, state.or.type = "type", cell.name = ctype)
    tryCatch(
      {
        edgerobj <- DGEList(counts = t(z), group = allbulks$mut)
        # edgerobj <- normLibSizes(edgerobj)
        mm <- model.matrix(~group, data = edgerobj$samples)
        # edgerobj <- calcNormFactors(edgerobj)
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
      },
      error = function(e) {
        print(paste0("The fit of celltype ", ctype, " in bayesPrism possesses at least one sample with zero library size and thus this celltype is excluded from downstream analysis."))
      }
    )
  }
  return(cse_all)
}

# BMIND

prepare_bmind_input <- function(bulk_counts, sclong, sc_library) {
  bulks <- bulk_counts %>%
    dplyr::select(-bulkid) %>%
    # bmind can not deal with genes identical to 0
    dplyr::select(where(~ n_distinct(.) > 1)) %>%
    as.matrix()
  rownames(bulks) <- bulk_counts %>% pull(bulkid)
  bulks <- t(bulks)

  sc_library <- sc_library %>% inner_join(sclong %>% dplyr::select(id, Patient) %>% dplyr::distinct(), by = "id") %>% filter(id %in% modified_ids)

  available_celltypes <- sc_library %>%
    pull(celltype) %>%
    unique()

  sc <- sc_library %>%
    dplyr::select(starts_with("ENSG")) %>%
    as.matrix()
  rownames(sc) <- sc_library %>% pull(id)
  sc <- t(sc)

  sc_meta <- sc_library %>%
    dplyr::select(Patient, celltype) %>%
    dplyr::rename(sample = Patient, cell_type = celltype) %>%
    data.frame()
  colnames(sc_meta) <- c('sample', 'cell_type')

    sc_meta_bisque <- sc_library %>%
      dplyr::select(id, celltype) %>%
      # note: subject in Bisque denotes the cell id NOT the patient
      dplyr::rename(subject = id, cell_type = celltype) %>%
      data.frame()
    row.names(sc_meta_bisque) <- sc_meta_bisque$subject

  return(list(bulk_counts = bulks, sc = sc, sc_meta = sc_meta, sc_meta_bisque = sc_meta_bisque))
}

fit_bulks_bmind <- function(train_counts, sclong, bulks_mut, bulks_nonmut, logfc, ngenes, irun, genelist) {
    # data preparation
    bulks_case_ctrl <- rbind(
      bulks_mut$bulks,
      bulks_nonmut$bulks
    )
    bulks_case_ctrl <- bulks_case_ctrl %>% mutate(bulkid = 1:(2 * nbulks))
    # case/ctrl pheno for bmind
    y <- c(rep(1, 50), rep(0, 50))
    input_bmind <- prepare_bmind_input(bulks_case_ctrl, sclong, train_counts)

    # on our data this produces only negative definite covariance matrices (for all genes),
    # thus it can't be used
    # prior <- get_prior(sc = input_bmind$sc, meta_sc = input_bmind$sc_meta)


    # infer proportions with Bisque from sc data
    bisque_proportions <- est_frac_sc(
      bulk = input_bmind$bulk_counts,
      frac_method = "Bisque",
      sc_count =input_bmind$sc,
      sc_meta = input_bmind$sc_meta_bisque
    )

    decon <- bMIND(bulk = input_bmind$bulk_counts,
      frac = bisque_proportions,
      y = y,
      ncore = 50)
    
    # direct compute pvalues with built in method
    bmind_pval <- t(decon$pval)
  return(bmind_pval)
}

benchmark_bmind <- function(train_counts, sclong, bulks_mut, bulks_nonmut, logfc, ngenes, irun, genelist) {
  cse_all <- tibble()
  bmind_pval <- fit_bulks_bmind(train_counts, sclong, bulks_mut, bulks_nonmut, logfc, ngenes, irun, genelist)
  for (ctype in colnames(bmind_pval)) {
    computed_genelist <- bmind_pval %>%
      as_tibble(rownames = "gene") %>%
      dplyr::select(gene, all_of(ctype)) %>%
      mutate(score = 1 - !!as.symbol(ctype)) %>%
      mutate(mutated = gene %in% genelist)
    precrec_obj <- evalmod(scores = computed_genelist$score, labels = computed_genelist$mutated)
    ggsave(paste0("data/plots/roc_cse_bmind_sim_", suffix, "_", ngenes, "genes_", logfc, "lfc_", ctype, "_run_", irun, ".pdf"), autoplot(precrec_obj))
    aucprc <- attr(precrec_obj, "aucs")
    cse_all <- rbind(cse_all, tibble(aucprc) %>% add_column(celltype = ctype, ereg = "expression") %>% add_column(algo = "bmind"))
  }
  return(cse_all)
}

# ISLET

# the following functions are helpers to determine ground truth proportions
# that can be used in the benchamrking of ISLET if use_true_proportions=TRUE 

# the total bulk count for a single bulk
# is the sum over all gene counts within that bulk
compute_bulk_total <- function(bulk_counts) {
  bulk_total <- bulk_counts %>%
    pivot_longer(cols = starts_with("ENSG"), names_to = "gene", values_to = "expression") %>%
    group_by(bulkid) %>%
    mutate(bulk_sum = sum(expression)) %>%
    ungroup() %>%
    dplyr::select(bulkid, bulk_sum) %>%
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
    inner_join(scwide, by = "id") %>%
    pivot_longer(cols = starts_with("ENSG"), names_to = "gene", values_to = "expression") %>%
    # the cse equals each weighted sc profile cumulated over the respective cell type
    mutate(expression = expression * weight) %>%
    group_by(bulkid, celltype, gene) %>%
    mutate(cse = sum(expression)) %>%
    ungroup()
  cse <- celltype_comp %>%
    dplyr::select(bulkid, gene, cse, celltype) %>%
    dplyr::distinct()
  return(cse)
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
    mutate(bulkid = as.character(bulkid))
  cse_total <- cse %>%
    group_by(bulkid, celltype) %>%
    # compute total cse
    mutate(cse_sum = sum(cse)) %>%
    ungroup() %>%
    dplyr::select(bulkid, celltype, cse_sum) %>%
    dplyr::distinct()

  bulk_total <- compute_bulk_total(bulk_counts)

  # compute fraction of cse over total bulk
  proportion <- cse_total %>%
    inner_join(bulk_total, by = "bulkid") %>%
    mutate(prp = cse_sum / bulk_sum)

  # now we need to include also celltypes that
  # were not available in the composition libraray
  # these then obtain weight 0
  bulk_ids <- proportion %>%
    pull(bulkid) %>%
    unique()
  all_celltypes <- sclong %>%
    dplyr::select(celltype) %>%
    unique()
  proportion_full <- tibble()
  for (b_id in bulk_ids) {
    prop <- proportion %>%
      filter(bulkid == b_id) %>%
      full_join(all_celltypes, by = "celltype") %>%
      replace_na(list(bulkid = b_id, prp = 0))
    proportion_full <- rbind(proportion_full, prop)
  }
  proportion_full <- proportion_full %>%
    dplyr::select(-cse_sum, -bulk_sum)

  return(list(proportions = proportion_full, cse = cse_output))
}

# this function reads the proportions pre-determined by CIBERSORT
read_cibersort_proportions <- function(cibersortfile, logfc, pctchange, irun) {
  # this is the folder structure defined by the bash script that runs the cibersort/hires DOCKER container
  input_folder <- paste0("HiRes_fc_", logfc, "_pct_", pctchange, "_run_", irun, "_", suffix)
  cibersort_prp <- read.table(
    file.path("cibersort", input_folder, paste0(cibersortfile)),
    sep = "\t",
    header = TRUE
  )
  cibersort_proportions <- cibersort_prp %>%
    as_tibble() %>%
    dplyr::select(-P.value, -Correlation, -RMSE) %>%
    dplyr::rename(bulkid = Mixture) %>%
    pivot_longer(!bulkid, names_to = "celltype", values_to = "prp")
  return(cibersort_proportions)
}

prepare_islet_input <- function(bulk_counts, proportions) {
  bulk_ids <- bulk_counts %>% pull(bulkid)
  bulk_mat <- bulk_counts %>%
    dplyr::select(starts_with("ENSG")) %>%
    # this is necessary in order to remove genes with constant expression
    # as ISLET can not deal with such genes
    dplyr::select(where(~ n_distinct(.) > 1)) %>%
    as.matrix()
  bulk_mat <- data.frame(bulk_mat %>% t())
  colnames(bulk_mat) <- bulk_ids %>% as.character()

  grouping <- bulk_counts %>% dplyr::select(bulkid, group)

  # remove celltypes which are 0 expressed more than 90 perc of cells
  # as ISLET can not deal with proportions that are 0 accross almost all samples
  proportions <- proportions %>%
    group_by(celltype) %>%
    filter(across(prp, ~sum( . == 0)) < 90) %>%
    ungroup()

  ct_proportions <- proportions %>%
    pivot_wider(names_from = "celltype", values_from = "prp")
  colData <- ct_proportions %>%
    inner_join(grouping, by = "bulkid") %>%
    dplyr::rename(sample_ID = bulkid) %>%
    relocate(group)

  colData <- data.frame(
    colData,
    row.names = bulk_ids
  )

  islet_se <- SummarizedExperiment(assays = list(counts = bulk_mat), colData = colData)
  islet_input <- dataPrep(dat_se = islet_se)

  return(list(islet_input = islet_input, islet_se = islet_se))
}

# here we provide the option to let run islet with actual ground truth proportions
# or (as intendend) with CIBERSORT proportions, i.e. use_true_proportions = FALSE
fit_bulks_islet <- function(train_counts, bulks_mut, bulks_nonmut, sclong, logfc, pctchange, ngenes, irun, genelist, use_true_proportions) {
  # data preparation
  bulks_case_ctrl <- rbind(
      bulks_mut$bulks %>% add_column(group = "case"),
      bulks_nonmut$bulks %>% add_column(group = "ctrl")
  )
  bulks_case_ctrl <- bulks_case_ctrl %>% mutate(bulkid = 1:(2 * nbulks))
  composition_case_ctrl <- rbind(
      bulks_mut$composition %>% add_column(group = "case"),
      bulks_nonmut$composition %>% add_column(group = "ctrl") %>% mutate(bulkid = bulkid + 50)
  )
  bulks_case_ctrl <- list(bulks = bulks_case_ctrl, composition = composition_case_ctrl)
  if(use_true_proportions) {
    proportions <- true_quantities(bulks_case_ctrl, sclong)$proportions
  } else {
    proportions <- read_cibersort_proportions("CIBERSORTxGEP_NA_Fractions-Adjusted.txt", logfc, pctchange, irun) %>% 
      pivot_wider(names_from = "celltype", values_from = "prp") %>%
      mutate(bulkid = 1:(2*nbulks)) %>%
      pivot_longer(-bulkid, names_to = "celltype", values_to = "prp")
  }

  islet_input <- prepare_islet_input(bulks_case_ctrl$bulks, proportions)$islet_input

  # here we just reproduce parts of the data prepeartion as the res.test
  # df does not store gene names
  bulk_mat <- bulks_case_ctrl$bulks %>%
    dplyr::select(starts_with("ENSG")) %>%
    dplyr::select(where(~ n_distinct(.) > 1)) %>%
    as.matrix()
  bulk_mat <- data.frame(bulk_mat %>% t())
  available_genes <- rownames(bulk_mat)

  res.test <- isletTest(input = islet_input, BPPARAM = MulticoreParam(workers = 50))
  rownames(res.test) <- available_genes

  return(res.test)
}

benchmark_islet <- function(train_counts, bulks_mut, bulks_nonmut, sclong, logfc, pctchange, ngenes, irun, genelist, use_true_proportions) {
  cse_all <- tibble()
  
  # here we use islet's custom function to determine p-values for differential expression
  tryCatch(
    {
      res.test <- fit_bulks_islet(train_counts, bulks_mut, bulks_nonmut, sclong, logfc, pctchange, ngenes, irun, genelist, use_true_proportions)
      for (ctype in colnames(res.test)) {
        computed_genelist <- res.test %>%
          as_tibble(rownames = "gene") %>%
          dplyr::select(gene, all_of(ctype)) %>%
          mutate(score = 1 - !!as.symbol(ctype)) %>%
          mutate(mutated = gene %in% genelist)
        precrec_obj <- evalmod(scores = computed_genelist$score, labels = computed_genelist$mutated)
        if (use_true_proportions) {
          ggsave(paste0("data/plots/roc_cse_islet_sim_", suffix, "_", ngenes, "genes_", logfc, "lfc_", ctype, "_run_", irun, ".pdf"), autoplot(precrec_obj))
          aucprc <- attr(precrec_obj, "aucs")
          cse_all <- rbind(cse_all, tibble(aucprc) %>% add_column(celltype = ctype, ereg = "expression") %>% add_column(algo = "islet"))
        } else {
          ggsave(paste0("data/plots/roc_cse_islet_cibersort_sim_", suffix, "_", ngenes, "genes_", logfc, "lfc_", ctype, "_run_", irun, ".pdf"), autoplot(precrec_obj))
          aucprc <- attr(precrec_obj, "aucs")
          cse_all <- rbind(cse_all, tibble(aucprc) %>% add_column(celltype = ctype, ereg = "expression") %>% add_column(algo = "islet_cibersort"))
        }
      }
    },
    error = function(e) {
      print(paste0("Singular Matrix in LU factorization for some celltypes"))
    }
  )
  return(cse_all)
}

# CIBERSORTx HiRes

benchmark_cibersort_hires <- function(train_counts, logfc, pctchange, ngenes, irun, genelist) {
  cse_all <- tibble()
  celltypes <- train_counts %>% pull(celltype) %>% unique()
  input_folder <- paste0("HiRes_fc_", logfc, "_pct_", pctchange, "_run_", irun, "_", suffix)
  for (ctype in celltypes) {
    # these are the matrices storing the cell type-specific expression in CIBERSORT
    cibersort_cse <- read.table(
      file.path("cibersort", input_folder, paste0("CIBERSORTxHiRes_NA_", ctype,"_Window36.txt")),
      sep = "\t",
      header = TRUE,
      as.is = TRUE,
      row.names = 1
    ) %>% as.matrix()
    cibersort_cse <- cibersort_cse[!rowSums(!is.finite(cibersort_cse)), ]
    # first 50 bulks are the mutated, last 50 the non-mutated
    sampleclasses <- c(rep(TRUE, 50), rep(FALSE, 50)) %>% as.factor()
    if (nrow(cibersort_cse) > 2) { # dispersion cannot be fitted from only two bulks.
      tryCatch(
        {
          edgerobj <- DGEList(counts = cibersort_cse, group = sampleclasses)
          mm <- model.matrix(~group, data = edgerobj$samples)
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
          ggsave(paste0("data/plots/roc_cse_cibersort_hires_sim_", suffix, "_", ngenes, "genes_", logfc, "lfc_", ctype, "_run_", irun, ".pdf"), autoplot(precrec_obj))
          aucprc <- attr(precrec_obj, "aucs")
          cse_all <- rbind(cse_all, tibble(aucprc) %>% add_column(celltype = ctype, ereg = "expression") %>% add_column(algo = "cibersort_hires"))
        },
        error = function(e) {
          print(paste0("There is no variance in the celltype specific expression for the celltype", ctype))
        }
      )
    }
  }
  return(cse_all)
}

# ACTUAL SIMULATION
cse_all <- tibble()
# stores both sc_library as well as bulk source separated via "train/test"
sclong <- readRDS(paste0("data/sclong_", suffix, ".rds"))

# We perform the  simulation along fixed "dimensions" simultaneously:
# pctgenes = 17 is fixed (38/255 genes) and we vary the logfc
# logfc = 0.8 is fixed and we vary the percentage of genes
for (logfc in log2foldchanges) {
  for (pctchange in pct_genes_to_change) {
    # stores the actually modified genes
    genelist <- readRDS(paste0("data/genelist_", pctchange, "_", suffix, ".rds"))
    ngenes <- length(genelist)
    if (logfc == 0.8 || pctchange == 17) {
      # stores the reference single cell library
      train_counts <- readRDS(paste0("data/train_counts_", pctchange, "_", logfc, "_", suffix, ".rds"))

      if ("bp_subtypes" %in% algos) {
          message("computing clustering")
          clustering <- compute_subclustering(train_counts)
      }
      for (irun in 1:nruns) {
        print(paste0("run ", irun, " / ", nruns, " logfc = ", logfc, ", ngenes = ", ngenes, " (", pctchange, " % of all)"))

        # simulated bulks
        bulks_mut <- readRDS(paste0("data/bulks_mut_", pctchange, "_", logfc, "_run_", irun, "_", suffix, ".rds"))
        bulks_nonmut <- readRDS(paste0("data/bulks_nonmut_", pctchange, "_", logfc, "_run_", irun, "_", suffix, ".rds"))

        if ("bp_nosubtypes" %in% algos) {
          cse_all <- rbind(cse_all, benchmark_bayesPrism(train_counts, bulks_mut, bulks_nonmut, NULL, logfc, ngenes, "nosub", irun, genelist) %>%
            add_column(run = irun, logfc = logfc, ngenes = ngenes, pct_genes = pctchange))
        }
        if ("bp_subtypes" %in% algos) {
          message("fitting")
          cse_all <- rbind(cse_all, benchmark_bayesPrism(train_counts, bulks_mut, bulks_nonmut, clustering, logfc, ngenes, "sub", irun, genelist) %>%
            add_column(run = irun, logfc = logfc, ngenes = ngenes, pct_genes = pctchange))
        }
        if ("tr" %in% algos) {
          cse_all <- rbind(
                      cse_all,
                            benchmark_tr(train_counts, bulks_mut, bulks_nonmut, logfc, ngenes, irun, genelist) %>% 
                              add_column(run = irun, logfc = logfc, ngenes = ngenes, pct_genes = pctchange))
        }
        if ("bmind" %in% algos) {
            cse_all <- rbind(
                    cse_all,
                    benchmark_bmind(train_counts, sclong, bulks_mut, bulks_nonmut, logfc, ngenes, irun, genelist) %>%
                      add_column(run = irun, logfc = logfc, ngenes = ngenes, pct_genes = pctchange))
        }
        if ("islet" %in% algos) {
            use_true_proportions <- TRUE
            cse_all <- rbind(
                    cse_all,
                    benchmark_islet(train_counts, bulks_mut, bulks_nonmut, sclong, logfc, pctchange, ngenes, irun, genelist, use_true_proportions) %>%
                      add_column(run = irun, logfc = logfc, ngenes = ngenes, pct_genes = pctchange))
        }
        if ("islet_cibersort" %in% algos) {
            use_true_proportions <- FALSE
            cse_all <- rbind(
                    cse_all,
                    benchmark_islet(train_counts, bulks_mut, bulks_nonmut, sclong, logfc, pctchange, ngenes, irun, genelist, use_true_proportions) %>%
                      add_column(run = irun, logfc = logfc, ngenes = ngenes, pct_genes = pctchange))
        }
        if ("cibersort_hires" %in% algos) {
          cse_all <- rbind(
              cse_all,
              benchmark_cibersort_hires(train_counts, logfc, pctchange, ngenes, irun, genelist) %>%
              add_column(run = irun, logfc = logfc, ngenes = ngenes, pct_genes = pctchange))
        }
      }
      saveRDS(cse_all, file = paste0("data/plots/sim_roc_fc_", logfc, "_pg_", pctchange, suffix, ".rds"))
    }
  }
}
# this file stores all AUC/ROC values reported in the publication
saveRDS(cse_all, file = paste0("data/plots/sim_roc_", suffix, ".rds"))