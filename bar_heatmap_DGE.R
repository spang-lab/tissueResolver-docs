# script to generate the bar/heatmap plot for differentially expressed genes
library(tissueResolver)
library(tidyverse)
library(edgeR)
library(ggplot2)
library(ggrepel)

# function to map gene ENSEMBL names to HGNC symbols
source("map_feature_names.R")


data_dir <- "data"
plotdir <- file.path(data_dir, "plots")
dir.create(plotdir, showWarnings = FALSE)


#################################################
############## data preparation ################
#################################################

# read bulks
if (!file.exists(file.path(data_dir, "bulks.rds"))) {
  stop("Bulk data missing. Please provide bulks.rds in your data folder")
} else {
  bulks <- readRDS(file.path(data_dir, "bulks.rds"))
}


# read CSRE
if (file.exists(file.path(data_dir, "csre.rds"))) {
  message(paste0("Reading csre from file csre.rds this may take a bit..."))
  csre <- readRDS(file.path(data_dir, "csre.rds"))
  message(paste0("Reading csre.rds done."))
} else {
  stop(paste0("No file csre.rds in data folder. Please compute it in tissue_resolver_vignette.Rmd or download it."))
}

# read qc
if (file.exists(file.path(data_dir, "qc.rds"))) {
  qc <- readRDS(file.path(data_dir, "qc.rds"))
} else {
  stop(paste0("No file qc.rds in data folder. Please compute it in tissue_resolver_vignette.Rmd"))
}

# join the HGNC symbols to qc for nice labeling
mapping <- map_feature_names(
  qc$genes %>% pull(gene),
  input.type = "ensembl_gene_id",
  output.type = "hgnc_symbol",
  n.tries = 10,
  undo.safety.check = TRUE
)
qc$genes <- qc$genes %>% inner_join(tibble(mapping) %>% rename(gene = ensembl_gene_id), by = "gene")


# read grouping
if (file.exists(file.path(data_dir, "grouping.csv"))) {
  grouping <- as_tibble(read.csv(file.path(data_dir, "grouping.csv")))
} else {
  stop(paste0("No file grouping.csv in data folder. Please compute it in tissue_resolver_vignette.Rmd or download it."))
}
groupA <- "GCB"
groupB <- "ABC"
grouping <- grouping %>% mutate(group = factor(group, levels = c(groupA, groupB)))


# consider only bulk samples which are not unclassified
actual_grouping <- grouping %>% drop_na()

#################################################
################## edgeR ########################
#################################################

# determine differentially expressed genes
edgerobj <- DGEList(counts=bulks$bulk.counts[,actual_grouping[["bulk_id"]]], group=actual_grouping[["group"]])
keep <- filterByExpr(y = edgerobj, min.count = 0.01*median(edgerobj$counts))
edgerobj <- edgerobj[keep, , keep.lib.sizes=FALSE]
edgerobj <- calcNormFactors(edgerobj)
edgerobj <- estimateDisp(edgerobj)
et <- exactTest(edgerobj)
tab <- topTags(et, n="Inf")
degs <- as_tibble(tab$table) %>% add_column(gene=rownames(tab$table))

# use quality control metrics for subsequent filtering
interesting_degs <- degs %>% inner_join(qc$genes, by="gene")

#################################################
################## filtering DEGs ###############
#################################################


fdr <- 0.05
message("filter by bulk deg fdr")
interesting_degs <- interesting_degs %>%
  filter(FDR < fdr)
message(paste0(nrow(interesting_degs), " genes remaining."))

maxrelres <- 0.5
message("filter by maximum relative residual")
interesting_degs <- interesting_degs %>% filter(relres <= maxrelres)
message(paste0(nrow(interesting_degs), " genes remaining."))

maxrelresmeanvar <- 0.03
message("filter by maximum bootstrap variance in residual")
interesting_degs <- interesting_degs %>% filter(relres_mean_var <= maxrelresmeanvar)
message(paste0(nrow(interesting_degs), " genes remaining."))


minabsfc <- 0.8
message("filter by fold change")
interesting_degs <- interesting_degs %>% filter(abs(logFC) > minabsfc)
message(paste0(nrow(interesting_degs), " genes remaining."))


message("sort genes by csre log ratio expression variance")
genes_sorted_by_csre_variance <- csre %>%
  inner_join(grouping %>% dplyr::select(bulk_id, group), by="bulk_id") %>%
  group_by(gene, celltype, group) %>% mutate(expr = mean(expression)) %>%
  dplyr::select(-bulk_id, -regulation, -expression) %>% ungroup() %>% unique() %>%
  pivot_wider(values_from=expr, names_from=group) %>%
  rename(groupA = groupA) %>% rename(groupB = groupB) %>%
  mutate(logratio=log((groupB+1)/(groupA+1))) %>%
  dplyr::select(celltype, gene, logratio) %>%
  group_by(gene) %>% mutate(varex=var(logratio)) %>% ungroup() %>%
  dplyr::select(-logratio, -celltype) %>% unique() %>% arrange(desc(varex)) %>% dplyr::select(gene)
interesting_degs <- genes_sorted_by_csre_variance %>% right_join(interesting_degs)

# save filtered DEGs
ngenes <- interesting_degs %>% nrow()
interesting_genes <- interesting_degs %>% pull(gene)

# filter DEGs from csre
csre <- csre %>% filter(gene %in% (interesting_genes))

# map csre genes to HGNC
mapvec <- mapping[["hgnc_symbol"]]
names(mapvec) <- mapping[["ensembl_gene_id"]]
csre <- csre %>% rename(ensembl=gene) %>% mutate(gene = mapvec[ensembl])

# join bulk groups
csre <- csre %>% inner_join(grouping %>% dplyr::select(bulk_id, group)) %>% drop_na()

#################################################
################## plotting ####################
#################################################


p <- plot_csre(csre, heatmapviz = "relchange", barplotviz = "relativeA", groupA = groupA, groupB = groupB, ctypes = NULL, addexpressionlevel = FALSE)
ggsave(file.path(plotdir, "bar_heatmap_DGE.pdf"),
  p,
  width = 12,
  height = 12,
  dpi = 150,
  units = "in"
)
