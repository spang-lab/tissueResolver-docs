# script to generate volcano plots for cluster 17 and 18
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


# read grouping
if (file.exists(file.path(data_dir, "grouping.csv"))) {
  grouping <- as_tibble(read.csv(file.path(data_dir, "grouping.csv")))
} else {
  stop(paste0("No file grouping.csv in data folder. Please compute it in tissue_resolver_vignette.Rmd or download it."))
}
groupA <- "GCB"
groupB <- "ABC"
grouping <- grouping %>% mutate(group = factor(group, levels = c(groupA, groupB)))

# convert csre tibble to matrix to be processed by edgeR
csre_to_matrix <- function(csdf, colname = "regulation", genename = "gene") {
  edf <- csdf %>% dplyr::select(-group) %>% pivot_wider(names_from=any_of(genename), values_from=any_of(colname))
  bulknames <- edf %>% pull(bulk_id)
  expr <- edf %>% dplyr::select(-bulk_id) %>% as.matrix()
  rownames(expr) <- bulknames
  return(t(expr))
}

# CLUSTER 17

#################################################
######## Retrieve cluster specific csre #########
#################################################


single_cluster <- "17"
message(paste0("Filtering single cluster regulation from csre this may take a bit ..."))
csre_cluster <- csre %>% filter(celltype == single_cluster) %>% ungroup() %>% dplyr::select(bulk_id, gene, regulation) %>% drop_na()
message(paste0("Done filtering."))


# consider only those bulk samples which intersect with the bulk samples of csre_cluster and which are not unclassified
actual_grouping <- grouping %>% drop_na() %>% filter(bulk_id %in% csre_cluster[["bulk_id"]])
csre_cluster_actual_bulks <- csre_cluster %>% filter(bulk_id %in% actual_grouping[["bulk_id"]]) %>% inner_join(actual_grouping, by = "bulk_id")
# convert cell cluster specifc regulation to matrix for edgeR
csre_matrix <- csre_to_matrix(csre_cluster_actual_bulks, colname = "regulation", genename = "gene")
# order columns consistent with grouping
csre_matrix <- csre_matrix[, actual_grouping[["bulk_id"]]]


#################################################
################## edgeR ########################
#################################################

edgerobj <- DGEList(counts = csre_matrix, group = actual_grouping[["group"]])
edgerobj <- calcNormFactors(edgerobj)
edgerobj <- estimateDisp(edgerobj)
et <- exactTest(edgerobj)
# for false discovery rate
tab <- topTags(et, n = "Inf")
degs <- as_tibble(tab$table) %>% add_column(gene = rownames(tab$table))
# store for each gene quality scores and exactTest results
interesting_degs <- degs %>% inner_join(qc$genes, by = "gene")
ngenes <- interesting_degs %>% nrow()
interesting_genes <- interesting_degs %>% pull(gene)
csre_filter <- csre_cluster_actual_bulks %>% filter(gene %in% (interesting_genes))


#################################################
############# map ENSEMBL to HGNC ################
#################################################


mapping <- map_feature_names(
  interesting_genes,
  input.type = "ensembl_gene_id",
  output.type = "hgnc_symbol",
  n.tries = 10,
  undo.safety.check = TRUE
)
mapvec <- mapping[["hgnc_symbol"]]
names(mapvec) <- mapping[["ensembl_gene_id"]]
# join the HGNC symbols
csre_filter_hgnc <- csre_filter %>%
  rename(ensembl = gene) %>%
  mutate(gene = mapvec[ensembl])

# join bulk groups for plotting
csre_plot <- csre_filter_hgnc %>%
  inner_join(grouping %>% dplyr::select(bulk_id, group)) %>%
  drop_na()


#################################################
################## plotting #####################
#################################################

# consider stromal signature genes for a posteriori labeling
additional_genes <- c(
              "ENSG00000148848",
              "ENSG00000182492",
              "ENSG00000245848",
              "ENSG00000197467",
              "ENSG00000084636",
              "ENSG00000108821",
              "ENSG00000164692",
              "ENSG00000130635",
              "ENSG00000204262",
              "ENSG00000142173",
              "ENSG00000163359",
              "ENSG00000171812",
              "ENSG00000198223",
              "ENSG00000011465",
              "ENSG00000172638",
              "ENSG00000213853",
              "ENSG00000078098",
              "ENSG00000166147",
              "ENSG00000115414",
              "ENSG00000136235",
              "ENSG00000142798",
              "ENSG00000115594",
              "ENSG00000138448",
              "ENSG00000160255",
              "ENSG00000049130",
              "ENSG00000112769",
              "ENSG00000172037",
              "ENSG00000196878",
              "ENSG00000129038",
              "ENSG00000119681",
              "ENSG00000139329",
              "ENSG00000117122",
              "ENSG00000157227",
              "ENSG00000100985",
              "ENSG00000145431",
              "ENSG00000122861",
              "ENSG00000133110",
              "ENSG00000169439",
              "ENSG00000113140",
              "ENSG00000140682",
              "ENSG00000137801",
              "ENSG00000035862",
              "ENSG00000196616",
              "ENSG00000181092",
              "ENSG00000105974",
              "ENSG00000105971",
              "ENSG00000125810",
              "ENSG00000107562",
              "ENSG00000172889",
              "ENSG00000024422",
              "ENSG00000157554",
              "ENSG00000170323",
              "ENSG00000115461",
              "ENSG00000144668",
              "ENSG00000128052",
              "ENSG00000091136",
              "ENSG00000116678",
              "ENSG00000173269",
              "ENSG00000189184",
              "ENSG00000261371",
              "ENSG00000138207",
              "ENSG00000154133",
              "ENSG00000095637",
              "ENSG00000152583",
              "ENSG00000164056",
              "ENSG00000120156",
              "ENSG00000229353",
              "ENSG00000236221",
              "ENSG00000231608",
              "ENSG00000229341",
              "ENSG00000206258",
              "ENSG00000233323",
              "ENSG00000236236",
              "ENSG00000168477",
              "ENSG00000110799"
)


volcano_data <- interesting_degs

# color label up and downregulated genes
volcano_data$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and  fdr < 0.05, consider as significantly upregulated
volcano_data$diffexpressed[volcano_data$logFC > 0.5 & volcano_data$FDR < 0.05] <- "UP"
# if log2Foldchange < -0.5 and fdr < 0.05, consider as significantly downregulated
volcano_data$diffexpressed[volcano_data$logFC < -0.5 & volcano_data$FDR < 0.05] <- "DOWN"

# name label differentially expressed genes if they belong to the stromal signature
volcano_data$de_label <- NA
# filter stromal genes
volcano_data$de_label[volcano_data$gene %in% additional_genes] <- volcano_data$gene[volcano_data$gene %in% additional_genes]
# if stromal genes are sufficiently differentially expressed give HGNC label
volcano_data$de_label[volcano_data$logFC > -0.5 & volcano_data$logFC < 0.5] <- NA
volcano_data$de_label[!is.na(volcano_data$de_label)] <- mapvec[volcano_data$de_label[!is.na(volcano_data$de_label)]]

# volcano plot log fold change vs. -log fdr
volcano_17 <- ggplot(data = volcano_data, aes(x = logFC, y = -log10(FDR), col = diffexpressed, label = de_label)) +
  geom_point() +
  theme_minimal() +
  geom_text(check_overlap = TRUE, nudge_y = 0.1, nudge_x = -0.3, show.legend = FALSE) +
  scale_color_manual(
    labels = c("DOWN", "NO", "UP"),
    values = c("red", "gray", "blue3"), name = "Regulation"
  ) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )


# CLUSTER 18

#################################################
######## Retrieve cluster specific csre #########
#################################################


single_cluster <- "18"
message(paste0("Filtering single cluster regulation from csre this may take a bit ..."))
csre_cluster <- csre %>%
  filter(celltype == single_cluster) %>%
  ungroup() %>%
  dplyr::select(bulk_id, gene, regulation) %>%
  drop_na()
message(paste0("Done filtering."))


# consider only those bulk samples which intersect with the bulk samples of csre_cluster and which are not unclassified
actual_grouping <- grouping %>%
  drop_na() %>%
  filter(bulk_id %in% csre_cluster[["bulk_id"]])
csre_cluster_actual_bulks <- csre_cluster %>%
  filter(bulk_id %in% actual_grouping[["bulk_id"]]) %>%
  inner_join(actual_grouping, by = "bulk_id")
# convert cell cluster specifc regulation to matrix for edgeR
csre_matrix <- csre_to_matrix(csre_cluster_actual_bulks, colname = "regulation", genename = "gene")
# order columns consistent with grouping
csre_matrix <- csre_matrix[, actual_grouping[["bulk_id"]]]


# #################################################
# ################## edgeR ########################
# #################################################

edgerobj <- DGEList(counts = csre_matrix, group = actual_grouping[["group"]])
edgerobj <- calcNormFactors(edgerobj)
edgerobj <- estimateDisp(edgerobj)
et <- exactTest(edgerobj)
# for false discovery rate
tab <- topTags(et, n = "Inf")
degs <- as_tibble(tab$table) %>% add_column(gene = rownames(tab$table))
# store for each gene quality scores and exactTest results
interesting_degs <- degs %>% inner_join(qc$genes, by = "gene")
ngenes <- interesting_degs %>% nrow()
interesting_genes <- interesting_degs %>% pull(gene)
csre_filter <- csre_cluster_actual_bulks %>% filter(gene %in% (interesting_genes))


# #################################################
# ############# map ENSEMBL to HGNC ################
# #################################################


mapping <- map_feature_names(
  interesting_genes,
  input.type = "ensembl_gene_id",
  output.type = "hgnc_symbol",
  n.tries = 10,
  undo.safety.check = TRUE
)
mapvec <- mapping[["hgnc_symbol"]]
names(mapvec) <- mapping[["ensembl_gene_id"]]
# join the HGNC symbols
csre_filter_hgnc <- csre_filter %>%
  rename(ensembl = gene) %>%
  mutate(gene = mapvec[ensembl])

# join bulk groups for plotting
csre_plot <- csre_filter_hgnc %>%
  inner_join(grouping %>% dplyr::select(bulk_id, group)) %>%
  drop_na()


#################################################
################## plotting #####################
#################################################

# consider ABC/GCB signature for a posteriori labeling
additional_genes <- c(
  "ENSG00000143727",
  "ENSG00000156127",
  "ENSG00000171791",
  "ENSG00000118971",
  "ENSG00000213923",
  "ENSG00000160213",
  "ENSG00000121966",
  "ENSG00000156136",
  "ENSG00000196937",
  "ENSG00000033170",
  "ENSG00000125166",
  "ENSG00000125245",
  "ENSG00000236418",
  "ENSG00000225890",
  "ENSG00000232062",
  "ENSG00000228284",
  "ENSG00000206305",
  "ENSG00000196735",
  "ENSG00000282657",
  "ENSG00000211899",
  "ENSG00000137265",
  "ENSG00000143772",
  "ENSG00000170421",
  "ENSG00000135363",
  "ENSG00000277443",
  "ENSG00000185697",
  "ENSG00000083454",
  "ENSG00000137193",
  "ENSG00000102096",
  "ENSG00000115956",
  "ENSG00000196396",
  "ENSG00000162924",
  "ENSG00000122026",
  "ENSG00000155926",
  "ENSG00000079263",
  "ENSG00000269404",
  "ENSG00000128040",
  "ENSG00000196628",
  "ENSG00000198467",
  "ENSG00000035403"
)

# cut off outlier genes with extreme low FDR, for nicer depiction
volcano_data <- interesting_degs %>% filter(-log10(FDR) < 10)

# color label up and downregulated genes
volcano_data$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and  fdr < 0.05, consider as significantly upregulated
volcano_data$diffexpressed[volcano_data$logFC > 0.5 & volcano_data$FDR < 0.05] <- "UP"
# if log2Foldchange < -0.5 and fdr < 0.05, consider as significantly downregulated
volcano_data$diffexpressed[volcano_data$logFC < -0.5 & volcano_data$FDR < 0.05] <- "DOWN"

# name label differentially expressed genes if they belong to the ABC/GCB signature
volcano_data$de_label <- NA
# filter stromal genes
volcano_data$de_label[volcano_data$gene %in% additional_genes] <- volcano_data$gene[volcano_data$gene %in% additional_genes]
# if stromal genes are sufficiently differentially expressed give HGNC label
volcano_data$de_label[volcano_data$logFC > -0.5 & volcano_data$logFC < 0.5 | volcano_data$FDR > 0.05] <- NA
volcano_data$de_label[!is.na(volcano_data$de_label)] <- mapvec[volcano_data$de_label[!is.na(volcano_data$de_label)]]

# volcano plot log fold change vs. -log fdr
volcano_18 <- ggplot(data = volcano_data, aes(x = logFC, y = -log10(FDR), col = diffexpressed, label = de_label)) +
  geom_point() +
  ggtitle("b") +
  theme_minimal() +
  geom_text(check_overlap = TRUE, nudge_y = 0.2, nudge_x = 0, show.legend = FALSE) +
  scale_color_manual(
    labels = c("DOWN", "NO", "UP"),
    values = c("red", "gray", "blue3"), name = "Regulation"
  ) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )




ggsave(file.path(plotdir, paste0("volcano_cluster_17_18.pdf")), arrangeGrob(volcano_17, volcano_18, ncol = 2), width = 20, height = 10)

