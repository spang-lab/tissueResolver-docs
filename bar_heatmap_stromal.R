# script to generate the bar/heatmap plot for stromal genes
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


# read grouping
if (file.exists(file.path(data_dir, "grouping.csv"))) {
  grouping <- as_tibble(read.csv(file.path(data_dir, "grouping.csv")))
} else {
  stop(paste0("No file grouping.csv in data folder. Please compute it in tissue_resolver_vignette.Rmd or download it."))
}
groupA <- "GCB"
groupB <- "ABC"
grouping <- grouping %>% mutate(group = factor(group, levels = c(groupA, groupB)))

# consider only stromal signature genes for comparison in csre
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
csre <- csre %>% filter(gene %in% (additional_genes))

# map csre genes to HGNC symbols
mapping <- map_feature_names(
  csre %>% pull(gene),
  input.type = "ensembl_gene_id",
  output.type = "hgnc_symbol",
  n.tries = 10,
  undo.safety.check = TRUE
)
mapvec <- mapping[["hgnc_symbol"]]
names(mapvec) <- mapping[["ensembl_gene_id"]]
csre <- csre %>% rename(ensembl=gene) %>% mutate(gene = mapvec[ensembl])

# join bulk groups (ABC vs. GCB)
csre <- csre %>% inner_join(grouping %>% dplyr::select(bulk_id, group)) %>% drop_na()

#################################################
################## plotting ####################
#################################################


p <- plot_csre(csre, heatmapviz = "relchange", barplotviz = "relative", groupA = groupA, groupB = groupB, ctypes = NULL, addexpressionlevel = FALSE)
ggsave(file.path(plotdir, "bar_heatmap_stromal.pdf"),
  p,
  width = 12,
  height = 12,
  dpi = 150,
  units = "in"
)
