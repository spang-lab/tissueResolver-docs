# script to generate quality score plots for stromal genes
library(tissueResolver)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(latex2exp)

# function to map gene ENSEMBL names to HGNC symbols
source("map_feature_names.R")


data_dir <- "data"
plotdir <- file.path(data_dir, "plots")
dir.create(plotdir, showWarnings = FALSE)


#################################################
############## data preparation ################
#################################################

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

# the genes from the stromal signature
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


#################################################
######## Retrieve cluster specific csre #########
#################################################

# filter cell cluster specific expression
single_cluster <- "17"
message(paste0("Filtering from csre this may take a bit ..."))
csre_type_expr <- csre %>%
  filter(celltype == single_cluster) %>%
  ungroup() %>%
  dplyr::select(bulk_id, gene, expression) %>%
  drop_na()
message(paste0("Done filtering."))

# compute mean expression over bulks
csre_type_mean_expr <- csre_type_expr %>%
  ungroup() %>%
  group_by(gene) %>%
  mutate(mean_expr = mean(expression)) %>%
  dplyr::select(-bulk_id, -expression) %>%
  unique()



#################################################
################## plotting #####################
#################################################

# genewise mean relative residuals vs. average bootstrap variance
p_gene_mean_var_stromal <- qc$genes %>%
  filter(gene %in% additional_genes) %>%
  ggplot(aes(x = relres, y = relres_mean_var, color = '#19c719')) +
  geom_point() +
  ggtitle("a") +
  xlab(TeX("$g_g$")) +
  ylab(TeX("$v_g$")) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 20)
  )+
  scale_color_manual(name = "Genes", labels = c("stromal", "other"), values = c("#19c719", "black"))+
  scale_fill_discrete(guide = "none")



# genewise mean expression in cluster 17 vs. genewise overall mean expression
qc_log_stromal <- qc$genes %>% inner_join(csre_type_mean_expr, by = "gene")

# label simultaneously highly expressed genes
qc_log_stromal$label <- NA
qc_log_stromal$label[log1p(qc_log_stromal$mean_expr) > 2] <- 
                  qc_log_stromal$hgnc_symbol[log1p(qc_log_stromal$mean_expr) > 2]

p_qc_log_stromal <- qc_log_stromal %>%
  ggplot(aes(x = log1p(bulk_expression), y = log1p(mean_expr), label = label, color = ifelse(gene %in% additional_genes,'#19c719','black'))) +
  geom_point() +
  geom_text_repel() +
  ggtitle("b") +
  ylab(TeX("$\\log(\\tilde{Y}^{17}_g+1)$")) +
  xlab(TeX("$\\log(Y_g+1)$")) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 20)
  ) +
  scale_color_manual(name = "Genes", labels = c("stromal", "other"), values = c("#19c719", "black"))+
  scale_fill_discrete(guide = "none")

p_qc_log_stromal_res <- qc_log_stromal %>%
  ggplot(aes(x = relres, y = log1p(mean_expr), label = label, color = ifelse(gene %in% additional_genes,'#19c719','black'))) +
  geom_point() +
  geom_text_repel() +
  ggtitle("b") +
  xlab(TeX("$\\g_g$")) +
  ylab(TeX("$\\log(\\tilde{Y}^{17}_g+1)$")) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 20)
  ) +
  scale_color_manual(name = "Genes", labels = c("stromal", "other"), values = c("#19c719", "black"))+
  scale_fill_discrete(guide = "none")

ggsave(file.path(plotdir, "stromal_res_new.pdf"), arrangeGrob(p_gene_mean_var_stromal, p_qc_log_stromal_res, p_qc_log_stromal, ncol = 3), width = 30, height = 10)


