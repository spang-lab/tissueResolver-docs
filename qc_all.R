# script to generate quality score plots for whole set of genes
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
  message(paste0("No file grouping.csv in data folder. Please compute it in tissue_resolver_vignette.Rmd or download it."))
}
groupA <- "GCB"
groupB <- "ABC"
grouping <- grouping %>% mutate(group = factor(group, levels = c(groupA, groupB)))


#################################################
################## plotting #####################
#################################################


# genewise mean relative residuals vs. avarage bootstrap variance
p_gene_mean_var_all <- qc$genes %>% ggplot(aes(x = relres, y = relres_mean_var)) +
  geom_point() +
  geom_point() +
  ggtitle("a") +
  xlab(TeX("$g_g$")) +
  ylab(TeX("$v_g$")) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 20, face = "bold")
  )

# mean relative residuals vs. average bulk_expression
qc_genes_plot_all <- qc$genes %>% mutate(uselabel = ifelse((log1p(bulk_expression) > 2.5 & relres > 0.65) | log1p(bulk_expression) > 4 & relres > 0.3, hgnc_symbol, ""))
p_qc_log_all <- ggplot(qc_genes_plot_all, aes(x = relres, y = log1p(bulk_expression), label = uselabel)) +
  geom_point() +
  geom_text_repel() +
  ggtitle("b") +
  xlab(TeX("$g_g$")) +
  ylab(TeX("$\\log(Y_g+1)$")) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 20, face = "bold")
  )

# bulkwise mean relative residuals vs. bootstrap variance
qc_bulks_plot_all <- qc$bulks %>% inner_join(grouping, by = "bulk_id")
p_bulks_mean_var <- qc_bulks_plot_all %>% ggplot(aes(x = relres, y = relres_mean_var, colour = group)) +
  geom_point() +
  ggtitle("c") +
  xlab(TeX("$b_s$")) +
  ylab(TeX("$v_s$")) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 20, face = "bold")
  ) +
  scale_color_manual(name = "Genetic sub-type", labels = c("GCB", "ABC", "Unclass"), values = c("red", "blue", "grey"))
  
ggsave(file.path(plotdir, "relres_all.pdf"), arrangeGrob(p_gene_mean_var_all, p_qc_log_all, p_bulks_mean_var, ncol = 3), width = 30, height = 10)