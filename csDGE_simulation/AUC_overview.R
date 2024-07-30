# reporduces the AUC curves in the publication

library(tidyverse)
library(ggplot2)
library(gridExtra)

log2foldchanges <- c(0.5, 0.8, 1.0, 1.5, 2)
pct_genes_to_change <- c(2, 5, 10, 17, 20, 30, 40, 50, 60, 70, 80, 90, 95)
# suffix <- "cd8"
suffix <- "B"

algos <- c("tr", "bp_sub", "bp_nosub", "bmind", "islet_cibersort", "cibersort_hires")

if (suffix == "cd8") {
    rds <- readRDS("data/sim_roc_cd8.rds") %>%
        filter(!(ngenes %in% c(29, 34))) %>% 
        filter(algo %in% algos) %>% 
        mutate(algo = ifelse(algo == "islet_cibersort", "islet", algo))
}

if (suffix == "B") {
    rds <- readRDS("data/sim_roc_B.rds") %>% 
        filter(algo %in% algos) %>% 
        mutate(algo = ifelse(algo == "islet_cibersort", "islet", algo))
}


rds_auc_mean <- rds %>%
    filter(ereg == "expression", curvetypes == "ROC") %>%
    group_by(celltype, logfc, ngenes, algo) %>%
    mutate(mean_auc = mean(aucs)) %>%
    mutate(sd_auc = sd(aucs)) %>%
    ungroup() %>%
    dplyr::select(celltype, algo, logfc, ngenes, mean_auc, sd_auc) %>%
    unique()

celltypes <- rds_auc_mean %>%
    pull(celltype) %>%
    unique()

if (suffix == "cd8") {
    rds_auc_mean_logfc <- rds_auc_mean %>% filter(ngenes == 38)
    celltypes <- c("T_CD8", "T_CD4", "B")
}

if (suffix == "B") {
    rds_auc_mean_logfc <- rds_auc_mean %>% filter(ngenes == 39)
    celltypes <- c( "B", "T_CD8", "T_CD4")
}

aucs <- list()

for (ctype in celltypes) {
    rds_aucs <- rds_auc_mean_logfc %>%
        filter(celltype == ctype)

    aucs[[ctype]] <- ggplot(data = rds_aucs, aes(x = logfc, y = mean_auc, colour = algo)) +
        geom_line() +
        geom_point() +
        geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc), width=.02) +
        # scale_y_continuous(limits = c(0.45, 1.0)) +
        coord_cartesian(ylim = c(0.5, 1.0)) +
        ggtitle(paste0(ctype))

}

ggsave(paste0("data/log_fc_auc_", suffix,".pdf"), grid.arrange(grobs = aucs, nrow = 1), width = 20, height = 2)

rds_auc_mean_ngenes <- rds_auc_mean %>% filter(logfc == 0.8) %>% filter(ngenes < 100)

aucs <- list()

for (ctype in celltypes) {
    rds_aucs <- rds_auc_mean_ngenes %>%
        filter(celltype == ctype)

    aucs[[ctype]] <- ggplot(data = rds_aucs, aes(x = ngenes, y = mean_auc, colour = algo)) +
        geom_line() +
        geom_point() +
        geom_errorbar(aes(ymin = mean_auc - sd_auc, ymax = mean_auc + sd_auc), width = 1) +
        # scale_y_continuous(limits = c(0.45, 1.0)) +
        coord_cartesian(ylim = c(0.5, 1.0)) +
        ggtitle(paste0(ctype))
}
ggsave(paste0("data/ngenes_auc_", suffix,".pdf"), grid.arrange(grobs = aucs, ncol = 1), width = 13, height = 7)
