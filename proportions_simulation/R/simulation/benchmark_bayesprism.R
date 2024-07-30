# BayesPrism

# this subclustering refines cell types into cell states in order
# to provide BayesPrism with more fine grained info
compute_subclustering <- function(train_counts, cluster_resolution = 0.5) {
    clustering <- train_counts %>% pull(cell_id)
    names(clustering) <- clustering
    for (cl in levels(factor(train_counts %>% pull(celltype)))) {
        message(paste0("clustering for cluster ", cl))
        exprmat <- train_counts %>%
            filter(celltype == cl) %>%
            dplyr::select(starts_with("ENSG")) %>%
            as.matrix()
        rownames(exprmat) <- train_counts %>%
            filter(celltype == cl) %>%
            pull(cell_id)
        # too few samples
        if (dim(exprmat)[1] < 2) {
            cellids <- colnames(seurobj)
            clustering[cellids] <- cl
            next
        }
        seurobj <- CreateSeuratObject(counts = t(exprmat), min.cells = 3, min.features = 5)
        # too small cluster size
        if (dim(seurobj)[2] < 100) {
            cellids <- colnames(seurobj)
            clustering[cellids] <- cl
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
        clusternames <- paste0(cl, "_", clusters %>% pull(seurat_clusters))
        clustering[cellids] <- clusternames
    }
    return(clustering)
}

# fit bulks using BayesPrism
fit_bulks_bayesPrism <- function(bulks, sccounts, clustering, test_samples) {
    # deconvolute only test_samples
    bulks <- bulks %>% filter(bulk_id %in% test_samples)
    bk.dat <- bulks %>%
        dplyr::select(-bulk_id) %>%
        as.matrix()
    rownames(bk.dat) <- bulks %>% pull(bulk_id)
    sc.dat <- sccounts %>%
        dplyr::select(starts_with("ENSG")) %>%
        as.matrix()
    rownames(sc.dat) <- sccounts %>% pull(cell_id)

    cell.type.labels <- sccounts %>% pull(celltype)
    names(cell.type.labels) <- sccounts %>% pull(cell_id)

    if (is.null(clustering)) {
        cell.state.labels <- cell.type.labels
    } else {
        cell.state.labels <- clustering[sccounts %>% pull(cell_id)]
    }


    myPrism <- new.prism(
        reference = sc.dat,
        mixture = bk.dat,
        input.type = "count.matrix",
        cell.type.labels = cell.type.labels,
        cell.state.labels = cell.state.labels,
        key = NULL,
        outlier.cut = 0.01,
        outlier.fraction = 0.1
    )

    message("Run bayesPrism")
    bp.res <- run.prism(prism = myPrism, n.cores = 50)
    return(bp.res)
}

get_bayesPrism_proportions <- function(bp.res, test_samples) {
    # compute and reshape cell proportions
    cell_proportions <- get.fraction(
        bp = bp.res,
        which.theta = "final",
        state.or.type = "type"
    )
    # just reshape the proportions to have format
    # consistent with other algorithms
    cell_proportions <- cell_proportions %>%
        as_tibble(rownames = "bulk_id") %>%
        pivot_longer(!bulk_id, names_to = "celltype", values_to = "prp")
    return(cell_proportions)
}

get_bayesPrism_cse <- function(bp.res) {
    # this is where the celltypes are stored in a prism object
    celltypes <- names(bp.res@prism@map)
    cse <- tibble()
    for (ct in celltypes) {
        cse_ct <- get.exp(
            bp = bp.res,
            state.or.type = "type",
            cell.name = ct
        )
        cse_ct <- cse_ct %>%
            as_tibble(rownames = "bulk_id") %>%
            add_column(celltype = ct)
        cse <- rbind(cse, cse_ct)
    }
    cse <- cse %>%
        pivot_longer(cols = starts_with("ENSG"), names_to = "gene", values_to = "cse") %>%
        pivot_wider(names_from = "celltype", values_from = "cse")
    
    return(cse)
}