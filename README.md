# tissueResolver Docs

`tissueResolver` is a package for converting bulk RNA-seq datasets into virtual tissues using information from similar single cell datasets by assigning weights to true single cells, maintaining their molecular integrity. Virtual tissues can be analyzed in a similar way as conventional single cell datasets.

<div style="text-align: center">
    <img src="schematics.png" width=500 alt="tR Pipeline" />
</div>

For theoretical details and interpretation of our results, see our preprint [bioarxiv](https://www.biorxiv.org/content/10.1101/2023.11.15.567357v1).

In this vignette we want to give a general overview over the most important functions
used when building a practical `tissueResolver` pipeline and provide reproducability of the results of our paper.

For detailed explanation and examples of the main package functions of `tissueResolver` 
we refer to the corresponding R package and its documentation available under [Spang-Lab GitHub]( https://github.com/spang-lab/tissueResolver).

# Data Accessibility

We provide all the necessary data for download via the Zenodo DOI: `10.5281/zenodo.10139154`.

The case study of our paper is based on the following raw data:
- **bulk RNA-seq**: From the publication ''Genetics and pathogenesis of diffuse large B-cell lymphoma'' by Schmitz et al. accessible under [GDC](https://gdc.cancer.gov/about-data/publications/DLBCL-2018)
- **single cell RNA-seq**: Here we combined two datasets without further harmonisation:
    - *Dissecting intratumour heterogeneity of nodal B-cell lymphomas at the transcriptional, genetic and drug-response levels* by **[Roider et al., 2020]** accessible under [heiDATA ID VRJUNV](https://heidata.uni-heidelberg.de/dataset.xhtml?persistentId=doi:10.11588/data/VRJUNV) or [EGAS00001004335](https://ega-archive.org/studies/EGAS00001004335)
    - *The landscape of tumor cell states and ecosystems in diffuse large B cell lymphoma* by **[Steen et al., 2021]**, accessibly under [GSE182436](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182436) and [GSE182434](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182434)


We note that the data comprised in `bulks.rds` and `sc.rds` is the result of the gene 
filtering procudure explained in supplemental section **Gene filtering** of our paper.

# Building the .Rmd File

In order to build the file `tissue_resolver_vignette.Rmd` invoke
```{r}
rmarkdown::render("tissue_resolver_vignette.Rmd")
```

If you want to run the code as a script you can simply extract it from the `.Rmd`` file via 

```{r}
knitr::purl("tissue_resolver_vignette.Rmd")
```


We intentionally set some computationally heavy R chunks to `eval = FALSE` and provided the output plots for a quick build.

Setting `eval = TRUE` and only storing `bulks.rds` and `sc.rds` in the `data` folder will compute everything from scratch.

If you run the `tissue_resolver_vignette.Rmd` for the first time, we provided all the data under Zenodo DOI: 10.5281/zenodo.10139154.
For compiling `tissue_resolver_vignette.Rmd` you only need `sc.rds` and `bulks.rds` to be stored in the `data` folder.
The more additional data you provide, the less will be computed by `tissueResolver`.

# Comparative Simulation of BayesPrism and tissueResolver

We additionally provide the script `simulations_paper.R` which allows to reproduce the benchmark simulation of our paper comparing `BayesPrism` with `tissueResolver`,
see section **Simulations** and in particular the pseudo algorithm we provide in the supplemental material of our publication.