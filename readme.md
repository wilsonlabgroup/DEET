
# Differential Expression Enrichment Tool (DEET)

### Dustin Sokolowski: dustin-dot-sokolowski-at-sickkids-dot-ca

### Date: 11/24/2022

![alt text](https://github.com/wilsonlabgroup/DEET/blob/master/vignettes/DEET_V2.png?raw=true)

Design by Sarah Watt

## Description

### DEET

*Description of DEET adapted from manuscript*. The primary use of the differential expression enrichment tool (DEET) is to allow researchers to perform gene set enrichment of their own gene list against the thousands of pairwise DE comparisons stored within DEET. Here, users input a list of genes with an associated p-value and summary statistic (i.e. fold-change), and the DEET_enrich() function will output which pairwise comparisons are associated with the inputted gene list at three levels: overlapping genes, pathways, and transcription factor (TF) target genes (see “Materials and Methods” for details). The gene-set enrichment within DEET_enrich() uses ActivePathways to test for significant overlap between gene sets. Unlike traditional pathway enrichment, the gene lists stored within DEET are also weighted by p-value and fold-change. DEET_enrich() leverages the summary statistics stored within DEET by correlating the fold-changes of overlapping DEGs between the researchers' inputted gene list and each of the enriched comparisons within DEET. This correlation analysis provides further evidence that the inputted gene list shares common biological underpinnings with the overlapping gene list because not only is there a significant overlap in DEGs, but those overlapping DEGs are changing in a similar pattern. Simultaneously, DEET_enrich() uses ActivePathways to enrich the researchers inputted gene list against the same biological pathway and transcription factor gene sets performed on all of the pairwise comparisons within DEET. Finally, DEET_enrich() leverages these enriched pathways and transcription factors to identify which comparisons significantly overlap with the researchers' inputted genes based on their shared biological pathways and transcription factors.
The secondary use of DEET is to interact with the results of the 3142 unique DE comparisons to find commonalities across different experimental designs. Specifically, the DEET_feature_extract() function identifies genes whose p-values drive metadata (e.g., cluster ID, number of DEGs in the study, study source, if TF is enriched etc.). DEET_feature_extract() uses an elastic net regression as well gene associations based on correlation, ANOVA, and Wilcoxon’s test depending on whether the inputted metadata is continuous, categorical, or binomial respectively.

### Shiny App

If you prefer a point-and-click equivalent with >90% of the capability and flexibility of the DEET R package, we encourage you check out our Shiny App!

https://wilsonlab-sickkids-uoft.shinyapps.io/DEET-shiny/

### Manuscript

The manuscript is currently in pre-print: 

https://www.biorxiv.org/content/10.1101/2022.08.29.505468v1.article-metrics

#### Running DEET_enrich

Users inputs a list of genes to identify which DE comparisons within DEET are significantly enriched based on overlapping DE genes, pathways, and TF information.


```{r DEET_enrich, eval=FALSE}

data("example_DEET_enrich_input")
data("DEET_example_data")
DEET_out <- DEET_enrich(example_DEET_enrich_input, DEET_dataset = DEET_example_data)

```

#### Visualization of DEET results

The `proccess_and_plot_DEET_enrich` function uses ggplot2 to generate barplots and dotplots of the User's gene list against the gene sets found within DEET. It uses the direct output of DEET_enrich to plot and relies on ggplot functions.

```{r proccess_and_plot_DEET_enrich, eval=FALSE}

plotting_example <- proccess_and_plot_DEET_enrich(DEET_out, text_angle = 45,
horizontal = TRUE, topn=4)

```

The `DEET_plot_correlation` function uses ggplot2 to generate scatterplots of the coefficients between the users inputted gene list and the log2FoldChanges associated with the DEGs in DEET. If the inputted list did not contain coefficients, then this function is not used.

```{r DEET_plot_correlation, eval=FALSE}

correlation_input <- DEET_out$DE_correlations
correlation_plots <- DEET_plot_correlation(correlation_input)

```


#### Identifying DE genes associated with metadata

The `DEET_feature_extract` uses an elastic net regression, wilcoxon test, ANOVA, and correlation (depending on the response variable) to identify genes that are associated with a particular response (e.g., whether a TF is enriched) across the pairwise comparisons stored within DEET.

```{r DEET_feature_extract, eval=FALSE}

data(DEET_feature_extract_example_matrix)
data(DEET_feature_extract_example_response)
single1 <- DEET_feature_extract(DEET_feature_extract_example_matrix,
DEET_feature_extract_example_response,"categorical")

```

## Downloading the DEG databases within DEET

All processed DEGs, metadata, and enriched pathways in formats compatible with this package as well as other methods such as gene set enrichment analysis are stored here: http://wilsonlab.org/public/DEET_data/

No functions within DEET automatically load data for the user, so the data either needs to be downloaded directly from the ftp, or using the downloader function.

The `DEET_data_download` function, with possible inputs "ALL", "metadata", "enrich", and "feature_extract" automatically downloads the data required to run `DEET_enrich` and/or `DEET_feature_extract`.

We reccomended using:
```{r download_data, eval=FALSE}

downloaded <- DEET_data_download("ALL")
metadata <- downloaded$metadata
DEET_feature_extract_input <- downloaded$DEET_feature_extract
DEET_enrich_input <- downloaded$DEET_enrich

```

Here:
`DEET_enrich_input` replaces `DEET_example_data` for `DEET_enrich()`.
`DEET_feature_extract_input` replaces `DEET_feature_extract_example_matrix` for `DEET_feature_extract()`
Lastly, `metadata` is not directly used in any of the function, but summarizes all of the pairwise comparisons using the following columns.

For the structure of these datatypes, please use `?DEET_data_download`, as they're explained in the `value` column.
`?DEET_example_data` also continuous to break down the structure of the DEET dataset used in `DEET_enrich_input`, if additional detail is required.

Once download, save these data and DEET can be used offline.

## Installation

 DEET relies on the following dependencies which should be downlaoded/updated with DEET automatically. Please ensure that these packages are not open when installing DEET. 

  * ggplot2 - CRAN
  * ActivePathways - CRAN
  * stats - CRAN
  * utils - CRAN
  * downloader - CRAN
  * glmnet - CRAN
  * ggrepel - CRAN
  * dplyr - CRAN
  * pbapply - CRAN

Because all dependencies are in CRAN, installation should be straightforward.

## Installation

1. Github (Development Version)


```{r install_developter, eval=FALSE}

devtools::install_github("wilsonlabgroup/DEET")

```


2. CRAN (Stable Release)

```{r install_cran, eval=FALSE}

install.packages("DEET")

```






