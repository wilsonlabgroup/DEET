#' @title DEET_enrich
#'
#' @description Core function of DEET where an input weighted human gene list
#' will be queried to DEETs library of studies.
#'
#' @param DEG_list Data frame or matrix of gene symbols with corresponding padj
#' and log2FC values (3 columns in total). Can also be a character vector of
#' gene symbols only. colnames of genes: c("gene_symbol", "padj", "coef")
#' The rownames of the dataframe are also the gene symbols.
#' @param ordered Boolean value specifying whether DEG_list is a character
#' vector of gene symbols that is ordered. Default value is FALSE.
#' @param background Character vector of human gene symbols showing all
#' possible genes. Default value is NULL.
#' @param example_gmt Boolean value specifying whether we're using all of 
#' data within the DEET database or whether we are using the example set 
#' of 100 comparisons.
#' 
#'
#' @return Named list where each element contains 6 objects. Each object will
#' contain the results (enrichment or correlation) and corresponding metadata.
#' \itemize{
#'   \item AP_INPUT_BP_output - Enriched BPs of input gene list.
#'   \item AP_INPUT_TF_output - Enriched TFs of input gene list.
#'   \item AP_DEET_DE_output  - Enrichment of input gene list on DEETs studies.
#'   \item AP_DEET_BP_output  - Enrichment of BPs of input gene list on DEETs
#'   BPs of studies.
#'   \item AP_DEET_TF_output  - Enrichment of TFs of input gene list on DEETs
#'   TFs of studies.
#'   \item DE_correlations    - Correlation values of input gene list to DEETs
#'   studies (both Pearson and Spearman).
#' }
#'
#' @author Dustin Sokolowski, Jedid Ahn
#'
#' @examples
#' '\dontrun{
#' 
#' 
#' GC_outrider_list <- readRDS("GC_outrider_outlier_gene_list.rds")
#' GC_outrider_DEET <- list()
#' background <- read.table("background_genes_202109.txt", header = F, as.is = T, sep = "\t")$V1
#' for(i in names(GC_outrider_list)){
#'   # Preprocessing of input gene list.
#'   DEG_list <- GC_outrider_list[[i]]
#'   DEG_list <- DEG_list[ , c("genename", "pValue", "zScore")]
#'   colnames(DEG_list) <- c("gene_symbol", "padj", "coef")
#'   DEG_list$padj[DEG_list$padj > 0.049] <- 0.049
#'   DEG_list <- DEG_list[!duplicated(DEG_list$gene_symbol), ]
#'   rownames(DEG_list) <- DEG_list$gene_symbol
#'
#'   # Core function.
#'   GC_outrider_DEET[[i]] <- DEET::DEET_enrich(DEG_list, background = background)
#'   }
#' }
#'
#' @references
#' Paczkowska M, Barenboim J, Sintupisut N, et al. Integrative pathway
#' enrichment analysis of multivariate omics data. Nat Commun. 2020;11(1):735.
#' doi:10.1038/s41467-019-13983-9
#'
#' @export
#' @importFrom ActivePathways read.GMT makeBackground ActivePathways
#' @importFrom pbapply pblapply
#' @importFrom stats cor.test p.adjust var
#'
DEET_enrich <- function(DEG_list, ordered = FALSE, background = NULL, example_gmt = FALSE){
  message(paste("Query start date and time:", Sys.time()))

  # Internal data loaded through sysdata.rda.
  # # 1) Start by loading the necessary files.
  # # a) DEET files.
   #load("DEET_metadata.rda")  # DEET_metadata: DEET's metadata of studies.
   #rownames(DEET_metadata) <- DEET_metadata$DEET.ID
   #load("DEET_DE_final.rda")  # DEET_DE: DEET's DE of studies.
  #
  # # b) ActivePathway (AP) files.
   #gmt_BP <- ActivePathways::read.GMT("Human_GO_AllPathways_with_GO_iea_June_01_2021_symbol.gmt")
   #gmt_TF <- ActivePathways::read.GMT("Human_TranscriptionFactors_MSigdb_June_01_2021_symbol.gmt")
   #DEET_gmt_DE <- ActivePathways::read.GMT("DEET_DE.gmt")
   #DEET_gmt_BP <- ActivePathways::read.GMT("DEET_BP.gmt")
   #DEET_gmt_TF <- ActivePathways::read.GMT("DEET_TF.gmt")

  # ============================================================================

  # Official code starts here.
  if(!(is.data.frame(DEG_list) | is.matrix(DEG_list) | is.character(DEG_list))){
    stop("Input DEG list must be of class data.frame, matrix, or character.")
  }

  # Case 1: DEG_list is a character vector.
  if (is.character(DEG_list)){
    message(
      paste("Input is a character vector of genes: Converting to dataframe",
            "for downstream analysis.")
      )

    DEG_processed <- data.frame(gene_symbol = DEG_list)

    if (ordered){
      message(
        paste("Input gene list is considered ORDERED: Spearman correlation",
              "is interpretable, but Pearson correlation is not.")
      )

      padj <- 0.049
      for(i in 2:nrow(DEG_processed)) {
        padj[i] <- padj[i-1] * 0.95
      }
      padj <- rev(padj)
      log2FC <- rev(seq(1, 1 + 0.1*(nrow(DEG_processed) - 1), 0.1))

      DEG_processed$padj <- padj
      DEG_processed$coef <- log2fc
      colnames(DEG_processed) <- c("gene_symbol", "padj", "coef")
    } else {
      message(
        paste("Input gene list is considered UNORDERED: Correlation analysis",
              "will not be run and pathway enrichment will be unordered.")
      )

      DEG_processed$padj <- 0.049
      DEG_processed$coef <- 1
      colnames(DEG_processed) <- c("gene_symbol", "padj", "coef")
    }
  }

  # Case 2: DEG_list is a data frame or matrix.
  if (is.data.frame(DEG_list) | is.matrix(DEG_list)){
    message("Assuming input DEG list is ordered.")

    if (!all( c("gene_symbol", "padj", "coef") %in% colnames(DEG_list)) ){
      stop(
        paste("Data frame or matrix needs columns named 'gene_symbol', 'padj',",
              "and 'coef' (case sensitive) in order. 'coef' in the context of",
              "differential expression represents the fold-change.")
      )
    } else{
      # No processing required.
      DEG_processed <- DEG_list
    }
  }

  # Also check background input.
  if (!is.null(background)){
    if (!is.character(background)){
      stop("Background must be a character vector.")
    }
  }
  if(example_gmt) {
    data("DEET_example_data")
    DEET_DE <- DEET_example_data$DEET_DE
    DEET_gmt_BP <- DEET_example_data$DEET_gmt_BP
    DEET_gmt_TF <- DEET_example_data$DEET_gmt_TF
    DEET_gmt_DE <- DEET_example_data$DEET_gmt_DE
    DEET_metadata <- DEET_example_data$DEET_metadata
    gmt_BP <- DEET_example_data$gmt_BP
    gmt_TF <- DEET_example_data$gmt_TF
    
    
  }

  # ============================================================================

  # Enrichment algorithms starts here.
  comp <- as.matrix(DEG_processed[ , "padj"])
  rownames(comp) <- toupper(rownames(DEG_processed))

  # 1) Find enriched BPs of input gene list.
  # If background variable is NULL, make background using AP.
  message("Running pathway enrichment of input gene list to find enriched BPs.")
  if(is.null(background)) {
    background <- ActivePathways::makeBackground(gmt_BP)
  }
  AP_INPUT_BP <- ActivePathways::ActivePathways(scores = comp,
                                                gmt = gmt_BP,
                                                background = background,
                                                geneset.filter = c(15, 2000),
                                                merge.method = "Brown",
                                                correction.method = "fdr",
                                                significant = 1,
                                                cutoff = 0.05)
  AP_INPUT_BP_output <- AP_INPUT_BP[AP_INPUT_BP$adjusted.p.val < 0.05, ]

  # ----------------------------------------------------------------------------

  # 2) Find enriched TFs of input gene list.
  message("Running motif enrichment of input gene list to find enriched TFs.")
  if (is.null(background)){
    background <- ActivePathways::makeBackground(gmt_TF)
  }
  AP_INPUT_TF <- ActivePathways::ActivePathways(scores = comp,
                                                gmt = gmt_TF,
                                                background = background,
                                                geneset.filter = c(15, 5000),
                                                merge.method = "Brown",
                                                correction.method = "fdr",
                                                significant = 1,
                                                cutoff = 0.05)
  AP_INPUT_TF_output <- AP_INPUT_TF[AP_INPUT_TF$adjusted.p.val < 0.05, ]

  # ----------------------------------------------------------------------------

  # 3) Gene set enrichment of input gene list with DEET studies.
  message("Running gene set enrichment with DEET studies.")

  if (is.null(background)){
    background <- ActivePathways::makeBackground(DEET_gmt_DE)
  }

  AP_DEET_DE <- ActivePathways::ActivePathways(scores = comp,
                                               gmt = DEET_gmt_DE,
                                               background = background,
                                               geneset.filter = c(15, 10000),
                                               merge.method = "Brown",
                                               correction.method = "fdr",
                                               significant = 1,
                                               cutoff = 0.05)
  AP_DEET_DE_sig <- AP_DEET_DE[AP_DEET_DE$adjusted.p.val < 0.05, ]

  # ----------------------------------------------------------------------------

  # 4) Find enriched BPs of input gene list on DEET’s BPs of studies.
  comp_bp <- as.matrix(AP_INPUT_BP$adjusted.p.val)
  rownames(comp_bp) <- AP_INPUT_BP$term.id

  if(min(comp_bp) >= 0.05) {
    AP_DEET_BP_sig <- "Internal pathway enrichment of input gene list did not
    discover biological pathways matching cutoff."
    warning(AP_DEET_BP_sig)
  } else {
    message(
      paste("Internal pathway enrichment of input gene list discovered",
            "biological pathways matching cutoff.")
      )

    AP_DEET_BP <- ActivePathways::ActivePathways(scores = comp_bp,
                                                 gmt = DEET_gmt_BP,
                                                 geneset.filter = c(15, 10000),
                                                 merge.method = "Brown",
                                                 correction.method = "fdr",
                                                 significant = 1,
                                                 cutoff = 0.05)
    AP_DEET_BP_sig <- AP_DEET_BP[ AP_DEET_BP$adjusted.p.val < 0.05, ]
  }

  # ----------------------------------------------------------------------------

  # 5) Find enriched TFs of input gene list on DEET’s TFs of studies.
  comp_tf <- as.matrix(AP_INPUT_TF$adjusted.p.val)
  rownames(comp_tf) <- AP_INPUT_TF$term.id

  if(min(comp_tf) >= 0.05) {
    AP_DEET_TF_sig <- "Internal motif enrichment of input gene list did not
    discover transcription factors matching cutoff."
    warning(AP_DEET_TF_sig)
  } else {
    message(
      paste("Internal motif enrichment of input gene list discovered",
            "transcription factors matching cutoff.")
      )

    AP_DEET_TF <- ActivePathways::ActivePathways(scores = comp_tf,
                                                 gmt = DEET_gmt_TF,
                                                 geneset.filter = c(15, 10000),
                                                 merge.method = "Brown",
                                                 correction.method = "fdr",
                                                 significant = 1,
                                                 cutoff = 0.05)
    AP_DEET_TF_sig <- AP_DEET_TF[ AP_DEET_TF$adjusted.p.val < 0.05, ]
  }

  # ----------------------------------------------------------------------------

  # 6) Perform a correlation test of input gene list to DEET’s studies
  # (both Pearson and Spearman).
  if (var(DEG_processed$coef) <= 0){
    cor_results_sig <- "No variance in coefs. Cannot proceed with correlation."
    warning(cor_results_sig)
  } else{
    cor_results_sig <- NULL # Initialize variable.
    AP_DEET_DE_sig_sub <- AP_DEET_DE_sig[lengths(AP_DEET_DE_sig$overlap) > 2, ]

    if (nrow(AP_DEET_DE_sig_sub) > 0){
      message(
        paste("There are significantly enriched studies with 3 or more",
              "overlapping DEGs. Measuring correlation in distribution of",
              "overlapping DEGs.")
      )

      DEET_DE_sub <- DEET_DE[AP_DEET_DE_sig_sub$term.id]
      DEET_metadata_sub <- DEET_metadata[AP_DEET_DE_sig_sub$term.id, ]

      tst <- pbapply::pblapply(1:length(DEET_DE_sub),
                               FUN = single_gene_set_cor_test,
                               DEG_processed = DEG_processed,
                               DEET_DE = DEET_DE_sub,
                               AP_DEET_DE_sig = AP_DEET_DE_sig_sub,
                               DEET_metadata = DEET_metadata_sub)
      names(tst) <- DEET_metadata_sub$DEET.ID

      # Retrieve correlation statistics and multiple-test correct.
      cor_results <- as.data.frame(
        do.call("rbind", lapply(tst, function(x) return (x$summary_cor)))
      )
      cor_results$FDR_Pear <- p.adjust(cor_results$P_Pear)
      cor_results$FDR_Spear <- p.adjust(cor_results$P_Spear)

      cor_results$Pear <- as.numeric(cor_results$Pear)
      cor_results$Spear <- as.numeric(cor_results$Spear)

      cor_results$P_Pear <- as.numeric(cor_results$P_Pear)
      cor_results$P_Spear <- as.numeric(cor_results$P_Spear)

      # Retrieve distributions of DE genes.
      cor_mats <- lapply(tst, function(x) return (x$matrix))

      # Make sure that correlations are significant with Pearson or Spearman.
      cor_results_sig <- cor_results[cor_results$FDR_Pear < 0.05 |
                                       cor_results$FDR_Spear < 0.05, ]

      cor_mats_sig <- cor_mats[cor_results$FDR_Pear < 0.05 |
                                 cor_results$FDR_Spear < 0.05]

      rownames(cor_results_sig) <- cor_results_sig$DEET.ID
      names(cor_mats_sig) <- cor_results_sig$DEET.ID
    }
  }

  # ----------------------------------------------------------------------------

  # Generate function output elements by matching metadata to results.

  # 3) AP_DEET_DE_output
  if (class(AP_DEET_DE_sig) == "data.table" && nrow(AP_DEET_DE_sig) > 0){
    meta_match <- DEET_metadata[AP_DEET_DE_sig$term.id, ]

    AP_DEET_DE_output <- list(results = AP_DEET_DE_sig,
                              metadata = meta_match)
  } else{
    AP_DEET_DE_output <- paste("Your gene list did not significantly enrich",
                               "any studies.")
    warning(AP_DEET_DE_output)
  }

  # 4) AP_DEET_BP_output
  if (class(AP_DEET_BP_sig) == "data.table" && nrow(AP_DEET_BP_sig) > 0){
    meta_match <- DEET_metadata[AP_DEET_BP_sig$term.id, ]

    AP_DEET_BP_output <- list(results = AP_DEET_BP_sig,
                              metadata = meta_match)
  } else{
    AP_DEET_BP_output <- paste("Your gene list converted to BP enrichment did",
                               "not significantly enrich any studies.")
    warning(AP_DEET_BP_output)
  }

  # 5) AP_DEET_TF_output
  if (class(AP_DEET_TF_sig) == "data.table" && nrow(AP_DEET_TF_sig) > 0){
    meta_match <- DEET_metadata[AP_DEET_TF_sig$term.id, ]

    AP_DEET_TF_output <- list(results = AP_DEET_TF_sig,
                              metadata = meta_match)
  } else{
    AP_DEET_TF_output <- paste("Your gene list converted to TF enrichment did",
                               "not significantly enrich any studies.")
    warning(AP_DEET_TF_output)
  }

  # 6) DE_correlations
  if(is.data.frame(cor_results_sig) && nrow(cor_results_sig) > 0) {
    meta_match <- DEET_metadata[rownames(cor_results_sig), ]

    # Output consists of the studies with significant correlations, the log2FC
    # of genes DE in one study, and the metadata for those correlations.
    DE_correlations <- list(results = cor_results_sig,
                            distributions = cor_mats_sig,
                            metadata = meta_match)
  } else{
    DE_correlations <- "No studies significantly correlate to your gene list."
    warning(DE_correlations)
  }

  # ----------------------------------------------------------------------------

  message(paste("Query end date and time:", Sys.time()))

  output <- list(
    AP_INPUT_BP_output = AP_INPUT_BP_output,
    AP_INPUT_TF_output = AP_INPUT_TF_output,
    AP_DEET_DE_output = AP_DEET_DE_output,
    AP_DEET_BP_output = AP_DEET_BP_output,
    AP_DEET_TF_output = AP_DEET_TF_output,
    DE_correlations = DE_correlations
    )

  return(output)
}


#' @title single_gene_set_cor_test
#'
#' @description Calculate the Pearson and Spearman correlation between the
#' input DEGs and an individual set of DEGs from DEET.
#'
#' @param index Value between 1 and length(DEET_DE).
#' @param DEG_processed Processed data frame consisting of the input genes.
#' The 3 columns are gene_symbol, padj, and coef.
#' @param DEET_DE List containing DEET's DE and library of studies. This
#' contains a subset of studies with 3 or more overlapping DEGs.
#' @param AP_DEET_DE_sig Data frame consisting of the enrichment of input gene
#' list on DEET’s studies. This contains a subset of studies with 3 or more
#' overlapping DEGs.
#' @param DEET_metadata Corresponding metadata of studies in DEET. This
#' contains a subset of studies with 3 or more overlapping DEGs.
#'
#' @return A list containing 2 objects.
#' \itemize{
#'   \item summary_cor - Summary of the correlations.
#'   \item matrix - Scatterplot of genes that are DE in the study.
#' }
#'
#' @author Dustin Sokolowski, Jedid Ahn
#'
single_gene_set_cor_test <- function(index, DEG_processed, DEET_DE,
                                     AP_DEET_DE_sig, DEET_metadata){
  # Get comparison in DEET.
  comp <- DEET_DE[[index]]
  DEET_comp_id <- DEET_metadata$DEET.ID[index]
  DEET_comp_name <- DEET_metadata$DEET.Name[index]

  # Get genes DE in either study.
  genes <- unique(c(rownames(DEG_processed), rownames(comp)))
  genes_cor <- AP_DEET_DE_sig$overlap[[index]]

  # Build matrix of fold-changes of genes in either study.
  mat <- as.data.frame(matrix(0, nrow = length(genes_cor), ncol = 2))
  colnames(mat) <- c("input", "DEET")
  rownames(mat) <- genes_cor

  # Populate input fold-change.
  mat$input <- DEG_processed[DEG_processed$gene_symbol %in% genes_cor, ]$coef

  # Populate DEET fold-change.
  mat$DEET <- comp[genes_cor, ]$log2FoldChange


  # Perform Pearson correlation.
  pearson_cor <-  cor.test(mat[, 1],mat[, 2], method = "pearson")
  Pear <- pearson_cor$estimate[[1]]
  P_Pear <- pearson_cor$p.value[[1]]

  # Perform Spearman correlation.
  spearman_cor <-  cor.test(mat[, 1],mat[, 2], method = "spearman")
  Spear <- spearman_cor$estimate[[1]]
  P_Spear <- spearman_cor$p.value[[1]]

  # Get overlapping DEGs in study.
  mat <- as.data.frame(matrix(0, nrow = length(genes), ncol = 2))
  colnames(mat) <- c("input", "DEET")
  rownames(mat) <- genes

  # Populate input fold-change
  mat[DEG_processed$gene_symbol, "input"] <- DEG_processed$coef
  # Populate DEET fold-change
  mat[rownames(comp), "DEET"] <- comp$log2FoldChange

  col <- rep("grey", nrow(mat))
  col[(mat[,1] > 0 & mat[,2] > 0) | (mat[,1] < 0 & mat[,2] < 0) ] <- "purple"
  col[(mat[,1] < 0 & mat[,2] > 0) | (mat[,1] > 0 & mat[,2] < 0) ] <- "orange"
  mat$color <- col

  sameDir <- length(mat$color == "purple")
  oppositeDir <- length(mat$color == "purple")
  genes <- rownames(mat)[mat$color != "grey"]

  summary_cor <- c(DEET_comp_id, DEET_comp_name, Pear, P_Pear, Spear, P_Spear)
  names(summary_cor) <- c("DEET.ID", "DEET.Name", "Pear", "P_Pear",
                          "Spear", "P_Spear")

  l <- list(summary_cor = summary_cor, matrix = mat)
  return (l)
}

# [END]
