#' @title DEET_Input_as_Reference
#'
#' @description Alternative function to DEET enrich for when the inputted gene
#' list is unordered. Here, we can increase the statistical rigour of enrichment
#' by levaraging the p-values of the DEGs within DEET. Specifically, the inputted
#' DE list is used as the reference and we test each DE list against your reference.
#' Specifically. We convert your reference into a gmt file before inputting each pairwise
#' DE list into ActivePathways. This function does not complete correlations or pathway-level
#' analysis.
#'
#' @param genes A character vector of gene symbols within `DEET_dataset`
#' @param DEET_dataset The databank of the differential expression enrichment tool.
#' Appropriate inputs here are "DEET_example_data" stored within DEET, the "DEET_combined.rda" file
#' from the DEET stable repositoy found at X, and the DEET database developmental repository found at Y.
#' The DEET_dataset is a named list where details of it's structure can be found ?DEET_example_data.
#' Unlike in DEET_enrich, this dataset does not require the pathway-relevant elements
#' of the DEET_dataset list, namely "gmt_BP", or "gmt_TF" "DEET_gmt_BP",  "DEET_gmt_TF".
#' It also does not need DEET_gmt_DE.
#' @param background Character vector of human gene symbols showing all
#' possible genes. Default value is NULL and the background is generated
#' as all detected DEGs across any comparison.
#'
#'
#' @return Named list containing the ActivePathways enrichment of each comparison
#' on the user's inputted gene list, as well as the associated metadata of 
#' each enriched comparison.
#'
#' @author Dustin Sokolowski
#'
#' @examples
#'
#' data("example_DEET_enrich_input")
#' genes <- rownames(example_DEET_enrich_input)
#' data("DEET_example_data")
#' DEET_out_ref <- DEET_Input_as_Reference(example_DEET_enrich_input, DEET_dataset = DEET_example_data)
#'
#'
#' @references
#' Paczkowska M, Barenboim J, Sintupisut N, et al. Integrative pathway
#' enrichment analysis of multivariate omics data. Nat Commun. 2020;11(1):735.
#' doi:10.1038/s41467-019-13983-9
#'
#' @export
#' @importFrom utils data
#' @importFrom ActivePathways makeBackground ActivePathways
#' @importFrom pbapply pblapply
#' @importFrom stats p.adjust
#'

DEET_Input_as_Reference <- function(genes, DEET_dataset, background = NULL) {
 
  if(!is.character(genes)) {
    stop("Gene list must be a character vector of gene symbols")
  }
  if(length(genes) < 3) {
    stop("Gene list must contain at least 3 genes.")
  }
  
  
  
  # Also check background input.
  if (!is.null(background)){
    if (!is.character(background)){
      stop("Background must be a character vector.")
    }
  }
  nameRight <- all(  c("DEET_DE", "DEET_metadata") %in% names(DEET_dataset) )
  
  if(!nameRight) {
    stop("DE Dataset is missing element 'DEET_DE' and/or 'DEET_metadata')")
  }
  
  DEET_DE <- DEET_dataset$DEET_DE

  DEET_metadata <- DEET_dataset$DEET_metadata
  rownames(DEET_metadata) <- DEET_metadata$DEET.ID
  
  DEET_metadata <- DEET_metadata[names(DEET_DE),]
  
  Ref <- list()
  Ref[[input_name]] <- list(id=input_name,
                            name=input_name,
                            genes=genes)
  Ref[["Dummy"]] <- list(id="Dummy",
                         name="Dummy",
                         genes="AaaAaaAa1", "AaaAaaAa2","AaaAaaAa3")
  
  class(Ref) <- "GMT"
  
  input_name = "input"
  
  if(is.null(background)) {
    message("Generating background from all genes detected as DE in any study.")
    
    background <- unique(unname(unlist(lapply(DEET_DE, rownames))))
  }  


  single_AP <- function(x) {
    
    comp <- as.matrix(x$padj)
    rownames(comp) <- rownames(x)
    colnames(comp) <- "comp"
    AP <- suppressMessages(ActivePathways::ActivePathways(comp,gmt = Ref,background = background,geneset.filter = c(1,1e5),significant = 1))
    return(AP)
  }
  message("Running ActivePathways of each gene list using input as a reference.")
  AP_out <- pbapply::pblapply(DEET_DE,single_AP)

  message("Postprocessing Outputs")
  
  enriched_studies <- !unname(unlist(lapply(AP_out,is.null)))
  
  if(!("TRUE" %in% as.character(unique(enriched_studies)))) {
    return("No gene lists eriched for your inputted gene set.")
  }

  
  AP_out1 <- AP_out[!unname(unlist(lapply(AP_out,is.null)))]
  ID_name <- DEET_metadata[!unname(unlist(lapply(AP_out,is.null))),]
  
  AP_out2 <- do.call("rbind",AP_out1)
  AP_out2 <- AP_out2[AP_out2$term.id == input_name & AP_out2$term.name == input_name,]
  AP_out2$term.id <- ID_name$DEET.ID
  AP_out2$term.name <- ID_name$DEET.Name
  AP_out2$p.val <- AP_out2$adjusted.p.val 
  
  AP_out2$adjusted.p.val <-stats::p.adjust(AP_out2$p.val, method = "fdr")
  AP_out2_sig <- AP_out2[AP_out2$adjusted.p.val < 0.05,]
  ID_name_sig <- ID_name[AP_out2$adjusted.p.val < 0.05,]
  AP_DEET_DE_output <- list(results=AP_out2_sig,
                            metadata=ID_name_sig)
  
  return(AP_DEET_DE_output)
}
