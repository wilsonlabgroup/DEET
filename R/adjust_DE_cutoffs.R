#' @title adjust_DE_cutoffs
#'
#' @description Utility function to adjust mean expression, FDR, and log2FC cutoffs of the database of DEGs inputted into DEET.
#'
#' @param DEET_combined The databank of the differential expression enrichment tool.
#' Appropriate inputs here are "DEET_example_data" stored within DEET, the "DEET_combined.rda" file
#' from the DEET stable repositoy found at X, and the DEET database developmental repository found at Y.
#' The DEET_dataset is a named list where details of it's structure can be found ?DEET_example_data.
#' @param redo_pathways_instructions Boolean value specifying whether to print the instructions
#' required to update all pathway enrichments based on new DE cutoffs.
#' @param baseMean Change the mean-expression cutoff.
#' @param abslog2FoldChange Change the log2 Fold-change cutoff. 
#' @param padj Change the FDR-adjusted p-value cutoff.
#'
#' @return The DEET_combined object but with the user-inputted 
#' expression, log2FC, and FDR-adjusted p-value cutoffs. DEET_gmt_DE
#' is also updated to the new cutoffs. 
#' 
#' 
#'
#' @author Dustin Sokolowski
#'
#' @examples
#'
#' 
#' data("DEET_example_data")
#' DEET_cutoff <- adjust_DE_cutoffs(DEET_example_data, abslog2FoldChange = 1, padj = 0.01)
#'
#'
#'
#' @export
#' @importFrom utils data
#'

adjust_DE_cutoffs <- function(DEET_combined, redo_pathways_instructions = FALSE, baseMean = 1, abslog2FoldChange = 0, padj = 0.05) {
  
  if(redo_pathways_instructions) {
    
    print("To recompute all of the pathwy enrichments for yourself, download >3Gb database of DEs from http://wilsonlab.org/public/DEET_data/intermediary_files/DEET_short_combined.rda,
    Then recompute pathways using ActivePathways (or any pathway enrichment tool of your choice).
    Our parameters are in https://figshare.com/articles/software/Code_and_data_describing_how_to_reproduce_the_Differential_Expression_Enrichment_Tool_DEET_consortium_/20425464/.
    Go to 4_postprocessing_scripts.R 5_DEET_AP.R Once this is finished, overlap the studies from the output with the final set of DEs as a few studies were flagged in the final metadata processing step.
    Since you'll be running over 6,000 pathway enrichments, this will take a number of hours on a cluster and we do not reccomend running it locally.")
    
  }
  
  message("This function only influences the Gene Level Analysis as  it is too computationally burdensome to re-compute over 6,000 pathway enrichments. 
          To get instructions on how to redo_pathway analysis, set redo_pathway instructions to 'TRUE'")
  
  if(all(baseMean == 1 & abslog2FoldChange == 0 & padj == 0.05)) {
    warning("You did not change any cutoffs from the default parameters in DEET, returning input.")
    return(DEET_combined)
  }
  if(!(all(c("DEET_DE", "DEET_gmt_DE") %in% names(DEET_combined) ))) {
    message("Check command returned false: all(c('DEET_DE', 'DEET_gmt_DE') %in% names(DEET_combined) ")
    stop("Missing key parts changed in DEET database by this function: DEET_DE, DEET_gmt_DE, or both. Check input database.")
  }
  if(padj > 0.05) {
    warning("Smallest available p-value cutoff is 0.05.
    A higher p-value threshold would lead to too many false positives. 
            Our Zendodo package has all p-values and fold-changes." )
  }
  
  DEET_DE <- DEET_combined$DEET_DE
  DEET_gmt <- DEET_combined$DEET_gmt_DE
  #DEET_DE1 <- list()
  for(i in names(DEET_DE)) {
    DE1 <- DEET_DE[[i]]
    DE1 <- DE1[(DE1[,1] > baseMean & abs(DE1[,2]) > abslog2FoldChange & DE1[,3] < padj),]
    DEET_DE[[i]] <- DE1
    
    DEET_gmt[[i]]$genes <- rownames(DE1)
    
  }
  
  lens <- unlist(lapply(DEET_DE, nrow))
  studies <- table(lens > 15)[["TRUE"]]
  percent <- round(studies/3162*100,2)
  message(paste0("With your new cut-offs, you now have ", studies, " comparisons with more than 15 DEGs from the original 3162 (", percent, "% comparisons)."))
  DEET_combined$DEET_DE <- DEET_DE
  DEET_combined$DEET_gmt_DE <- DEET_gmt
  return(DEET_combined)
}
