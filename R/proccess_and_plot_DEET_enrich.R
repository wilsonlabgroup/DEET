#' @title proccess_and_plot_DEET_enrich
#'
#' @description Generates barplots and dotplots based on the output of
#' the DEET_enrich function.
#'
#' @param DEET_output Direct output of the DEET_enrich function.
#' A list with all of the same names as DEET_output.
#' @param colour_barplot Pick dotplot or barplot colours.
#' It can be NULL, in which all bars are the same or it can
#' be a (case sensitive) column within the metadata. Defaults
#' to "source".
#' @param width The number of inches in the barplot or dotplot.
#' @param text_angle The angle of the enriched studies.
#' @param horizontal Whether the output barplot is vertical or horizontal
#' @param topn the top number of studies (by p-value) to be plotted.
#' @param ol_size the minimum number of overlapping genes (or paths) in an enriched study.
#' @param exclude_domain Exclude studies enriched based on DEGs, Paths, or TF if
#'  the user happened to aggregate the results into a single DF, generally unused.
#' @param cluster_order Factor to group studies based on the researchers custom annotation.
#' @param colors Type of color pallete to input into 'scale_fill_brewer' of ggplot.
#'
#' @return Named list where each element is a ggplot object plotting the output of
#' the enrichment tests within DEET. The final element is the output of ActivePathways
#' (in DEET) that is directly compatible with the DEET_enrichment_barplot function.
#'
#' \itemize{
#'   \item DEET_DotPlot - ggplot object of Dotplot of enrichment of enrichment of DEET studies
#'   based on DE, BP, and TF information. Only plotted if 2/3 levels contain
#'   at least one significant study.
#'   \item Pathway_barplot - ggplot object of Barplot of standard gene set enrichment based on
#'   gene ontology and TFs. Only plotted if there is at least one enriched
#'   significant pathway/TF.
#'   \item individual_barplot - ggplot object of Barplot of the top enriched pathways or studies
#'   (depending on the input list).Barplot is only generated if each list has at least one pathway (or study)
#'   is enriched.
#'   \item DEET_output_forplotting - output of Activepathways with "domain", "overlap.size", and "p.value"
#'   columns added to be compatible with the DEET_enrichment_barplot function.
#'
#' }
#'
#'
#' @author Dustin Sokolowski, Hauyun Hou PhD
#'
#' @examples
#'
#' data("example_DEET_enrich_input")
#' data("DEET_example_data")
#' DEET_out <- DEET_enrich(example_DEET_enrich_input, DEET_dataset = DEET_example_data)
#' plotting_example <- proccess_and_plot_DEET_enrich(DEET_out, text_angle = 45,
#' horizontal = TRUE, topn=4)
#'
#'
#' @export
#' @importFrom utils data
#'
proccess_and_plot_DEET_enrich <- function(DEET_output,colour_barplot = "Source", width=8, text_angle=0, horizontal =F, topn=5, ol_size=1, exclude_domain="", cluster_order=NULL, colors = "Set2") {

  message("Removing DE_correlations element from output")
  DEET_output <- DEET_output[names(DEET_output) != "DE_correlations"]

  message("Removing non-significant DE lists")

  DEET_output <- DEET_output[!("character" == unname(unlist(lapply(DEET_output,function(x) {
    return(class(x)[1])
  }))))]

  DEET_plot_processed <- list()
  for(i in names(DEET_output)) {
    DEET_list <- DEET_output[[i]]
    if(is.character(DEET_list)) {
      message(paste0("No significant results for ",i, " were detected and therefore cannot be plotted"))

    } else {
      if(grepl("AP_INPUT_",i)) {
        if((nrow(DEET_list) == 0)[1] ) {
          message(paste0("No significant results for ",i, " were detected and therefore cannot be plotted"))
          next
        }
      }

      domain <- gsub("AP_","",i)
      domain <- gsub("_output","",i)
      domain1 <- domain
      if(grepl("AP_INPUT",i)) {
        DEET_list$domain <- domain
        DEET_list$overlap.size <- lengths(DEET_list$overlap)
        DEET_list$p.value <- DEET_list$adjusted.p.val
        DEET_plot_processed[[domain]] <- DEET_list
      } else {
        DEET_out <- DEET_output[[i]]
        DEET_list <- DEET_out$results
        DEET_metadata <- DEET_out$metadata

        if(!is.null(colour_barplot)) {

        outMeta <- !(colour_barplot %in% colnames(DEET_metadata))[1]
        if(outMeta) {
          stop("'colour_barplot' variable is not in DEET's metadata (case sensitive). If you need custom colours then use 'DEET_enrichment_plot()' with a custom 'domain' value. ")
        }

        domain <- DEET_metadata[,colour_barplot]

        }

        DEET_list$term.name <- DEET_out$metadata$DEET.Name
        DEET_list$domain <- domain
        DEET_list$overlap.size <- lengths(DEET_list$overlap)
        DEET_list$p.value <- DEET_list$adjusted.p.val
        DEET_plot_processed[[domain1]] <- DEET_list
      }
    }
  }


message("Genarating batplot of traditional pathway enrichments")
each_barplot <- list()
for(i in names(DEET_plot_processed)) {
  each_barplot[[i]] <- DEET_enrichment_plot(DEET_plot_processed[i], i, width=width, text_angle=text_angle, horizontal =horizontal, topn=topn, ol_size = ol_size, exclude_domain = exclude_domain, cluster_order=cluster_order, dot=FALSE, colors = colors)
}

message("Genarating barplot of pathway enrichment (BP + TF) if available")
if(length(grep("INPUT", names(DEET_plot_processed))) == 2) {
  Traditional_Barplot <- DEET_enrichment_plot(DEET_plot_processed[grep("INPUT", names(DEET_plot_processed))], "Pathway_Enrichment", width=width, text_angle=text_angle, horizontal =horizontal, topn=topn, ol_size = ol_size, exclude_domain = exclude_domain, cluster_order=cluster_order, dot=FALSE, colors = colors)
} else {
  Traditional_Barplot <- "There was not significant enrichment of at least one gene set from both BPs or TFs, please see individual_barplot for results."
}

message("Genarating barplot of DEET enrichment (DE + BP + TF) if available")
if(length(grep("_DEET_", names(DEET_plot_processed))) > 1) {
  tmp <- DEET_plot_processed[grep("_DEET_", names(DEET_plot_processed))]
  for(i in names(tmp)) {
    tmp1 <- tmp[[i]]
    #tmp1$domain <- "DEET"
    tmp[[i]] <- tmp1
  }
  DEET_DotPlot <- DEET_enrichment_plot(tmp, "DEET_Enrichment", width=width, text_angle=text_angle, horizontal =TRUE, topn=topn, ol_size = ol_size, exclude_domain = exclude_domain, cluster_order=cluster_order, dot=TRUE, colors = colors)
} else {
  DEET_DotPlot <- "There was only enrichment of the DEET database from only on level (DE, BP, or TF). See result in 'individual_barplot`."
}

return(list(DEET_DotPlot = DEET_DotPlot, Pathway_barplot = Traditional_Barplot, individual_barplot = each_barplot, DEET_output_forplotting = DEET_plot_processed))

}

