#' @title DEET_plot_correlation
#'
#' @description Take significant correlation outputs and generate
#' scatterplots of the genes DE in one or the other.
#'
#' @param correlation_input The "DE_correlations" element of the
#' output of the DEET_enrich function. This function only works if
#' there is at least one significantly correlated study.
#'
#' @return Named list of ggplot objects with the correlation between
#' the input study and the study within DEET
#'
#' @author Dustin Sokolowski, Jedid Ahn
#'
#' @examples
#' \donttest{
#'
#' data("example_DEET_enrich_input")
#' data("DEET_example_data")
#' DEET_out <- DEET_enrich(example_DEET_enrich_input, DEET_dataset = DEET_example_data)
#' correlation_input <- DEET_out$DE_correlations
#' correlation_plots <- DEET_plot_correlation(correlation_input)
#'
#'}
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_bw xlab ylab theme element_text
#' @importFrom ggrepel geom_text_repel
#'



DEET_plot_correlation <- function(correlation_input ) {


  # dummy variables for ggplot

  input <- DEET <- color <- label <- ""

if(!is.list(correlation_input)) {
  stop("correlation_input is not of class list.")
}
if(nrow(correlation_input$results) < 1) {
  stop("There were no significant correlations, ")
}

  # get DEET IDs
  distNames <- names(correlation_input$distributions)
  message(paste0("There are ", length(distNames), " comparisons with a significant correlation."))

  rownames(correlation_input$metadata) <- correlation_input$metadata$DEET.ID

  output_plots <- list()

  for(i in distNames) {
    metaColumn <- correlation_input$metadata[i,]

    plts <- distribution <- correlation_input$distributions[[i]]

    name <- rep("", nrow(plts))
    name[plts$color != "grey"] <- rownames(plts)[plts$color != "grey"]
    plts$label <- name
    plts$color <- as.factor(plts$color)
    ttl <- metaColumn[i,"DEET.Name"][[1]]
    cor_out <- correlation_input$results[i,]
    cor_out1 <- cor_out[i,"Spear"]
    corLab <- round(cor_out1,3)

    p <- ggplot2::ggplot(plts, ggplot2::aes(x=input, y=DEET, color=color, label=label)) + ggplot2::geom_point() +
      ggplot2::scale_color_manual(values=c("grey", "orange", "purple")) +
      ggrepel::geom_text_repel(size=2,max.overlaps = 30) + ggplot2::theme_bw() +
      ggplot2::xlab(paste0("Input coefficient Rho: ", corLab)) + ggplot2::ylab(paste0("DEET-",ttl)) + ggplot2::theme(axis.title.y = ggplot2::element_text(size = 8))

    output_plots[[i]] <- p

  }

  return(output_plots)

}
