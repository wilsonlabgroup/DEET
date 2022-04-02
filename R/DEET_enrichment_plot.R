#' @title DEET_enrichment_plot
#'
#' @description Generate barplots or dotplots from the output of DEET
#'
#' @param enrich_list A list of enrichments from DEET, with each element
#'  post-processed with the barplot enrichment function.
#' @param outname A character giving the title of the barplot or dotplot.
#' @param width The number of inches in the barplot or dotplot.
#' @param text_angle The angle of the enriched studies.
#' @param horizontal Whether the output barplot is vertical or horizontal
#' @param topn the top number of studies (by p-value) to be plotted.
#' @param ol_size the minimum number of overlapping genes (or paths) in an enriched study.
#' @param exclude_domain Exclude studies enriched based on DEGs, Paths, or TF if
#'  the user happened to aggregate the results into a single DF, generally unused.
#' @param cluster_order Factor to group studies based on the researchers custom annotation.
#' @param dot logical (T/F) of whether to produce a dotplot or a barplot
#' @param colors Type of color pallete to input into 'scale_fill_brewer' of ggplot.
#' @param split_domain logical (T/F) of whether to plot the "topn" studies for
#'  each "domain" (default is source) or to plot the topn pathwys regardless of domain.
#'  default is set to FALSE, meaning it plots the topn pathways regardless of domain.
#'
#'
#' @return A ggplot2 object (barplot or dotplot) of enrichment identified within DEET.
#'

#' @author Dustin Sokolowski, Hauyun Hou PhD
#'
#' @examples
#'
#' data("example_DEET_enrich_input")
#' data("DEET_example_data")
#' DEET_out <- DEET_enrich(example_DEET_enrich_input, DEET_dataset = DEET_example_data)
#'
#' # converting output to format compatible with DEET_enrichment plot
#' DE_example <- DEET_out$AP_DEET_DE_output$results
#' DE_example$term.name <- DEET_out$AP_DEET_DE_output$metadata$DEET.Name
#' DE_example$domain <- "DE"
#' DE_example$overlap.size <- lengths(DE_example$overlap)
#' DE_example$p.value <- DE_example$adjusted.p.val
#'
#' DE_example_plot <- DEET_enrichment_plot(list(DE_example = DE_example), "DE_example")
#'
#'
#' @export
#' @importFrom utils data
#' @importFrom dplyr arrange desc mutate row_number %>%
#' @importFrom ggplot2 ggplot geom_bar aes geom_text coord_flip theme_bw theme element_text scale_x_continuous scale_fill_brewer facet_grid theme geom_point
#' @importFrom stats na.omit
#'


DEET_enrichment_plot <- function(enrich_list, outname, width=8, text_angle=0, horizontal =FALSE, topn=5, ol_size=1, exclude_domain="", cluster_order=NULL, dot=FALSE, colors = "Set2", split_domain = FALSE){

  #empty variables to deal with column-name note fo dplyr
  n <- cluster <- domain <- p.value <- overlap.size <- term.id <- term.name <- ""


  # Internal functions used
  count_intersection <- function(x){
    return(length(unlist(strsplit(x, ","))))
  }

  # Get plot data ready for ggplot
  process_plotdata <- function(plotdata, exclude_domain="", testclusters=""){
    plotdata <- plotdata[order(plotdata$domain, plotdata$p.value),]
    plotdata$term.name <- factor(plotdata$term.name, levels=unique(as.character(plotdata$term.name)))
    #plotdata$domain <- factor(plotdata$domain, levels=c("BP", "MF", "CC", "keg", "hp"))
    plotdata <- subset(plotdata, !domain %in% exclude_domain)
    plotdata$n <- paste0("ol=", plotdata$`overlap.size`)
    plotdata$cluster <- factor(plotdata$cluster, levels=testclusters)

    # This is just so that the "R code for possible problems" check in
    # devtools is passed, as it has a hard time knowing col names
    # in tidyr are not internal variables



    plotdata <- plotdata %>%
      dplyr::arrange(cluster, dplyr::desc(domain), -p.value) %>%
      dplyr::mutate(order = dplyr::row_number())
    return(plotdata)
  }

  testclusters <- names(enrich_list[sapply(enrich_list, nrow) >0])
  # organize data frame for enrichment group
  enrich_data_list <- lapply(testclusters, function(x){
    print(x)
    enrich_data <- enrich_list[[x]]
    enrich_data$domain <- factor(enrich_data$domain, levels = as.character(unique(enrich_data$domain)))
    enrich_data <- do.call("rbind", lapply(split(enrich_data, f=enrich_data$domain), function(x) {
      x <- subset(x, `overlap.size` >= ol_size)
      x <- stats::na.omit(x[order(x$p.value),][1:min(nrow(x), topn),])
      return(x)
    }))

    enrich_data$term.name <- factor(enrich_data$term.name, levels=rev(unique(as.character(enrich_data$term.name))))
    enrich_data$cluster <- rep(x, nrow(enrich_data))
    return(enrich_data)
  })
  # Get data for plotting together into a df
  plotdata <- do.call("rbind", enrich_data_list)
  plotdata <- process_plotdata(plotdata, exclude_domain=exclude_domain, testclusters = testclusters)

  if(nrow(plotdata) > topn) {
  if(!split_domain) {

    plotdata <- plotdata[order(plotdata$p.value),]
    plotdata <- plotdata[1:topn,]
    #plotdata$order <- plotdata$order - min(plotdata$order)
    plotdata$order <- order(plotdata$order)

  }
  }

  plotdata_all <- do.call("rbind", lapply(testclusters, function(x) {
    enrich_list[[x]]$cluster <- x
    return(enrich_list[[x]])}))
  plotdata_all <- process_plotdata(plotdata_all,exclude_domain=exclude_domain, testclusters = testclusters)
  if(!is.null(cluster_order)){
    plotdata$cluster <- factor(plotdata$cluster, levels=cluster_order)
    plotdata_all$cluster <- factor(plotdata_all$cluster, levels=cluster_order)
  }
  #pdf(paste0(outname, "_gprofiler_enrichment_top", topn, "_", Sys.Date(), ".pdf"), width=width, height = nrow(plotdata)*0.15+1)
  textpos <- 0.80*(-log10(min(plotdata$p.value)))
  # make the barplot
  p <- ggplot2::ggplot(plotdata) +
    #geom_bar(aes(x=`term.name`, y=-log10(p.value), fill=domain), stat="identity") +
    ggplot2::geom_bar(ggplot2::aes(x=order, y=-log10(p.value), fill=domain), stat="identity") +
    ggplot2::geom_text(ggplot2::aes(x=order, y=0.7*(-log10(p.value)), label=n)) +
    ggplot2::coord_flip()+
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text = ggplot2::element_text(color="black", size=10)) +
    ggplot2::scale_x_continuous(
      breaks = plotdata$order,
      labels = plotdata$`term.name`,
      expand = c(0,0)) +
    ggplot2::scale_fill_brewer(palette = colors)
  p <- p +
    ggplot2::facet_grid(cluster~., scales = "free", space = "free") +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle=text_angle))
  if(horizontal & !(dot)){
    used_data <- subset(plotdata_all, `term.id` %in% plotdata$term.id)
    used_data <-  dplyr::arrange(used_data, dplyr::desc(domain),-p.value, cluster)
    idorder <- c(1:length(unique(used_data$term.name)))
    names(idorder) <- unique(used_data$term.name)
    used_data$order <- idorder[used_data$term.name]
    used_data$term.name <- factor(used_data$term.name, levels = unique(used_data$term.name))
    p <-  ggplot2::ggplot(used_data) +
      #geom_bar(aes(x=`term.name`, y=-log10(p.value), fill=domain), stat="identity") +
      ggplot2::geom_bar( ggplot2::aes(x=`term.name`, y=-log10(p.value), fill=domain), stat="identity") +
      ggplot2::geom_text( ggplot2::aes(x=`term.name`, y=0.7*(-log10(p.value)), label=n)) +
      ggplot2::coord_flip()+
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text =  ggplot2::element_text(color="black", size=10),
            axis.text.x =  ggplot2::element_text(angle=text_angle)) +
      ggplot2::facet_grid(~cluster, scales = "free_x") +
      ggplot2::scale_fill_brewer(palette = colors)
    }
    if(dot){
      # make the dotplot
      used_data <- subset(plotdata_all, `term.id` %in% plotdata$term.id)
      used_data <-  dplyr::arrange(used_data, dplyr::desc(domain),-p.value, cluster)
      idorder <- c(1:length(unique(used_data$term.name)))
      names(idorder) <- unique(used_data$term.name)
      used_data$order <- idorder[used_data$term.name]
      used_data$term.name <- factor(used_data$term.name, levels = unique(used_data$term.name))
      if(horizontal) {
      p <-  ggplot2::ggplot(used_data) +
        ggplot2::geom_point( ggplot2::aes(x=`term.name`, y=cluster, size=-log10(p.value), color=domain)) +
        ggplot2::scale_color_brewer(palette = colors) +
        ggplot2::theme_bw() +
        ggplot2::coord_flip() +
        ggplot2::theme(axis.text =  ggplot2::element_text(color="black", size=10),
              axis.text.x =  ggplot2::element_text(angle=text_angle))
      } else {
        p <-  ggplot2::ggplot(used_data) +
          ggplot2::geom_point( ggplot2::aes(x=`term.name`, y=cluster, size=-log10(p.value), color=domain)) +
          ggplot2::scale_color_brewer(palette = colors) +
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text =  ggplot2::element_text(color="black", size=10),
                         axis.text.x =  ggplot2::element_text(angle=text_angle + 90))
      }


    }

  #print(p)
  #dev.off()
  return(p)
}

