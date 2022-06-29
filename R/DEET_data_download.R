#' @title DEET_data_download
#'
#' @description Function to automatically download the files within
#' the DEET database that are required for the DEET_enrich
#'  and DEET_feature_extractfunctions.
#'
#' @param x categorical variable containing options "ALL", "enrich",
#' "metadata" or "feature_matrix".
#'
#'
#' @return Named list with the neccesary data required to input into
#' DEET_feature_extract or DEET_enrich. The metadata within DEET
#' can also be downloaded.
#'
#' \itemize{
#'   \item feature_matrix - A gene by comparison matrix populated
#'    with the log2FC of gene expression for all genes, regardless
#'    of DE status.
#'   \item metadata - a comparison - by - explanatory piece of data dataframe
#'    providing important details to contextualize each study.
#'    For every pairwise comparison, the study name,
#'    source (SRA, TCGA, GTEx and SRA-manual),
#'    description from the DRA compendium,
#'    the number of samples (total, up-condition,
#'    and down-condition), samples (total ,up-condition, down-condition),
#'    tissue (including tumour from TCGA),
#'    number of DEs (total, up-condition, down-condition),
#'    age (mean +- sd), sex, top 15 DEGs - up, top 15 DEGs - down,
#'    top 5 enriched pathways, and top 5 enriched TFs.
#'    PMID are also available for studies selected from SRA.
#'    Lastly, each pairwise comparison was given an overall
#'    category based on those decided in Crow et al., 2019.
#'   \item DEET_enrich - A named list of seven objects containing
#'    the data frames summarizing the DEGs from comparisons within
#'    DEET, GMT objects of comparisons within DEET for enrichment
#'    through ActivePathways, GMT objects for basic pathway and TF
#'    enrichment, and a dataframe for the metadata of each study.
#'    For more detail on each element of the list, please consult
#'    the vignette or "?DEET_example_data", as it is a subset
#'    of this object
#' }
#'
#' @author Dustin Sokolowski, Jedid Ahn
#'
#' @examples
#' \donttest{
#'
#' # Download the metadata. Downloading other
#' # files within DEET are larger and take
#' # a bit more time.
#' downloaded <- DEET_data_download(x = "metadata")
#'
#' # extract metadata from the list
#' metadata <- downloaded[["metadata"]]
#'
#'}
#'
#' @references
#' Engebretsen, S., & Bohlin, J. (2019). Statistical predictions with glmnet.
#' Clinical epigenetics, 11(1), 1-3.
#'
#' @export
#' @importFrom downloader download
#'
DEET_data_download <- function(x = "enrich") {
  #dummy_downloads
  DEET_metadata <- DEET_combined <- DEET_log2FC_matrix <- ""
  inputCheck <- c(x %in% c("ALL", "enrich", "metadata", "feature_matrix"))
  if(!inputCheck) {
  stop("Inputted variable was not in options: ALL, enrich, metadata, feature_matrix.
  Case sensitive. Please input one of these options or download manually from
  http://wilsonlab.org/public/DEET_data/.")
  }

  l <- list()

  if(x == "ALL") {
    x <- c("enrich", "metadata", "feature_matrix")
  }

  for(i in x) {
    if(i %in% "metadata") {
      datafile <- "DEET_metadata.rda"
      metafile <- paste0(datafile)
      url <- paste0("http://wilsonlab.org/public/DEET_data/",
                    metafile, "?raw=true")
      destfile <- file.path(tempdir(), metafile)
      downloader::download(url, destfile = destfile, mode = "wb")
      load(destfile)
      l[["metadata"]] <- DEET_metadata
    }
    if(i %in% "enrich") {
      datafile <- "DEET_combined.rda"
      metafile <- paste0(datafile)
      url <- paste0("http://wilsonlab.org/public/DEET_data/",
                    metafile, "?raw=true")
      destfile <- file.path(tempdir(), metafile)
      downloader::download(url, destfile = destfile, mode = "wb")
      load(destfile)
      l[["DEET_enrich"]] <- DEET_combined
    }
    if(i %in% "feature_matrix") {
      datafile <- "DEET_log2FC_matrix.rda"
      metafile <- paste0(datafile)
      url <- paste0("http://wilsonlab.org/public/DEET_data/",
                    metafile, "?raw=true")
      destfile <- file.path(tempdir(), metafile)
      downloader::download(url, destfile = destfile, mode = "wb")
      load(destfile)
      l[["DEET_feature_extract"]] <- DEET_log2FC_matrix
    }
  }
  return(l)
}

