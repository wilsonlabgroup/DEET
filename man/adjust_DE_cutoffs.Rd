% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjust_DE_cutoffs.R
\name{adjust_DE_cutoffs}
\alias{adjust_DE_cutoffs}
\title{adjust_DE_cutoffs}
\usage{
adjust_DE_cutoffs(
  DEET_combined,
  redo_pathways_instructions = FALSE,
  baseMean = 1,
  abslog2FoldChange = 0,
  padj = 0.05
)
}
\arguments{
\item{DEET_combined}{The databank of the differential expression enrichment tool.
Appropriate inputs here are "DEET_example_data" stored within DEET, the "DEET_combined.rda" file
from the DEET stable repositoy found at X, and the DEET database developmental repository found at Y.
The DEET_dataset is a named list where details of it's structure can be found ?DEET_example_data.}

\item{redo_pathways_instructions}{Boolean value specifying whether to print the instructions
required to update all pathway enrichments based on new DE cutoffs.}

\item{baseMean}{Change the mean-expression cutoff.}

\item{abslog2FoldChange}{Change the log2 Fold-change cutoff.}

\item{padj}{Change the FDR-adjusted p-value cutoff.}
}
\value{
The DEET_combined object but with the user-inputted 
expression, log2FC, and FDR-adjusted p-value cutoffs. DEET_gmt_DE
is also updated to the new cutoffs.
}
\description{
Utility function to adjust mean expression, FDR, and log2FC cutoffs of the database of DEGs inputted into DEET.
}
\examples{


data("DEET_example_data")
DEET_cutoff <- adjust_DE_cutoffs(DEET_example_data, abslog2FoldChange = 1, padj = 0.01)



}
\author{
Dustin Sokolowski
}
