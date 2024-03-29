% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/res_enrich.R
\name{enrich}
\alias{enrich}
\title{Function to perform biological enrichment analysis.}
\usage{
enrich(
  mydata,
  pathways = NULL,
  showCategory = 10,
  showLeadingGenes = FALSE,
  fdr = TRUE,
  multilevel = FALSE,
  ICAbased = FALSE
)
}
\arguments{
\item{mydata}{T matrix.}

\item{pathways}{List of gene sets to check}

\item{showCategory}{Maximum of enriched terms to display}

\item{showLeadingGenes}{Display or not leading genes for each significant term}

\item{fdr}{Whether correction for multiple testing should be considered}

\item{multilevel}{Whether fgseamult should be used}

\item{ICAbased}{If ICA based deconvolution algorithm was used to obtain the T matrix.}
}
\value{
return a list with :
enrichment results for each column (cell type) of the T matrix.
}
\description{
This function perform biological enrichment analysis for the different component of the T matrix
}
