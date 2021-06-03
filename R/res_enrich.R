# function to reorient the components if ICA-based deconvolution algorithm was performed
.orient_funct <- function(S) {
  orient <-
    apply(S, 2, function(x) {
      if (min(x) < -3 & max(x) > 3) {
        ifelse(sum(x > 3) < sum(x < -3), -1, 1)
      } else {
        ifelse(sum(x > 2) < sum(x < -2), -1, 1)
      }
    })
  S <- as.matrix(S) %*% diag(orient)
  return(S)
}


# function to run the fgsea function and return a dataframe of the results
.enrichfun <- function(pathways,
                       genescores,
                       showCategory = 10,
                       showLeadingGenes = FALSE,
                       fdr = TRUE,
                       multilevel = FALSE) {
  genes <- sort(genescores, decreasing = TRUE)
  if (multilevel == FALSE) {
    if ("nperm" %in% names(formals(fgsea::fgsea))) {
      fgseaRes <- data.frame(fgsea(
        pathways,
        genes,
        minSize = 2,
        maxSize = 200,
        nperm = 2e4
      ))
    } else {
      fgseaRes <- data.frame(fgsea(pathways, genes,
        minSize = 2, maxSize =
          200
      ))
    }
  } else if (multilevel == TRUE) {
    if ("eps" %in% names(formal(fgsea::fgseaMultilevel()))) {
      fgseaRes <- data.frame(fgseaMultilevel(
        pathways,
        genes,
        minSize = 2,
        maxSize = 200,
        eps = 0
      ))
    } else {
      fgseaRes <- data.frame(fgseaMultilevel(pathways, genes,
        minSize = 2, maxSize =
          200
      ))
    }
  }
  pv <- ifelse(fdr == TRUE, "padj", "pval")
  res <- fgseaRes[, c("pathway", "padj", "pval", "ES")]
  if (showLeadingGenes == TRUE) {
    res$leadingEdge <- unlist(lapply(fgseaRes$leadingEdge, paste, collapse = "_"))
  }
  res <- res[which(res[, pv] < 0.05), ]
  res <- res[order(res[, pv], decreasing = FALSE), ]
  respos <- res[res[, "ES"] > 0, ]
  if (nrow(respos) == 0) {
    p <- print("No significant cell-type enrichment")
  } else {
    respos <- respos[1:min(showCategory, nrow(respos)), ]
    resf <- respos
  }
}

#' Function to perform biological enrichment analysis.
#'
#' This function perform biological enrichment analysis for the different component of the T matrix
#'
#' @param mydata T matrix.
#' @param pathways List of gene sets to check
#' @param showCategory Maximum of enriched terms to display
#' @param showLeadingGenes Display or not leading genes for each significant term
#' @param fdr Whether correction for multiple testing should be considered
#' @param multilevel Whether fgseamult should be used
#' @param ICAbased If ICA based deconvolution algorithm was used to obtain the T matrix.
#'
#' @return return a list with :
#' enrichment results for each column (cell type) of the T matrix.
#'
#'
#' @export
#'
#'
#'
enrich <- function(mydata,
                   pathways = NULL,
                   showCategory = 10,
                   showLeadingGenes = FALSE,
                   fdr = TRUE,
                   multilevel = FALSE,
                   ICAbased = FALSE) {
  if (ICAbased == TRUE) {
    mydata <- .orient_funct(mydata)
  }
  apply(
    mydata, 2,
    function(x) {
      .enrichfun(
        pathways = pathways,
        genescores = x,
        showCategory = showCategory,
        showLeadingGenes = showLeadingGenes,
        fdr = fdr,
        multilevel = multilevel
      )
    }
  )
}
