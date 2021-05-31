#' Function to run RNAseq feature selection
#'
#' This function run a RNAseq feature selection methods and takes as input a gene expression matrix
#'
#' available methods are cv and none
#'
#' #'pipes from ( ... %>% )  : [data.matrix] , [run_norm] , [run_trans]
#'
#' pipes to ( %>% ... ) : [run_deconv]
#'
#' returns a matrix
#'
#' @param mix_matrix The gene expression matrix samples *x* genes.
#' @param method The feature selection method to be used
#'
#' @return This function return a matrix  samples *x* genes
#'
#' @export


run_featsel <- function(mix_matrix, method) {
  cv_sel <- function(mat, nmarker) {
    mm <- apply(mat,1,mean)
    vv <- apply(mat,1,var)
    cv <- sqrt(vv) / (mm + 1)
    cv[is.na(cv)] <- 0
    ix <- sort(cv, dec = TRUE, index = TRUE)$ix
    index.select <- rownames(mat)[ix[1:nmarker]]
    return(mat[index.select, ])
  }
  if (!{
    method %in% c("cv1000", "cv5000", "none",  "None")
  }) {
    print("Unknown method argument, please specify a method within the following list: cv, none")
  }
  if (method %in%  c("none", "None")) {
    sel.mat <- mix_matrix
  } else if (method == "cv1000") {
    nmarker <- 1000
    sel.mat <- cv_sel(mix_matrix, nmarker)
  } else if (method == "cv5000") {
    nmarker <- 5000
    sel.mat <- cv_sel(mix_matrix, nmarker)
  }

  return(sel.mat)
}