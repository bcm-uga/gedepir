#' Function to run RNAseq count transformation
#'
#' This function run a RNAseq count transformation methods and takes as input a raw gene expression matrix
#'
#' available methods are log2 pseudoLog (pseudolog)
#'
#' log2 : f(x)=log2(x+1)
#'
#' pseudoLog : f(x)=asinh(x)
#'
#' linear: f(x) = x
#'
#' pipes from ( ... %>% )  : [data.matrix] , [run_norm]
#'
#' pipes to ( %>% ... ) : [run_deconv]
#'
#' returns a matrix
#'
#' @param mix_matrix The gene expression matrix samples *x* genes.
#' @param method The transformation method to be used
#'
#' @return This function return a matrix  samples *x* genes
#'
#'
#' @export


run_trans <- function(mix_matrix, method) {
  if (!{
    method %in% c("log2", "pseudoLog", "pseudolog", "linear")
  }) {
    print("Unknown method argument, please specify a method within the following list: log2, pseudoLog, pseudolog")
  }
  if (method == "linear") {
    trans.mat <- mix_matrix
  } else if (method == "log2") {
    trans.mat <- log2(1 + mix_matrix)
  }
  else if (method %in% c("pseudolog", "pseudoLog")) {
    trans.mat <- asinh(mix_matrix)
  }

  return(trans.mat)
}