#'Function to run RNAseq feature selection
#'
#'This function run a RNAseq feature selection methods and takes as input a gene expression matrix
#'
#'available methods are cv and none
#'
#'#'pipes from ( ... %>% )  : [data.matrix] , [run_norm] , [run_trans]   
#'
#'pipes to ( %>% ... ) : [run_deconv]
#'
#'returns a matrix
#'
#'@param mix_matrix The gene expression matrix samples *x* genes.
#'@param method The feature selection method to be used 
#'
#'@return This function return a matrix  samples *x* genes
#'
#'@importFrom PREDE  "select_feature"
#'
#'@export


run_featsel <- function(mix_matrix, method) {
  if ( !{ method %in% c("cv1000","cv5000",  "none") } ) {
    print("Unknown method argument, please specify a method within the following list: cv, none")
  }
  if (method == "none") sel.mat= mix_matrix
  else if (method == "cv1000") {
    feat = PREDE::select_feature(mat = mix_matrix ,method = "cv",nmarker = 1000,startn = 0)
    sel.mat = mix_matrix[feat,]
  } else if (method == "cv5000") {
    feat = PREDE::select_feature(mat = mix_matrix ,method = "cv",nmarker = 5000,startn = 0)
    sel.mat = mix_matrix[feat,]
  }
  
  return(sel.mat)
}