#'Function to run deconvolution methods.
#'
#'This function run a deconvolution  methods and takes as input a pre-treated gene expression matrix
#'
#'@param mix_matrix The convoluted pre-treated matrix patients*probes.
#'@param k The number of cell types to estimate
#'@param method The deconvolution method to be used 
#'
#'@return This function return a list with the A and the T matrices computed.
#'
#'@importFrom NMF "nmf"
#'@importFrom fastICA "fastICA"
#'@importFrom PREDE "PREDE"
#'@importFrom deconica "generate_markers"
#'@importFrom deconica "get_scores"
#'
#'@export
#'
run_deconv <- function(mix_matrix, k =5, method) {

  
  if ( !{ method %in% c("NMF", "ICA", "CDSeq", "PREDE", "DSA") } ) {
    print("Unknown method argument, please specify a method within the following list: NMF, ICA, CDSeq, PREDE, DSA")
  }
  
  if (method == "NMF") {
    res <- NMF::nmf(x = mix_matrix, rank = k, method = "snmf/r", seed = 1)
    A   <- apply(
      X      = res@fit@H
      , MARGIN = 2
      , FUN    = function( x ) {
        x / sum( x )
      }
    )
    
    A_matrix <- A
    T_matrix <- res@fit@W
    remove(list = "res")
  }
  
  if (method == "ICA") {
    
    rna_data = mix_matrix[!duplicated(mix_matrix), ]
    
    ICA_deconv = fastICA::fastICA(X = rna_data, n.comp = k, maxit = 1000, tol = 1e-09)
    ICA_deconv$names = row.names(rna_data)
    weighted.list <- deconica::generate_markers(df = ICA_deconv,
                                                n = 30,
                                                return = "gene.ranked"
    )
    
    # Use the most important genes to weight the components score
    ICA_scores_weighted = deconica::get_scores(ICA_deconv$X, weighted.list, summary = "weighted.mean", na.rm = TRUE)
    
    
    
    # Extract the proportion from the weighted scores
    tmp_dat = t(ICA_scores_weighted)
    colnames(tmp_dat) = colnames(rna_data)
    # tmp_rna = deconica::stacked_proportions_plot(tmp_dat)
    A_rna = abs(tmp_dat) %*% diag(1/colSums(abs(tmp_dat)))
    #A_rna = matrix(tmp_rna)
    
    A_matrix <- A_rna
    T_matrix <- ICA_deconv$S
  }
  
  if (method == "PREDE") {
    
    mat = as.matrix(mix_matrix)
    
    feat = PREDE::select_feature(mat = mat ,method = "cv",nmarker = 1000,startn = 0)
    
    pred <- PREDE::PREDE(mat[feat,], W1=NULL,type = "GE",K=k,iters = 100,rssDiffStop=1e-5)
    
    A_matrix <- pred$H
    T_matrix <- pred$W
    
  }
 
  
  return( list(A_matrix = A_matrix,T_matrix = T_matrix) )
  
}
