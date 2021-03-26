#'Function to run RNAseq count normalisation methods.
#'
#'This function run a RNAseq count normalisation methods and takes as input a raw gene expression matrix
#'
#'available methods are "DESeq2", "TPM", "RPM", "CPM", "edgeR"
#'
#'TPM  and RKPM methods require gene length argument
#'
#'#'pipes from ( ... %>% )  :  [data.matrix]
#'
#'pipes to ( %>% ... ) :  [run_trans] , [run_deconv]   
#'
#'returns a matrix 
#'
#'
#'@param mix_matrix The raw gene count matrix samples *x* genes.
#'@param meta The count matrix metadata
#'@param method The normalization method to be used 
#'
#'@return This function return a matrix  samples *x* genes
#'
#'
#'@importFrom DESeq2 "estimateSizeFactorsForMatrix"
#'
#'@export

run_norm <- function(mix_matrix, method, gene_length_bp=NULL) {
  if ( !{ method %in% c("DESeq2", "TPM", "RPM", "CPM", "edgeR") } ) {
    print("Unknown method argument, please specify a method within the following list: DSEq2, edgeR, TPM")
  }
  if (method == "DESeq2") {
  size.factor <- DESeq2::estimateSizeFactorsForMatrix(mix_matrix) ## first calculate the size factors
  norm.counts <- sweep(mix_matrix, 2, size.factor, "/")  ### divide by tje SF
  #norm.counts.pseudoc.count.log2 <- log2(norm.counts + 1) ## If you nedd the pseudocounts ...
  # estimsF is  DESq2::estimateSizeFactorsForMatrix
  # estimSf <- function (count_mtx){
  #   # Compute the geometric mean
  #   geomMean <- function(x) prod(x)^(1/length(x))
  #   # Compute the geometric mean over the line
  #   gm.mean <- apply(count_mtx, 1, geomMean)
  #   # Zero values are set to NA (avoid subsequentcdsdivision by 0)
  #   gm.mean[gm.mean == 0] <- NA
  #   # Divide each line by its corresponding geometric mean
  #   cts <- sweep(count_mtx, 1, gm.mean, FUN="/")
  #   # Compute the median over the columns
  #   sFactor <- apply(cts, 2, median, na.rm=TRUE)
  #   # Return the scaling factor
  #   return(sFactor)
  # }
  }
  
  else if (method %in%  c("RPM", "CPM")) {
    norm.counts <- sweep(mix_matrix*1e6 , 2, colSums(mix_matrix), FUN="/") 
    ## CPM(RPM) = TotalReadPerGene*1e6 / TotalMappedRead
    
  }
  
  else if (method %in%  c("TPM")) {
    A= mix_matrix*1e3 / gene_length_bp  
    ## A = TotalReadPerGene*1e3/GeneLengthBp
    norm.counts <- sweep(A*1e6 , 2, colSums(mix_matrix), FUN="/")      
    ## TPM = A * 1/sum(A) *1e6
  }
  else if (method %in%  c("RKPM")) {
    # norm.counts <- mix_matrix*1e9 / (sum(mix_matrix)*gene_length_bp) 
    norm.counts <- sweep(mix_matrix*1e9 , 2, colSums(mix_matrix)*gene_length_bp, FUN="/")
    ## RKPM =TotalReadPerGene*1e9 / (TotalMappedRead*GeneLengthBp)
    
  }
  
  return(norm.counts)
  }


