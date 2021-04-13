#' Function to run RNAseq count normalisation methods.
#'
#' This function run a RNAseq count normalisation methods and takes as input a raw gene expression matrix
#'
#' available methods are "DESeq2", "TPM", "RPM", "CPM", "edgeR"
#'
#' TPM  and RKPM methods require gene length argument
#'
#' #'pipes from ( ... %>% )  :  [data.matrix]
#'
#' pipes to ( %>% ... ) :  [run_trans] , [run_deconv]
#'
#' returns a matrix
#'
#'
#' @param mix_matrix The raw gene count matrix samples *x* genes.
#' @param method The normalization method to be used
#' @param gene_length_bp The gene length in bp
#'
#' @return This function return a matrix  samples *x* genes
#'
#'
#' @importFrom DESeq2 "estimateSizeFactorsForMatrix"
#' @importFrom edgeR "DGEList"
#' @importFrom edgeR "estimateCommonDisp"
#' @importFrom edgeR "calcNormFactors"
#'
#' @export

run_norm <- function(mix_matrix, method, gene_length_bp = NULL, group = NULL) {
  if (!{
    method %in% c("DESeq2", "TPM", "RPM", "CPM", "edgeR", "MR")
  }) {
    print("Unknown method argument, please specify a method within the following list: DSEq2, edgeR, TPM")
  }

  if (is.null(group)) group <- rep("a", ncol(mix_matrix))

  if (method == "DESeq2") {
    size.factor <- DESeq2::estimateSizeFactorsForMatrix(mix_matrix) ## first calculate the size factors
    norm.counts <- sweep(mix_matrix, 2, size.factor, "/") ### divide by tje SF
  }
  if (method == "MR") {
    count_mtx <- mix_matrix
    # estimsF is  DESq2::estimateSizeFactorsForMatrix
    # estimSf <- function (count_mtx){
    #   # Compute the geometric mean
    geomMean <- function(x) prod(x)^(1 / length(x))
    # Compute the geometric mean over the line
    gm.mean <- apply(count_mtx, 1, geomMean)
    # Zero values are set to NA (avoid subsequentcdsdivision by 0)
    gm.mean[gm.mean == 0] <- NA
    # Divide each line by its corresponding geometric mean
    cts <- sweep(count_mtx, 1, gm.mean, FUN = "/")
    # Compute the median over the columns
    sFactor <- apply(cts, 2, median, na.rm = TRUE)
    # Return the scaling factor
    #   return(sFactor)
    # }
    # size.factor <- estimSf(mix_matrix) ## first calculate the size factors
    norm.counts <- sweep(mix_matrix, 2, sFactor, "/") ### divide by tje SF
  }

  else if (method %in% c("RPM", "CPM")) {
    norm.counts <- sweep(mix_matrix * 1e6, 2, colSums(mix_matrix), FUN = "/")
    ## CPM(RPM) = TotalReadPerGene*1e6 / TotalMappedRead
  }

  else if (method %in% c("TPM")) {
    A <- mix_matrix * 1e3 / gene_length_bp
    ## A = TotalReadPerGene*1e3/GeneLengthBp
    norm.counts <- sweep(A * 1e6, 2, colSums(mix_matrix), FUN = "/")
    ## TPM = A * 1/ TotalMappedRead*1e6
  }
  else if (method %in% c("RKPM")) {
    # norm.counts <- mix_matrix*1e9 / (sum(mix_matrix)*gene_length_bp)
    norm.counts <- sweep(mix_matrix * 1e9, 2, colSums(mix_matrix) * gene_length_bp, FUN = "/")
    ## RKPM =TotalReadPerGene*1e9 / (TotalMappedRead*GeneLengthBp)
  }
  else if (method %in% c("edgeR")) {
    d <- edgeR::DGEList(counts = mix_matrix, group = group)
    d <- edgeR::estimateCommonDisp(d)
    d <- edgeR::calcNormFactors(d, method = "TMM")
    norm.counts <- edgeR::cpm(d)
  }
  return(norm.counts)
}

# # from https://doi.org/10.1186/s12859-018-2246-7 supp File 4 :
# # geneID from Ensembl  used as row names in data matrix (x)
# # first column in (x) holds the gene length in kb and the
# # remaining columns contain read counts of each sample.
# #
# # calculate RPK
# rpk <- (x[,2:ncol(x)]/x[,1])
# # remove length col in x
# x <- x[,-1]
# # for normalization purposes, no grouping of samples
# group <- c(rep("A",ncol(x)))
# #EdgeR
# x.norm.edger <- DGEList(counts=x,group=group)
# x.norm.edger <- calcNormFactors(x.norm.edger)
# norm.counts.edger <- cpm(x.norm.edger)
#
# #GeTMM
# rpk.norm <- DGEList(counts=rpk,group=group)
# rpk.norm <- calcNormFactors(rpk.norm)
# norm.counts.rpk_edger <- cpm(rpk.norm)
#
# #TPM
# tpm = rpk
# for (i in 1:ncol(rpk) ) {
#   tpm[,i] <- rpk[,i]/(sum(rpk[,i])/1e6)
# }

# #DESeq2
# # no group & no design implemented
# colData = data.frame(group)
# rownames(colData)=colnames(x)
# dds<-DESeqDataSetFromMatrix(countData=x,colData=colData, design=~ 1)
# dds <- estimateSizeFactors(dds)
# sizefact <- sizeFactors(dds)
# norm.counts.deseq <- counts(dds, normalized=TRUE)
