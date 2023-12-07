#' Function to run deconvolution methods.
#'
#' This function run a deconvolution  methods and takes as input a pre-treated gene expression matrix
#'
#' available methods are "NMF", "ICA", "CDSeq", "PREDE", "DSA
#'
#' pipes from ( ... %>% )  : [data.matrix], [run_norm] , [run_trans]
#'
#' pipes to ( %>% ... ) : [comp_res], [score_res]
#'
#' returns a list(A_matrix, T_matrix)
#'
#'
#' @param mix_matrix The convoluted (pre-treated) matrix patients *x* genes.
#' @param k The number of cell types to estimate
#' @param method The deconvolution method to be used
#'
#' @return return a list with :
#' A_matrix a proportion matrix   cellType *x* samples and
#' T_matrix  a cellType reference matrix  genes *x* cellType .
#'
#' @importFrom NMF "nmf"
#' @importFrom fastICA "fastICA"
#' @importFrom parallel "detectCores"
#' @importFrom MASS "ginv"
#' @export
#'
run_deconv <-
  function(mix_matrix,
           k = 5,
           method = "NMF",
           gene_length = NULL,
           gene_id = NULL,
           cpu_number = NULL,
           ICA_orient=FALSE) {
    if (!{
      method %in% c("NMF", "ICA", "ICA-deconica", "CDSeq", "PREDE","debCAM")
    }) {
      print(
        "Unknown method argument, please specify a method within the following list: NMF, ICA, ICA-deconica, debCAM, CDSeq, PREDE"
      )
    }

    if (method == "NMF") {
      # library(NMF)
      # detach("package:DelayedArray")
      res <-
        NMF::nmf(
          x = mix_matrix,
          rank = k,
          method = "snmf/r",
          seed = 1
        )
      A <- apply(
        X = res@fit@H,
        MARGIN = 2,
        FUN = function(x) {
          x / sum(x)
        }
      )
      colnames(A)=colnames(mix_matrix)
      A_matrix <- A
      #T_matrix <- res@fit@W
      remove(list = "res")
      row.names(A_matrix) <- paste("NC", 1:k, sep = "")
      T_matrix=data.matrix(mix_matrix) %*% MASS::ginv(data.matrix(A))
      T_matrix[T_matrix<0]=0
    }

    if (method == "ICA") {
      # library(NMF)
      # detach("package:DelayedArray")
      rna_data <- mix_matrix[!duplicated(mix_matrix), ]

      ICA_deconv <-
        fastICA::fastICA(
          X = rna_data,
          n.comp = k,
          maxit = 1000,
          tol = 1e-09
        )
      ICA_deconv$names <- row.names(rna_data)
      ## generate_markers
      ##  deconica::generate_markers(df = ICA_deconv,
      ## deconica::orient_func
      sel.comp = paste("IC", 1:ncol(ICA_deconv$S),sep="")
      n=30 # default deconica value
      thr= Inf# default deconica value
      S <- ICA_deconv$S
      orient <-
        apply(S, 2, function(x) {
          if (min(x) < -3 & max(x) > 3) {
            ifelse(sum(x > 3) < sum(x < -3), -1, 1)
          } else {
            ifelse(sum(x > 2) < sum(x < -2), -1, 1)
          }
        })
      S_or <- as.matrix(S) %*% diag(orient)
      colnames(S_or) <- paste("IC", 1:ncol(S_or), sep = "")
      row.names(S_or) <- ICA_deconv$names
      S_or <- S_or[, sel.comp]
      metagenes <-
        apply(S_or, 2, function(col) {
          data.frame(
            GENE = row.names(S_or),
            col
          )
        })
      weight.list <- lapply(metagenes, function(x) {
        x <- x[order(-x[, 2]), ]
        return(x[which(x[, 2] < thr), ][1:n, ])
      })
      #  Use the most important genes to weight the components score
      #  deconica::get_scores
      ICA_scores_weighted <- sapply(weight.list, function(metagene) {
        apply(ICA_deconv$X[metagene[, 1], ], 2, stats::weighted.mean, w = metagene[, 2])
      })
      tmp_dat <- t(ICA_scores_weighted)
      colnames(tmp_dat) <- colnames(rna_data)
      #  deconica::stacked_proportions_plot(tmp_dat)
      A_rna <- abs(tmp_dat) %*% diag(1 / colSums(abs(tmp_dat)))
      # A_rna = matrix(tmp_rna)
      colnames(A_rna)=colnames(mix_matrix)
      A_matrix <- A_rna
      if(ICA_orient) T_matrix <- S_or
      else  T_matrix <- S

      # OTHER APPROACH with NMF
      # Ap=ICA_deconv$A ; Ap[ICA_deconv$A<0] =0
      # Am=-ICA_deconv$A ; Am[ICA_deconv$A>0] =0
      #
      # Aa=rbind(Ap,Am)
      # Aa=Aa[apply(Aa, 1, var)>0 ,]
      # res <-
      #   NMF::nmf(
      #     x = Aa,
      #     rank = k,
      #     #method = "snmf/r",
      #     method = "lee",
      #     seed = 1
      #   )
      # A <- apply(
      #   X = res@fit@H,
      #   MARGIN = 2,
      #   FUN = function(x) {
      #     x / sum(x)
      #   }
      # )
      #
      # A_matrix <- A
      # Ainv=MASS::ginv(A_matrix)
      # T_matrix= mix_matrix %*% Ainv
      # #T_matrix <- res@fit@W
      # remove(list = "res")
      # row.names(A_matrix) <- paste("NC", 1:k, sep = "")
      # colnames(T_matrix) <- paste("NC", 1:k, sep = "")
    }

    if (method == "ICA-deconica") {
      if (!requireNamespace("deconica", quietly = TRUE)) {
        stop("Package \"deconica\" needed for this function to work. Please install it.",
          call. = FALSE
        )
      }
      rna_data <- mix_matrix[!duplicated(mix_matrix), ]

      ICA_deconv <-
        fastICA::fastICA(
          X = rna_data,
          n.comp = k,
          maxit = 1000,
          tol = 1e-09
        )
      ICA_deconv$names <- row.names(rna_data)
      weighted.list <- deconica::generate_markers(
        df = ICA_deconv,
        n = 30,
        return = "gene.ranked"
      )

      # Use the most important genes to weight the components score
      ICA_scores_weighted <-
        deconica::get_scores(ICA_deconv$X,
          weighted.list,
          summary = "weighted.mean",
          na.rm = TRUE
        )



      # Extract the proportion from the weighted scores
      tmp_dat <- t(ICA_scores_weighted)
      colnames(tmp_dat) <- colnames(rna_data)
      # tmp_rna = deconica::stacked_proportions_plot(tmp_dat)
      A_rna <- abs(tmp_dat) %*% diag(1 / colSums(abs(tmp_dat)))
      # A_rna = matrix(tmp_rna)
      colnames(A_rna)=colnames(mix_matrix)
      A_matrix <- A_rna
      T_matrix <- ICA_deconv$S
    }

    if (method == "PREDE") {
      if (!requireNamespace("PREDE", quietly = TRUE)) {
        stop("Package \"PREDE\" needed for this function to work. Please install it.",
          call. = FALSE
        )
      }
      mat <- as.matrix(mix_matrix)

      # feat = PREDE::select_feature(mat = mat ,method = "cv",nmarker = 1000,startn = 0)

      # pred <- PREDE::PREDE(mat[feat,], W1=NULL,type = "GE",K=k,iters = 100,rssDiffStop=1e-5)
      pred <- PREDE::PREDE(
        mat,
        W1 = NULL,
        type = "GE",
        K = k,
        iters = 100,
        rssDiffStop = 1e-5
      )

      A_matrix <- pred$H
      T_matrix <- pred$W
      row.names(A_matrix) <- paste("PREDE_C", 1:k, sep = "")
    }

    if (method == "CDSeq") {
      if (!requireNamespace("CDSeq", quietly = TRUE)) {
        stop("Package \"CDSeq\" needed for this function to work. Please install it.",
          call. = FALSE
        )
      }
      # library(CDSeq)
      # samples_id= colnames(test_data_rna[[1]])
      # row.names(mix_matrix) <- gene_id
      # colnames(mix_matrix) <- samples_id
      nsz <- ceiling(nrow(mix_matrix) * 1e-3 / 8)
      nblock <- ceiling(nrow(mix_matrix) / (nsz * 1e3))
      redFact <- 2^(1 + (median(log2(1 + mix_matrix[mix_matrix > 0])) %/% 5))
      if (nblock>1) {
      print(
        sprintf(
          "%d var in %d blocks of size %d with reduce factor %d",
          nrow(mix_matrix),
          nblock,
          nsz * 1e3,
          redFact
          
        )
      )
        gene_subset_size = nsz * 1e3
      } else {
        nblock = NULL
        gene_subset_size = NULL
      }
      result1 <- CDSeq::CDSeq(
        bulk_data = mix_matrix,
        cell_type_number = k,
        beta = 0.5,
        alpha = 5,
        mcmc_iterations = 300,
        cpu_number = cpu_number,
        dilution_factor = redFact,
        block_number = nblock,
        gene_subset_size =  gene_subset_size,
        gene_length = as.vector(gene_length)
      )

      A_matrix <- result1$estProp
      T_matrix <- result1$estGEP
    }
    
    if (method == "debCAM") {
      if (!requireNamespace("debCAM", quietly = TRUE)) {
        stop("Package \"debCAM\" needed for this function to work. Please install it.",
             call. = FALSE
        )
      }
      cluster.num = min(5*k, ncol(mix_matrix) - 1 )
      if (nrow(mix_matrix) < 200) {
        dim.rdc =max(cluster.num, nrow(mix_matrix)/10 ) 
      }else dim.rdc= 10
      rCAM <- debCAM::CAM(data =  mix_matrix,
                          K = k,
                          cluster.num =cluster.num,
                          MG.num.thres = 1,
                          lof.thres =0,
                          dim.rdc =  dim.rdc
                          )
                          #thres.low = 0, thres.high = 1)

      A_matrix <- t(debCAM::Amat(rCAM, k))
      T_matrix <- debCAM::Smat(rCAM, k)
    }
    
    # compute soft param for CDSeq
    nsz <- ceiling(nrow(mix_matrix) * 1e-3 / 8)
    nblock <- ceiling(nrow(mix_matrix) / (nsz * 1e3))
    redFact <- 2^(1 + (median(log2(1 + mix_matrix[mix_matrix > 0])) %/% 5))
    if (nblock>1) {
      print(
        sprintf(
          "%d var in %d blocks of size %d with reduce factor %d",
          nrow(mix_matrix),
          nblock,
          nsz * 1e3,
          redFact
          
        )
      )
      gene_subset_size = nsz * 1e3
    } else {
      nblock = NULL
      gene_subset_size = NULL
    }
    #compute soft param for debCAM
    cluster.num = min(5*k, ncol(mix_matrix) - 1 )
    if (nrow(mix_matrix) < 200) {
      dim.rdc =max(cluster.num, nrow(mix_matrix)/10 ) 
    }else dim.rdc= 10
    
    
    # return res object 
    return(list(A_matrix = A_matrix, 
                T_matrix = T_matrix, 
                run_params=list(
                  "NMF"=list(method = "snmf/r",
                             seed = 1),
                  "ICA"=list(maxit = 1000,
                             tol = 1e-09),
                  "ICA-deconica"= list (score.meth = weighted.list,
                                        summary.score = "weighted.mean"),
                                         
                  "PREDE" = list(W1 = NULL,
                                 type = "GE",
                                 K = k,
                                 iters = 100,
                                 rssDiffStop = 1e-5),
                  "CDSeq"=list (  beta = 0.5,
                                  alpha = 5,
                                  mcmc_iterations = 300,
                                  cpu_number = cpu_number,
                                  dilution_factor = redFact,
                                  block_number = nblock,
                                  gene_subset_size =  gene_subset_size,
                                  gene_length = as.vector(gene_length)),
                  "debCAM"= list( cluster.num =cluster.num,
                                  MG.num.thres = 1,
                                  lof.thres =0,
                                  dim.rdc =  dim.rdc)  
                  )))
  }
