#' Function to plot deconvolution and reference matrices.
#'
#' This function compares a deconvolution result with its reference
#'
#' #'pipes from ( ... %>% )  :  [run_deconv]nts <- sweep(mix_matrix*1e6 , 2, colSums(mix_matrix), FUN="/")
#'
#' pipes to ( %>% ... ) :  [score_res]
#'
#' with pipe=TRUE, returns a list(A_matrix, T_matrix)
#'
#'
#' @param res res, a list(A_matrix, T_matrix)
#' @param ref, a list(A_matrix, T_matrix)
#'
#' @return With pipe=TRUE (default) this function return res, a list(A_matrix, T_matrix)
#'
#' @importFrom clue "solve_LSAP"
#' @importFrom pheatmap "pheatmap"
#'
#' @export

comp_res <- function(res = res_, ref = ref_, pipe = TRUE) {
  A_est <- res$A_matrix
  A_ref <- ref$A_matrix
  T_est <- res$T_matrix
  T_ref <- ref$T_matrix
  cmat_A <- cor(
    t(A_est),
    t(A_ref)
  )
  cor_plot(cmat_A, main = "A cor")
  cmat_T <- cor(
    T_est,
    T_ref
  )
  row.names(cmat_T) <- row.names(cmat_A)
  cor_plot(cmat_T, main = "T cor")
  if (pipe) {
    return(res_ = res)
  }
}

cor_plot <- function(c_mat = c_mat_, ...) {
  ord_c <- c(clue::solve_LSAP((1 + c_mat)^2, maximum = TRUE))
  # ord_c=clue::solve_LSAP(c_mat-min(c_mat)
  #                      ,maximum = TRUE
  # )
  pheatmap::pheatmap(c_mat[, ord_c],
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    cellwidth = 20,
    cellheight = 20,
    ...
  )
}

#' Function to evaluate deconvolution and reference matrices.
#'
#' This function evaluates a deconvolution result with its reference
#'
#' #'pipes from ( ... %>% )  :  [run_deconv]
#'
#' pipes to ( %>% ... ) :  [comp_res]
#'
#' with pipe=TRUE, returns a list(A_matrix, T_matrix)
#'
#'
#' @param res res, a list(A_matrix, T_matrix)
#' @param ref, a list(A_matrix, T_matrix)
#'
#' @return With pipe=TRUE (default) this function return res, a list(A_matrix, T_matrix)
#'
#'
#'
#' @export
score_res <- function(res = res_, ref = ref_, pipe = TRUE) {
  sc <- scoring_function(
    Aest = res$A_matrix,
    Aref = ref$A_matrix,
    Tref = ref$T_matrix
  )
  if (pipe) {
    print(sc)
    return(res_ = res)
  }
  else {
    return(sc)
  }
}


##########################
### SCORING FUNCTIONS ###
#########################

#########################
# merging fonction when too many estimated cell types are estimated

#' Title
#'
#' @param A
#' @param k
#' @param plot
#' @param unknown_est
#'
#' @return
#' @export
#'
merge_hc <- function(A, k, T = NULL, plot = NULL, unknown_est = NULL) {
  if (k >= nrow(A)) k <- nrow(A)

  if (is.null(T)) {
    hc <- hclust(d = as.dist(1 - cor(t(A))))
  } else {
    hc <- hclust(d = dist(scale(t(T))))
  }

  dend <- as.dendrogram(hc)

  if (!is.null(unknown_est)) {
    bestk_clust <- cutree(hc, k - 1)
    bestk_clust[unknown_est] <- k
  } else {
    bestk_clust <- cutree(hc, k)
  }

  print(c(plot, bestk_clust))
  dend <- as.dendrogram(hc)

  Amerge <- data.matrix(do.call(
    rbind,
    lapply(
      unique(bestk_clust),
      function(icol) {
        idx <- bestk_clust == icol
        if (sum(idx) > 1) {
          return(apply(
            A[idx, ],
            2,
            sum
          ))
        } else {
          return(A[idx, ])
        }
      }
    )
  ))

  return(Amerge)
}

#########################
# homogenization function to find the best match between real and estimated A matrix

#' Title
#'
#' @param A_r
#' @param A_est
#'
#' @return
#'
#' @importFrom clue "solve_LSAP"
#'
#' @export
#'
homgeneized_cor_mat <- function(A_r, A_est) {
  cmat <- cor(t(A_r), t(A_est))
  pvec <- c(clue::solve_LSAP((1 + cmat)^2, maximum = TRUE))
  return(A_est[pvec, ])
}

# Estimated A pre-treatment
#' Title
#'
#' @param A_r
#' @param A_est
#' @param T_ref
#'
#' @return
#' @export
#'
prepare_A <- function(A_r, A_est, T_ref) {
  N <- ncol(x = A_r)
  K <- nrow(x = A_r)
  stopifnot(K > 1)

  stopifnot(ncol(x = A_est) == N)
  stopifnot(!anyNA(x = A_est))

  ### STEP 1 : matching the number of estimated component to real number of cell types K

  ## remove rows at 0
  idx_to_keep <- which(rowSums(A_est) > 0)
  A_est <- A_est[idx_to_keep, ]


  ## if not supplying enough types (make sure that {nrow(A_est) >= K})

  if (nrow(x = A_est) < K) {
    # set positive random values closed to 0 for missing rows
    set.seed(1)
    random_data <- abs(jitter(matrix(data = 0, nrow = K - nrow(x = A_est), ncol = N), factor = 0.01))
    A_est <- rbind(A_est, random_data)
    print("Add rows of 0 to match the exact number of K")
  }

  ## if supplying too manycell types

  if (nrow(x = A_est) > K) {
    print("Number of cell type exceded K, filtering and clustering are applies")
    filteredout <- rowMeans(A_est) > 3e-2
    print(which(!filteredout))
    tokeep <- setdiff(which(filteredout), nrow(A_est))

    ### Case 1 : Nb of filtered component are lower that K -> we add positive random values closed to 0 for missing rows
    if (sum(filteredout) < K) {
      A_est <- A_est[order(rowSums(A_est), decreasing = TRUE)[1:K], ]
      print("Number of cell type exceded K, only the K most contributing cell types were kept for scoring.")
    }

    ### Case 2 : Nb of filtered component are equel to K
    if (sum(filteredout) == K) {
      print("Nb of kept component is equal K")
      A_est <- A_est[filteredout, ]
    }

    ### Case 3 : Nb of filtered component are greater than K -> we apply hierarchical clustering to aggregate similar components based of the reference profiles (T matrix)
    if (sum(filteredout) > K) {
      print("Nb of kept component greater than equal K, we apply a clustering step")
      T_ref <- apply(T_ref, 2, as.numeric)
      nb_clust <- nrow(A_r)
      tokeep <- setdiff(which(filteredout), nrow(A_est))
      unchar <- NULL
      if (nrow(A_est) - ncol(T_ref) > 0) {
        unchar <- sum(filteredout)
      } # handle the case in which an uncharacterized cell type is estimated by a reference based method
      A_est <- merge_hc(
        A = A_est[filteredout, ],
        k = nb_clust,
        T = T_ref[, tokeep],
        unknown_est = unchar
      )
    }
  }

  ### STEP 2 : reordering estimated components to find the best match

  A_est <- t(t(A_est) / colSums((A_est))) # set col sum to 1 if not already
  A_est_best_cor <- as.matrix(homgeneized_cor_mat(A_r, A_est))

  ### Return the pre-treated estimated A matrix, an ordered matrix with correct number of estimated components

  return(A_est_best_cor)
}

#########################
# Coputation of column correlation score

#' Title
#'
#' @param A_r
#' @param Aest_p
#'
#' @return
#' @export
#'
correlation_col <- function(A_r, Aest_p) {
  res <- c()
  for (i in 1:ncol(A_r)) {
    if (sd(Aest_p[, i]) > 0 & sd(A_r[, i]) > 0) {
      res[i] <- cor(A_r[, i], Aest_p[, i], method = "pearson")
    }
  }

  res <- res[!is.na(res)]
  print(paste0(length(res), " cols are kept for correlation analysis"))
  res[which(res < 0)] <- 0
  COR <- sum(res) / length(res)

  return(COR)
}

#########################
# Coputation of row correlation score

correlation_row <- function(A_r, Aest_p) {
  res <- c()
  for (i in 1:nrow(A_r)) {
    if (sd(Aest_p[i, ]) > 0 & sd(A_r[i, ]) > 0) {
      res[i] <- cor(A_r[i, ], Aest_p[i, ], method = "pearson")
    }
  }

  res <- res[!is.na(res)]
  print(paste0(length(res), " rows are kept for correlation analysis"))
  res[which(res < 0)] <- 0
  COR <- sum(res) / length(res)

  return(COR)
}

#########################
# Computation of Mean absolute error score

eval_MAE <- function(A_r, Aest_p) {
  #' Title
  #'
  #' @param M1
  #' @param M2
  #'
  #' @return
  #' @export
  #'
  MAE <- function(M1, M2) {
    return(mean(x = abs(x = M1 - M2)))
  }
  return(MAE(A_r, Aest_p))
}

#########################
# Scoring function

#' Title
#'
#' @param Aref
#' @param Aest
#' @param Tref
#'
#' @return
#' @export
#'
scoring_function <- function(Aref, Aest, Tref) {
  #  pretreatment of estimated A
  Aest_p <- prepare_A(A_r = Aref, A_est = Aest, T_ref = Tref)
  #  scoring
  mae <- eval_MAE(Aref, Aest_p)
  cr <- correlation_row(Aref, Aest_p)
  cc <- correlation_col(Aref, Aest_p)
  #  scoring agregation (using a derivative of the maxmin approach)
  rd_mae <- c()
  set.seed(1)
  random_col <- c(1, rep(0, (nrow(Aref) - 1)))
  random_base <- matrix(rep(random_col, ncol(Aref)), nrow(Aref), ncol(Aref))
  for (j in 1:1000) {
    rd <- random_base[sample(nrow(Aref)), ]
    rd_mae[j] <- eval_MAE(Aref, rd)
  }
  max_mae <- max(rd_mae)
  min_mae <- 0
  max_cor <- 1
  min_cor <- 0
  mae_maxmin <- (mae - min_mae) / (max_mae - min_mae)
  cr_maxmin <- (cr - min_cor) / (max_cor - min_cor)
  cc_maxmin <- (cc - min_cor) / (max_cor - min_cor)
  score_combine <- (cr_maxmin + cc_maxmin + (1 - mae_maxmin)) / 3
  return(list(
    score_combine = score_combine,
    mae = mae,
    cr = cr,
    cc = cc,
    mae_maxmin = mae_maxmin,
    cr_maxmin = cr_maxmin,
    cc_maxmin = cc_maxmin
  ))
}
