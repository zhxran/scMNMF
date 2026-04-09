#' @import Rcpp
#' @import scInt
#' @importFrom stats model.matrix setTxtProgressBar txtProgressBar
#' @useDynLib scMNMF, .registration = TRUE
NULL

#' Compute one-hot matrix for given data frame and variable(s)
#'
#' @param x Input data frame.
#' @param ivar Variable(s) for one-hot computation.
#' @return A one-hot encoded matrix.
#' @export
group_onehot <- function(x, ivar) {
  if (length(unique(x[, ivar])) == 1) {
    matrix(1, nrow = length(x[, ivar]), ncol = 1)
  } else {
    x_df <- data.frame(ivar = x[, ivar])
    x_reduced <- Reduce(paste0, x_df)
    model.matrix(~ 0 + x_reduced)
  }
}

#' Single-cell Modular Non-negative Matrix Factorization (scMNMF)
#'
#' A modular NMF-based framework for scRNA-seq integration, decoupling shared biological signals 
#' from condition-specific effects.
#'
#' @param X A normalized gene-by-cell expression matrix.
#' @param meta A data frame containing cell metadata.
#' @param batch Column name in meta for batch info.
#' @param cell_type Column name in meta for cell type identity.
#' @param condition Column name in meta for biological conditions.
#' @param k1 Number of shared factors. Default is 20.
#' @param k2 Number of condition-specific factors. Default is 5.
#' @param lambda_1 Hyperparameter for cell-type graph regularization. Default is 5.
#' @param lambda_2 Hyperparameter for condition graph regularization. Default is 5.
#' @param eta Hyperparameter for batch effect size control. Default is 5.
#' @param gamma Hyperparameter for sparsity. Default is 5.
#' @param thresh Convergence threshold. Default is 1e-6.
#' @param max.iters Maximum number of iterations. Default is 100.
#'
#' @return A list containing matrices: W, H, C, and D.
#' @export
run_scMNMF <- function(
    X, meta, batch, cell_type, condition,
    k1 = 20, k2 = 5, 
    lambda_1 = 5, lambda_2 = 5, eta = 5, gamma = 5,
    thresh = 1e-6, max.iters = 100
) {
  # --- Step 0: Handle Internal scInt Dependencies ---
  # We use the internal functions from scInt for high-performance matrix operations
  eigenMapMatMult <- scInt:::eigenMapMatMult
  eigenMapMatcrossprod <- scInt:::eigenMapMatcrossprod
  eigenMapMattcrossprod <- scInt:::eigenMapMattcrossprod

  # --- Step 1: Initialization ---
  set.seed(123)
  tmp <- gc() 
  M <- nrow(X) 
  N <- ncol(X) 
  
  W <- matrix(data = abs(runif(M * k1, 0, 2)), nrow = M, ncol = k1)
  C <- matrix(data = abs(runif(M * k2, 0, 2)), nrow = M, ncol = k2)
  rownames(W) <- rownames(C) <- rownames(X)
  
  H <- matrix(data = abs(runif(k1 * N, 0, 2)), nrow = k1, ncol = N)
  D <- matrix(data = abs(runif(k2 * N, 0, 2)), nrow = k2, ncol = N)
  
  # --- Step 2: Graph Construction ---
  OH_t <- group_onehot(meta, cell_type)
  OH_bt <- group_onehot(meta, c(cell_type, batch))
  OH_c <- group_onehot(meta, condition)
  OH_bc <- group_onehot(meta, c(condition, batch))
  
  # createK is provided by the compiled C++ code in src/runK.cpp
  K_bt_plus <- createK(OH_bt)
  K_bc_plus <- createK(OH_bc)
  K_t <- createK(OH_t)
  K_c <- createK(OH_c)
  
  E <- matrix(1, nrow = M, ncol = M)
  X <- as.matrix(X)
  
  # --- Step 3: Iterative Optimization ---
  delta <- 1
  iters <- 0
  pb <- txtProgressBar(min = 0, max = max.iters, style = 3)
  
  # Initial Objective function
  obj0 <- norm(X - eigenMapMatMult(W, H) - eigenMapMatMult(C, D), type = "F") ^ 2 +
    sum(sapply(1:ncol(OH_bt), function(i) {
      which_i = which(OH_bt[, i] == 1)
      lambda_1 * sum(rowSums(H[, which_i, drop = FALSE])^2) / length(which_i)
    })) +
    sum(sapply(1:ncol(OH_bc), function(i) {
      which_i = which(OH_bc[, i] == 1)
      lambda_2 * sum(rowSums(D[, which_i, drop = FALSE])^2) / length(which_i)
    })) -
    sum(sapply(1:ncol(OH_t), function(i) {
      which_i = which(OH_t[, i] == 1)
      lambda_1 * sum(rowSums(H[, which_i, drop = FALSE])^2) / length(which_i)
    })) -
    sum(sapply(1:ncol(OH_c), function(i) {
      which_i = which(OH_c[, i] == 1)
      lambda_2 * sum(rowSums(D[, which_i, drop = FALSE])^2) / length(which_i)
    })) + 
    eta * norm(eigenMapMatMult(C, D), type = "F") ^ 2 +
    gamma * sum(apply(abs(C), 2, function(x) sum(x)^2))
  
  while (delta > thresh & iters < max.iters) {
    # Updates
    W <- W * eigenMapMattcrossprod(X, H) / 
      (eigenMapMatMult(W, eigenMapMattcrossprod(H, H)) + eigenMapMattcrossprod(eigenMapMatMult(C, D), H) + 1e-8)
    
    H <- H * (eigenMapMatcrossprod(W, X)  + lambda_1 * eigenMapMatMult(H, K_t)) / 
      (eigenMapMatMult(eigenMapMatcrossprod(W, W), H) + eigenMapMatcrossprod(W, eigenMapMatMult(C, D)) + lambda_1 * eigenMapMatMult(H, K_bt_plus) + 1e-8)
    
    C <- C * eigenMapMattcrossprod(X, D) / 
      (eigenMapMatMult(W, eigenMapMattcrossprod(H, D)) + (eta + 1)*eigenMapMatMult(C, eigenMapMattcrossprod(D, D)) + gamma * eigenMapMatMult(E,C) + 1e-8)
    
    D <- D * (eigenMapMatcrossprod(C, X) + lambda_2 * eigenMapMatMult(D, K_c)) / 
      (eigenMapMatMult(eigenMapMatcrossprod(C, W), H) + (eta + 1)*eigenMapMatMult(eigenMapMatcrossprod(C, C), D) + lambda_2 * eigenMapMatMult(D, K_bc_plus) + 1e-8)
    
    # Calculate Objective
    obj <- norm(X - eigenMapMatMult(W, H) - eigenMapMatMult(C, D), type = "F") ^ 2 +
      sum(sapply(1:ncol(OH_bt), function(i) {
        which_i = which(OH_bt[, i] == 1)
        lambda_1 * sum(rowSums(H[, which_i, drop = FALSE])^2) / length(which_i)
      })) +
      sum(sapply(1:ncol(OH_bc), function(i) {
        which_i = which(OH_bc[, i] == 1)
        lambda_2 * sum(rowSums(D[, which_i, drop = FALSE])^2) / length(which_i)
      })) -
      sum(sapply(1:ncol(OH_t), function(i) {
        which_i = which(OH_t[, i] == 1)
        lambda_1 * sum(rowSums(H[, which_i, drop = FALSE])^2) / length(which_i)
      })) -
      sum(sapply(1:ncol(OH_c), function(i) {
        which_i = which(OH_c[, i] == 1)
        lambda_2 * sum(rowSums(D[, which_i, drop = FALSE])^2) / length(which_i)
      })) + 
      eta * norm(eigenMapMatMult(C, D), type = "F") ^ 2 +
      gamma * sum(apply(abs(C), 2, function(x) sum(x)^2))
    
    delta <- abs(obj0 - obj) / (mean(c(obj0, obj)))
    obj0 <- obj
    iters <- iters + 1
    setTxtProgressBar(pb = pb, value = iters)
  }
  
  cat("\nMax iterations: ", max.iters, "\nFinal delta: ", delta, ".\n", sep = "")
  return(list("W" = W, "H" = H, "C" = C, "D" = D))
}