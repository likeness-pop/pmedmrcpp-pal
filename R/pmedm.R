#' P-MEDM Wrapper
#'
#' Wrapper for `PMEDMrcpp::pmedm_solve`. Adapted from various examples by Nagle
#' (2013, 2015). See `PMEDMrcpp::??pmedm_solve`.
#'
#'@param wt sample weights
#'@param cind individual constraints
#'@param cg2 block group (or target zone) constraints
#'@param cg1 (or aggregate zone) constraints
#'@param sg2 group (or target zone) constraint standard errors
#'@param sg1 tract (or aggregate zone) constraint standard errors
#'@param include_cg0 whether to include Level 0 (PUMA) constraints (default TRUE)
#'@param g1 an optional vector of aggregate zone IDs
#'@param g2 an optional vector of target zone IDs
#'@param tr_iter maximum number of inner trust region interations (optional)
#'@param tr_tol tolerance for stopping the solver (optional)
#'@return A list containing the P-MEDM allocation matrix, PUMS serials, and block group IDs.
#'@export
pmedm <- function(wt, cind, cg2, cg1, sg2, sg1, include_cg0 = TRUE, g1 = NULL, g2 = NULL, tr_iter = NULL, tr_tol = NULL){

  cind <- Matrix::Matrix(as.matrix(cind), sparse = TRUE)
  pX <- list(cind, cind)

  ## Topology
  if(!is.null(g1) & !is.null(g2)){
	  topo = data.frame(g2, g1)
  }else{
    g1_ids = sapply(rownames(cg2), function(g) substr(g, 1, nchar(g)-1))
          topo = data.frame(g2 = rownames(cg2), g1 = g1_ids)
  }

  ## Geographies
  if (include_cg0){
    A0 <- Matrix::t(rep(1, nrow(topo)))
  }else{
    A0 <- NULL
  }

  A1 <- do.call('rbind', lapply(unique(topo[,2]), function(g){
    nbt <- topo[topo[,2]==g,][,1]
    ifelse(topo[,1] %in% nbt,1,0)
  }))
  rownames(A1) <- unique(topo[,2])
  colnames(A1) <- topo[,1]
  A1 <- Matrix::Matrix(A1, sparse = TRUE)

  A2 <- do.call('rbind', lapply(topo[,1], function(g){
    ifelse(topo[,1] %in% g,1,0)
  }))
  rownames(A2) <- topo[,1]
  colnames(A2) <- topo[,1]
  A2 <- Matrix::Matrix(A2, sparse = TRUE)

  A <- list(A1, A2)
  names(A) <- c("lv1", "lv2")
  if (!is.null(A0)){
    A <- c(list("lv0" = A0), A)
  }

  ## Geographic Constraints
  if (include_cg0){
    cg0 <- data.frame(t(colSums(cg2)))
  }else{
    cg0 <- NULL
  }

  Y <- list(as.matrix(cg1), as.matrix(cg2))
  names(Y) <- c("lv1", "lv2")
  if(!is.null(cg0)){
    Y <- c(list("lv0" = as.matrix(cg0)), Y)
  }


  ## Error variances
  if (include_cg0){
    sg0 <- do.call(cbind, lapply(sg2, function(x) sqrt(sum(x^2))))
    sg0 <- sg0 * 0.1 # TODO: parameterize scaling factor
  }else{
    sg0 <- NULL
  }

  V <- list(as.matrix(sg1^2),as.matrix(sg2^2))
  names(V) <- c("lv1", "lv2")
  if (!is.null(sg0)){
    V <- c(list("lv0" = as.matrix(sg0^2)), V)
  }

  ## Generate PUMS Solver Inputs
  N <- sum(wt) # Population Size
  n <- nrow(pX[[1]]) # Sample Size

  # Since we are optimizing probabilities p, rather than weights w
  # Normalize Y (tract/bg data) by N (pop size) and V (tract/bg variances) by n/N^2
  Y_vec <- do.call('c', lapply(Y, function(x) as.vector(as.matrix(x)))) / N
  V_vec <- do.call('c', lapply(V, function(x) as.vector(as.matrix(x)))) * (n / N^2)

  # Will need a matrix V, not a vector V (variance-covariance matrix as sparse diagonal mat)
  sV <- Matrix::.sparseDiagonal(x = V_vec, shape = "g")

  # Solution space
  X <- rbind(
    kronecker(Matrix::t(pX[[1]]), A[["lv1"]]),
    kronecker(Matrix::t(pX[[2]]), A[["lv2"]])
  )
  if (!is.null(A0)){
    X <- rbind(kronecker(Matrix::t(pX[[1]]), A[["lv0"]]), X)
  }
  X <- Matrix::t(X)
  X <- Matrix::Matrix(X, sparse = TRUE)

  # Create design weights and normalize
  q <- matrix(wt, n, dim(A[["lv1"]])[2])
  q <- q / sum(as.numeric(q))
  q <- as.vector(Matrix::t(q))

  ## Solver options
  if ((!is.null(tr_iter)) | (!is.null(tr_tol))){
    pmedm_opts = PMEDMrcpp::get_PMEDM_opts()
    if (!is.null(tr_iter)){
      pmedm_opts[["tr_iter"]] <- tr_iter
    }
    if (!is.null(tr_tol)){
      pmedm_opts[["tr_tol"]] <- tr_tol
    }
  }else{
    pmedm_opts <- NULL
  }
                               
  ## Solve PMEDM Problem
  t <- PMEDMrcpp::PMEDM_solve(X, Y_vec, sV, q, opt = pmedm_opts)

  ## Output allocation matrix, response IDs, and zone IDs
  out <- list(
        almat = matrix(t$p * N, nrow(cind), dim(A[["lv1"]])[2], byrow = TRUE),
        serialno = rownames(cind),
        geoid = rownames(Y[[length(Y)]])
  )

  return(out)

}
