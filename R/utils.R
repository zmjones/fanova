#' (Internal) Build a centered 1D spline basis
#' @importFrom splines bs
build_basis_1d <- function(x, df, degree, w) {
  B_raw <- splines::bs(x, df = df, degree = degree, intercept = TRUE)
  # Orthogonalize wrt intercept (center)
  B_centered <- scale(B_raw, center = colSums(w * B_raw), scale = FALSE)
  # Drop any constant columns that result
  B_centered[, colSums(abs(B_centered)) > 1e-10, drop = FALSE]
}

#' (Internal) Tensor product for two basis matrices
tensor2 <- function(B1, B2) {
  if (nrow(B1) != nrow(B2)) stop("Basis matrices must have same number of rows")
  if (ncol(B1) == 0 || ncol(B2) == 0) return(matrix(0, nrow = nrow(B1), ncol = 0))
  
  # Efficient row-wise Kronecker product
  m1 <- ncol(B1)
  m2 <- ncol(B2)
  B1_rep <- B1[, rep(1:m1, each = m2), drop = FALSE]
  B2_rep <- B2[, rep(1:m2, times = m1), drop = FALSE]
  
  B1_rep * B2_rep
}

#' (Internal) Weighted Gram-Schmidt residualization
w_residualize <- function(Y, X, w) {
  if (ncol(X) == 0) return(Y)
  # Weighted QR decomposition of X
  qrx <- qr(sqrt(w) * X)
  Q <- qr.Q(qrx)
  # Project Y onto the orthogonal complement of the space spanned by X
  Y_w <- sqrt(w) * Y
  Y_proj <- Q %*% (t(Q) %*% Y_w)
  (Y_w - Y_proj) / sqrt(w)
}

#' (Internal) Stable weighted least squares solver
#' @importFrom MASS ginv
wls_solve <- function(X, y, w) {
  MASS::ginv(crossprod(X, w * X)) %*% crossprod(X, w * y)
}