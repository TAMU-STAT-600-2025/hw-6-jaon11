# [ToDo] Standardize X and Y: center both X and Y; scale centered X
# X - n x p matrix of covariates
# Y - n x 1 response vector
standardizeXY <- function(X, Y){
  X <- as.matrix(X); storage.mode(X) <- "double"
  Y <- as.numeric(Y)
  
  n <- nrow(X); p <- ncol(X)
  # [ToDo] Center Y
  Ymean <- mean(Y)
  Ytilde <- Y - Ymean
  
  
  # [ToDo] Center and scale X
  Xmeans <- colMeans(X)
  Xc     <- sweep(X, 2, Xmeans, FUN = "-")
  
  weights <- sqrt(colSums(Xc^2) / n)                 # length p
  weights[weights == 0] <- 1
  Xtilde  <- sweep(Xc, 2, weights, FUN = "/")
  
  
  
  
  # Return:
  # Xtilde - centered and appropriately scaled X
  # Ytilde - centered Y
  # Ymean - the mean of original Y
  # Xmeans - means of columns of X (vector)
  # weights - defined as sqrt(X_j^{\top}X_j/n) after centering of X but before scaling
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

# [ToDo] Soft-thresholding of a scalar a at level lambda 
# [OK to have vector version as long as works correctly on scalar; will only test on scalars]
soft <- function(a, lambda){
  return(sign(a) * pmax(abs(a) - lambda, 0))
}

# [ToDo] Calculate objective function of lasso given current values of Xtilde, Ytilde, beta and lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta - value of beta at which to evaluate the function
lasso <- function(Xtilde, Ytilde, beta, lambda){
  n = nrow(Xtilde)
  r <- Ytilde - Xtilde %*% beta
  fval = (0.5 / n) * sum(r * r) + lambda * sum(abs(beta))

  return(fval)
  
 
}

# [ToDo] Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1 (vector)
# lamdba - tuning parameter
# beta_start - p vector, an optional starting point for coordinate-descent algorithm
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.001){
  #[ToDo]  Check that n is the same between Xtilde and Ytilde
  n = nrow(Xtilde)
  p = ncol(Xtilde)
  if (n != length(Ytilde)) {
    stop("Number of rows in Xtilde must match length of Ytilde")
  }
  
  
  #[ToDo]  Check that lambda is non-negative
  if (lambda < 0) {
    stop("Lambda must be non-negative")
  }
  
  #[ToDo]  Check for starting point beta_start. 
  # If none supplied, initialize with a vector of zeros.
  # If supplied, check for compatibility with Xtilde in terms of p
  if (is.null(beta_start)) {
    beta = rep(0, p)
  } else {
    if (length(beta_start) != p) {
      stop("Length of beta_start must match number of columns in Xtilde")
    }
    beta = beta_start
  }
  
  
  #[ToDo]  Coordinate-descent implementation. 
  # Stop when the difference between objective functions is less than eps for the first time.
  # For example, if you have 3 iterations with objectives 3, 1, 0.99999,
  # your should return fmin = 0.99999, and not have another iteration
  # Fast coordinate-descent loop that uses lasso() for fnew/fmin
  # Assumes: Xtilde (n x p), Ytilde (n), lambda, beta (length p), eps
  
  col_sq_norm_over_n <- colSums(Xtilde^2) / n
  # Maintain residual: r = y - X beta  (avoid %*% inside the inner loop)
  r <- as.numeric(Ytilde - Xtilde %*% beta)
  fmin_vec <- numeric() # to store objective values for debugging if needed
  fmin <- lasso(Xtilde, Ytilde, beta, lambda)
  fmin_vec <- c(fmin_vec, fmin)
  repeat {
    for (j in 1:p) {
      xj     <- Xtilde[, j]
      bj_old <- beta[j]
      
      # zj = (1/n) * x_j^T (y - X_{-j} beta_{-j}) = (1/n) * x_j^T (r + x_j * b_j)
      zj <- as.numeric(crossprod(xj, r + xj * bj_old)) / n
      
      denom <- col_sq_norm_over_n[j]
      bj_new <- if (denom > 0) soft(zj, lambda) / denom else 0
      
      if (bj_new != bj_old) {
        db <- bj_new - bj_old
        r  <- r - xj * db          # update residual cheaply
        beta[j] <- bj_new
      }
    }
    
    # Use your lasso() to compute fnew (one call per sweep)
    fnew <- lasso(Xtilde, Ytilde, beta, lambda)
    
    # Stop the first time |fmin - fnew| < eps
    if (abs(fmin - fnew) < eps) {
      fmin <- fnew
      break
    }
    fmin <- fnew
    fmin_vec <- c(fmin_vec, fmin) # for debugging if needed
  }
  
  
  
  # Return 
  # beta - the solution (a vector)
  # fmin - optimal function value (value of objective at beta, scalar)
  # fmin_vec - vector of objective function values at each iteration (for debugging if needed)
  return(list(beta = beta, fmin = fmin, fmin_vec = fmin_vec))
}

# [ToDo] Fit LASSO on standardized data for a sequence of lambda values. Sequential version of a previous function.
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence,
#             is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSOstandardized_seq <- function(Xtilde, Ytilde, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Check that n is the same between Xtilde and Ytilde
  n = nrow(Xtilde)
  p = ncol(Xtilde)
  if (n != length(Ytilde)) {
    stop("Number of rows in Xtilde must match length of Ytilde")
  }
  
 
  # [ToDo] Check for the user-supplied lambda-seq (see below)
  # If lambda_seq is supplied, only keep values that are >= 0,
  # and make sure the values are sorted from largest to smallest.
  # If none of the supplied values satisfy the requirement,
  # print the warning message and proceed as if the values were not supplied.
  if (!is.null(lambda_seq)) {
    lambda_seq = lambda_seq[lambda_seq >= 0]
    n_lambda = length(lambda_seq)
    if (length(lambda_seq) == 0) {
      warning("No non-negative lambda values supplied. Proceeding to calculate lambda_max.")
      lambda_seq = NULL
    } else {
      lambda_seq = sort(lambda_seq, decreasing = TRUE)
    }
  } 
  
  
  # If lambda_seq is not supplied, calculate lambda_max 
  # (the minimal value of lambda that gives zero solution),
  # and create a sequence of length n_lambda as
  if (is.null(lambda_seq)) {
    lambda_max <- max(abs(t(Xtilde) %*% Ytilde)) / nrow(Xtilde)
    lambda_seq = exp(seq(log(lambda_max), log(0.01), length = n_lambda))
  }
  
  
  # [ToDo] Apply fitLASSOstandardized going from largest to smallest lambda 
  # (make sure supplied eps is carried over). 
  # Use warm starts strategy discussed in class for setting the starting values.
  beta_mat = matrix(0, nrow = p, ncol = length(lambda_seq))
  fmin_vec = numeric(length(lambda_seq))
  beta_start = NULL
  for (i in 1:length(lambda_seq)) {
    result = fitLASSOstandardized(Xtilde, Ytilde, lambda_seq[i], beta_start, eps)
    beta_mat[, i] = result$beta
    fmin_vec[i] = result$fmin
    beta_start = result$beta # warm start
  }
  
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value
  # fmin_vec - length(lambda_seq) vector of corresponding objective function values at solution
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, fmin_vec = fmin_vec))
}

# [ToDo] Fit LASSO on original data using a sequence of lambda values
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.001
fitLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, eps = 0.001){
  # [ToDo] Center and standardize X,Y based on standardizeXY function
  standardization = standardizeXY(X, Y)
  Xtilde = standardization$Xtilde
  Ytilde = standardization$Ytilde
  
 
  # [ToDo] Fit Lasso on a sequence of values using fitLASSOstandardized_seq
  # (make sure the parameters carry over)
  lasso_result = fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq, n_lambda, eps)

  
 
  # [ToDo] Perform back scaling and centering to get original intercept and coefficient vector
  # for each lambda
  lambda_seq = lasso_result$lambda_seq
  beta_mat_tilde = lasso_result$beta_mat
  p = ncol(X)
  beta_mat = matrix(0, nrow = p, ncol = length(lambda_seq))
  beta0_vec = numeric(length(lambda_seq))
  for (i in 1:length(lambda_seq)) {
    beta_mat[, i] = beta_mat_tilde[, i] / standardization$weights
    beta0_vec[i] = standardization$Ymean - sum(standardization$Xmeans * beta_mat[, i])
  }

  
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec))
}


# [ToDo] Fit LASSO and perform cross-validation to select the best fit
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lambda_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# k - number of folds for k-fold cross-validation, default is 5
# fold_ids - (optional) vector of length n specifying the folds assignment (from 1 to max(folds_ids)), if supplied the value of k is ignored 
# eps - precision level for convergence assessment, default 0.001
cvLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 60, k = 5, fold_ids = NULL, eps = 0.001, seed = NULL){
  # [ToDo] Fit Lasso on original data using fitLASSO
  lasso_full = fitLASSO(X, Y, lambda_seq, n_lambda, eps)
  lambda_seq = lasso_full$lambda_seq
  beta_mat = lasso_full$beta_mat
  beta0_vec = lasso_full$beta0_vec
  
 
  # [ToDo] If fold_ids is NULL, split the data randomly into k folds.
  # If fold_ids is not NULL, split the data according to supplied fold_ids.
  n = nrow(X)
  p = ncol(X)
  if (is.null(fold_ids)) {
    # For reproducibility (if needed)
    if (!is.null(seed)) {
      set.seed(seed)
    }
    fold_ids <- sample(rep(1:k, length.out = n))
  } else {
    if (length(fold_ids) != n) {
      stop("Length of fold_ids must match number of rows in X")
    }
    k = max(fold_ids)
  }
  
  
  # [ToDo] Calculate LASSO on each fold using fitLASSO,
  # and perform any additional calculations needed for CV(lambda) and SE_CV(lambda)
  m <- length(lasso_full$lambda_seq)
  fold_losses <- matrix(NA_real_, nrow = k, ncol = m)  # store MSE per fold × lambda
  
  for (fold in 1:k) {
    test_indices = which(fold_ids == fold)
    train_indices = setdiff(1:n, test_indices)
    X_train = X[train_indices, , drop = FALSE]
    Y_train = Y[train_indices]
    X_test = X[test_indices, , drop = FALSE]
    Y_test = Y[test_indices]
    
    lasso_fold = fitLASSO(X_train, Y_train, lasso_full$lambda_seq, n_lambda, eps)
    
    for (j in 1:m) {
      beta_fold = lasso_fold$beta_mat[, j]
      beta0_fold = lasso_fold$beta0_vec[j]
      predictions = X_test %*% beta_fold + beta0_fold
      mse <- mean((Y_test - predictions)^2)
      fold_losses[fold, j] <- mse
    }
  }
  cvm  <- colMeans(fold_losses)
  cvse <- apply(fold_losses, 2, function(z) sd(z) / sqrt(k))
  
  
  # [ToDo] Find lambda_min
  lambda_min = lasso_full$lambda_seq[which.min(cvm)]
  

  # [ToDo] Find lambda_1SE
  min_index = which.min(cvm)
  lambda_1se = max(lasso_full$lambda_seq[cvm <= cvm[min_index] + cvse[min_index]])
  
  
  # Return output
  # Output from fitLASSO on the whole data
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  # fold_ids - used splitting into folds from 1 to k (either as supplied or as generated in the beginning)
  # lambda_min - selected lambda based on minimal rule
  # lambda_1se - selected lambda based on 1SE rule
  # cvm - values of CV(lambda) for each lambda
  # cvse - values of SE_CV(lambda) for each lambda
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec, fold_ids = fold_ids, lambda_min = lambda_min, lambda_1se = lambda_1se, cvm = cvm, cvse = cvse))
}

