
# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)

# Source your C++ funcitons
sourceCpp("LassoInC.cpp")

# Source your LASSO functions from HW4 (make sure to move the corresponding .R file in the current project folder)
source("LassoFunctions.R")

# Do at least 2 tests for soft-thresholding function below. You are checking output agreements on at least 2 separate inputs
#################################################
# Two example inputs
a_vals <- c(-1.5, 0.2, 1.8)
lambda1 <- 0.5
lambda2 <- 1.0

r_out1 <- sapply(a_vals, soft, lambda = lambda1)
cpp_out1 <- sapply(a_vals, function(a) soft_c(a, lambda1))
stopifnot(all.equal(r_out1, cpp_out1, tolerance = 1e-10))

r_out2 <- sapply(a_vals, soft, lambda = lambda2)
cpp_out2 <- sapply(a_vals, function(a) soft_c(a, lambda2))
stopifnot(all.equal(r_out2, cpp_out2, tolerance = 1e-10))

cat("oft_c matches soft() for two test cases\n")




# Do at least 2 tests for lasso objective function below. You are checking output agreements on at least 2 separate inputs
#################################################
n <- 10; p <- 3
X <- scale(matrix(rnorm(n * p), n, p))
Y <- scale(rnorm(n))
beta <- rnorm(p)
lambda <- 0.1

fR <- lasso(X, Y, beta, lambda)
fC <- lasso_c(X, Y, beta, lambda)
stopifnot(all.equal(fR, fC, tolerance = 1e-10))

# another test
beta2 <- rep(0, p) + rnorm(p, sd = 0.5)
lambda2 <- 0.5
fR2 <- lasso(X, Y, beta2, lambda2)
fC2 <- lasso_c(X, Y, beta2, lambda2)
stopifnot(all.equal(fR2, fC2, tolerance = 1e-10))

cat("lasso_c matches lasso() for two test cases\n")



# Do at least 2 tests for fitLASSOstandardized function below. You are checking output agreements on at least 2 separate inputs
#################################################
lambda <- 0.1
fitR <- fitLASSOstandardized(X, Y, lambda)
fitC <- fitLASSOstandardized_c(X, Y, lambda, rep(0, ncol(X)))

stopifnot(all.equal(as.numeric(fitR$beta), as.numeric(fitC), tolerance = 1e-6))

# second test with different lambda
lambda2 <- 0.5
fitR2 <- fitLASSOstandardized(X, Y, lambda2)
fitC2 <- fitLASSOstandardized_c(X, Y, lambda2, rep(0, ncol(X)))
stopifnot(all.equal(as.numeric(fitR2$beta), as.numeric(fitC2), tolerance = 1e-6))

cat("fitLASSOstandardized_c matches R version for two test cases\n")


# Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################
mb1 <- microbenchmark(
  R = fitLASSOstandardized(X, Y, 0.1),
  Cpp = fitLASSOstandardized_c(X, Y, 0.1, rep(0, ncol(X))),
  times = 20
)
print(mb1)
cat("Speedup (median R / median C++) =", 
    median(mb1$time[mb1$expr=="R"]) / median(mb1$time[mb1$expr=="Cpp"]), "\n\n")


# Do at least 2 tests for fitLASSOstandardized_seq function below. You are checking output agreements on at least 2 separate inputs
#################################################
lambda_seq <- seq(0.5, 0.05, length.out = 5)

fitRseq <- fitLASSOstandardized_seq(X, Y, lambda_seq)
fitCseq <- fitLASSOstandardized_seq_c(X, Y, lambda_seq)


stopifnot(all.equal(as.numeric(fitRseq$beta_mat), 
                    as.numeric(fitCseq), tolerance = 1e-5))

# another small test
lambda_seq2 <- seq(0.3, 0.1, length.out = 4)
fitRseq2 <- fitLASSOstandardized_seq(X, Y, lambda_seq2)
fitCseq2 <- fitLASSOstandardized_seq_c(X, Y, lambda_seq2)
stopifnot(all.equal(as.numeric(fitRseq2$beta_mat), 
                    as.numeric(fitCseq2), tolerance = 1e-5))

cat("fitLASSOstandardized_seq_c matches R version for two test cases\n")


# Do microbenchmark on fitLASSOstandardized_seq vs fitLASSOstandardized_seq_c
######################################################################

mb2 <- microbenchmark(
  R = fitLASSOstandardized_seq(X, Y, lambda_seq),
  Cpp = fitLASSOstandardized_seq_c(X, Y, lambda_seq),
  times = 10
)
print(mb2)
cat("Speedup (median R / median C++) =", 
    median(mb2$time[mb2$expr=="R"]) / median(mb2$time[mb2$expr=="Cpp"]), "\n\n")

# Tests on riboflavin data
##########################
require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression

# Make sure riboflavin$x is treated as matrix later in the code for faster computations
class(riboflavin$x) <- class(riboflavin$x)[-match("AsIs", class(riboflavin$x))]

# Standardize the data
out <- standardizeXY(riboflavin$x, riboflavin$y)

# This is just to create lambda_seq, can be done faster, but this is simpler
outl <- fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, n_lambda = 30)

# The code below should assess your speed improvement on riboflavin data
microbenchmark(
  fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, outl$lambda_seq),
  fitLASSOstandardized_seq_c(out$Xtilde, out$Ytilde, outl$lambda_seq),
  times = 10
)
