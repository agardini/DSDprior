center_col <- function(mat) mat - rep(1, nrow(mat)) %*% spam::t(spam::colMeans(mat))


my_rank <- function(mat){
  res <- Matrix::rankMatrix(mat, method = "qr")
  res[[1]]
}


rMG <- function(n, alpha, beta, b, p, q){
  y <- GB2::rgb2(n = n, shape1 = 1, scale = b, shape2 = p, shape3 = q)
  x <- rgamma(n = n, shape = alpha, rate = beta/y)
  return(x)
}



#' @title Prior elicitation under B2 prior
#'
#' @description Allow to retrieve the B2 scale parameter b under the considered probabilistic statement
#'
#' @param pi probability related to the quantile.
#' @param c fixed quantile.
#' @param n data sample size.
#' @param p_DSD value of the parameter \code{p} of the DSD prior.
#' @param q_DSD value of the parameter \code{q} of the DSD prior.
#' @param n_sim number of Monte Carlo replicates to evaluate the quantile.
#'
#'
#' @return
#'
#' The parameter \code{b} that delivers a prior distribution for the random effect sample variance \eqn{V} satisfying the relation \deqn{
#' P[V>c]=pi
#' }
#'
#' @export
get_b_DSD <- function(pi_0 = 0.5, c, n, p_DSD, q_DSD, n_sim= 1e6){
  V_eta_sims <- rMG(n = n_sim, alpha = (n - 1) / 2,
                    beta = (n - 1) / 2, b = 1, p = p_DSD, q = q_DSD)
  q_star <- quantile(V_eta_sims, prob = pi_0)
  return(c / q_star)
}



#' @title Function to set an Inverse Gamma prior
#'
#' @description Allow the user to specify an Inverse Gamma prior for a variance parameter
#'
#' @param a shape parameter.
#' @param b scale parameter.
#'
#'
#' @return
#'
#' A list with the prior name and all the quantities required to implement it within the sampler.
#' @export
IG_pri <- function(a, b){
  list(Kind_pri = "InvGamma", a = a, b = b)
}

#' @title Function to set a B2 prior
#'
#' @description Allow the user to specify a B2 prior for a variance parameter
#'
#' @param b scale parameter.
#' @param p shape parameter.
#' @param q shape parameter.
#'
#'
#' @return
#'
#' A list with the prior name and all the quantities required to implement it within the sampler.
#' @export
B2_pri <- function(b, p, q){
  list(Kind_pri = "B2", b = b, p = p, q = q)
}

#' @title Function to set a Gamma prior
#'
#' @description Allow the user to specify a Gamma prior for a variance parameter
#'
#' @param a shape parameter.
#' @param b scale parameter.
#'
#'
#' @return
#'
#' A list with the prior name and all the quantities required to implement it within the sampler.
#' @export
G_pri <- function(a, b){
  list(Kind_pri = "Gamma", a = a, b = b)
}

#' @title Function to set a Weibull prior
#'
#' @description Allow the user to specify a Weibull prior for a variance parameter
#'
#' @param lambda shape parameter.
#' @param kappa scale parameter.
#'
#'
#' @return
#'
#' A list with the prior name and all the quantities required to implement it within the sampler.
#' @export
W_pri <- function(lambda, kappa){
  list(Kind_pri = "Weibull", lambda = lambda, kappa = kappa)
}


#' @title Function to set a DSD prior in a model
#'
#' @description Allow the user to specify a DSD prior for the scale parameter related to a Gaussian random effect
#'
#' @param alpha parameter related to the design.
#' @param beta parameter related to the design.
#' @param alpha_til parameter related to the precision matrix.
#' @param beta_til parameter related to the precision matrix.
#' @param b scale parameter of the underlying B2 prior.
#' @param p shape parameter of the underlying B2 prior.
#' @param q shape parameter of the underlying B2 prior.
#'
#'
#'
#' @return
#'
#' A list with the prior name and all the quantities required to implement it within the sampler.
#' @export
DSD_pri <- function(alpha, beta, alpha_til, beta_til, b, p, q){
  list(Kind_pri = "DSD", alpha = alpha,
       alpha_til = alpha_til,
       beta = beta, beta_til = beta_til,
       b = b, p = p, q = q)
}




#' @title Function to set a standard prior on beta
#'
#' @description Allow the user to specify a standard Gaussian prior for the regression coefficient beta.
#'
#'
#'
#'
#' @return
#'
#' A list with the prior name and all the quantities required to implement it within the sampler.
#' @export
Gaussian_pri_beta <- function(){
  list(Mixture = "No")
}




#' @title Function to set a DSD prior on the scale parameter of beta
#'
#' @description Allow the user to specify a DSD prior for the scale parameter for the regression coefficient beta prior.
#'
#' @param alpha parameter related to the design.
#' @param x covariate vector.
#' @param b scale parameter of the underlying B2 prior.
#' @param q shape parameter of the underlying B2 prior.
#'
#'
#' @return
#'
#' A list with the prior name and all the quantities required to implement it within the sampler.
#' @export
DSD_pri_beta <- function(alpha, x, b, q){
  pri <- B2_pri(b / var(x), alpha, q)
  pri$Mixture <- "Yes"
  pri
}



#' @title Computing useful quantities about random effects
#'
#' @description Starting from the design matrix and the precision matrix, some details about the random effect are returned.
#'
#' @param Z Design matrix.
#' @param K Precision matrix.
#' @param K_full_rank Logical. If TRUE the function assumes that \code{K} is full-rank
#'
#'
#'
#' @return
#'
#' A list with \code{Z}, \code{K}, the rank of \code{K}, the number of columns of \code{K}, the matrix A that contain a basis
#' of the null space, the vector e with the vector employed to set the linear constraint and \code{K_full_rank}.
#'
#' @export
get_reff <- function(Z, K, K_full_rank = FALSE){
  if (!K_full_rank) {
  m_j <- ncol(K)
  rank_K <- my_rank(K)

  # size
  egiendec_K <- eigen(K)
  U0 <- egiendec_K$vectors[ , (rank_K + 1):m_j]

  A <- as.matrix(spam::t(U0) %*% (spam::t(Z) %*% Z))
  e <- rep(0, nrow(A))
  }else{
    m_j <- ncol(K)
    rank_K <- m_j
    # size
    U0 <- rep(1, m_j)
    U0 <- U0 / sqrt(sum(U0^2))
    A <- as.matrix(spam::t(U0) %*% (spam::t(Z) %*% Z))
    e <- rep(0, nrow(A))
  }
  out <- list(Z = Z, K = K, rank_K = rank_K, m_j = m_j,
              A = A, e = e, K_full_rank = K_full_rank)
  class(out) <- "reff"
  out

}

#' @title Computing useful quantities for DSD prior
#'
#' @description Starting from a list produced by the function \link{get_reff}, the quantities required to implement a DSD prior are computed.
#'
#' @param reff An object of class \code{"reff"} produced by the function \link{get_reff}
#'
#'
#'
#' @return
#'
#' A list with \eqn{\tilde{\alpha}}, \eqn{\tilde{\beta}} and the mean of the non-null eigenvalues \eqn{\lambda}.
#'
#' @export
get_DSD_par <- function(reff){
 if(class(reff)!="reff"){
   stop("An object of class 'reff' must be provided.")
 }
  with(reff, {
    if (!K_full_rank) {
      n <- nrow(Z)

      egiendec_K <- eigen(K)
      Up <- egiendec_K$vectors[ , -((rank_K + 1):m_j)]

      mat <- spam::crossprod.spam(center_col(Z), Z) %*% spam::t(spam::t(Up) * 1/egiendec_K$values[-((rank_K + 1):m_j)]) %*%  spam::t(Up)
      lambdas_pri <- Re(eigen(mat)$val)[1:my_rank(mat)]

      alpha_til <- sum(lambdas_pri)^2 / (2 * sum(lambdas_pri^2))
      beta_til <- (n - 1)  * sum(lambdas_pri) / (2*sum(lambdas_pri^2))
    }else{
      n <- nrow(Z)
      mat <- spam::crossprod.spam(center_col(Z), Z) %*% spam::solve(K)
      lambdas_pri <- Re(eigen(mat)$val)

      alpha_til <- sum(lambdas_pri)^2 / (2 * sum(lambdas_pri^2))
      beta_til <- (n - 1)  * sum(lambdas_pri) / (2*sum(lambdas_pri^2))
    }
    list(alpha_til = alpha_til, beta_til = beta_til, mean_lambdas = mean(lambdas_pri))
  })

}


