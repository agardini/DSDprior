#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL



#' Burn-in and thinning of posterior drawns
#'
#'@description Providing as input an object created by \code{sample_modelDSD} and \code{sample_model_manual}, the function performs burn-in and thinning of the chain.
#'
#'@param out Object created by \code{sample_modelDSD} or \code{sample_model_manual}.
#'@param n_burn Length burn-in period.
#'@param n_thin Length thinning interval.
#'
#'
#'
#' @export
burn_thin <- function(out, n_burn, n_thin = 1){
  n_it <- nrow(out[[1]])
  sel_thin <- seq((n_burn + 1), n_it, by = n_thin)
  purrr::map(out, dplyr::slice, sel_thin)
}

#'
#' Computation of posterior summaries
#'
#'@description Providing as input an object created by \code{sample_modelDSD}, \code{sample_model_manual} or \code{burn_thin}, the function computes posterior summaries.
#'
#'@param model_out Object created by \code{sample_modelDSD}, \code{sample_model_manual} or \code{burn_thin}.
#'
#'
#' @export
posterior_summaries <- function(model_out){
  model_out %>%
    purrr::map(coda::as.mcmc) %>%
    purrr::map(MCMCvis::MCMCsummary, Rhat = FALSE) %>%
    purrr::map(dplyr::as_tibble, rownames = "parameter")
}

#'
#'Evaluate linear predictor
#'
#'@description Receiving as input the result of \code{sample_modelDSD}, \code{sample_model_manual} or \code{burn_thin} functions, computes the linear predictor in the points provided through the design matrices.
#'
#'@param out Object created by \code{sample_modelDSD}, \code{sample_model_manual} or \code{burn_thin}.
#'@param X Fixed effects design matrix.
#'@param Z_list List of design matrices related to random effects.
#'
#'
#'
#'@export

compute_linpred <- function(out, X, Z_list){
  n <- nrow(X)

  # random effects
  if (is.null(names(Z_list))) {
    names(Z_list) <- paste0("nu_",1:length(Z_list))
  }

  lp <- as.matrix(out[["beta"]]) %*% t(X)
  for (k in 1:length(Z_list)) {
    lp <- lp + Matrix::t(Z_list[[k]] %*% t(as.matrix(out[[names(Z_list)[k]]])))
  }
  lp <- as.matrix(lp)
  colnames(lp) <- paste0("yhat_",1:n)
  tibble::as_tibble(lp)

}

