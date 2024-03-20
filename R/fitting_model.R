sample_blmm <- function(y, X, Z_list, K_list, rank_K_g, A_list, e_list, S_beta, pri_s2e, pri_s2b, pri_s2g, beta_init, g_init_list, S2g_init, s2e_init, FFe, FFb, FFg, niter, pr, thin, ntuning, stop_tuning) {
  .Call(`_DSDprior_sample_blmm`, y, X, Z_list, K_list, rank_K_g, A_list, e_list, S_beta, pri_s2e, pri_s2b, pri_s2g, beta_init, g_init_list, S2g_init, s2e_init, FFe, FFb, FFg, niter, pr, thin, ntuning, stop_tuning)
}
sample_bprobm <- function(y, X, Z_list, K_list, rank_K_g, A_list, e_list, S_beta, pri_s2e, pri_s2b, pri_s2g, beta_init, g_init_list, S2g_init, s2e_init, FFe, FFb, FFg, niter, pr, thin, ntuning, stop_tuning) {
  .Call(`_DSDprior_sample_bprobm`, y, X, Z_list, K_list, rank_K_g, A_list, e_list, S_beta, pri_s2e, pri_s2b, pri_s2g, beta_init, g_init_list, S2g_init, s2e_init, FFe, FFb, FFg, niter, pr, thin, ntuning, stop_tuning)
}

sample_logitm <- function(y, X, Z_list, K_list, rank_K_g, A_list, e_list, S_beta, pri_s2e, pri_s2b, pri_s2g, beta_init, g_init_list, S2g_init, s2e_init, FFe, FFb, FFg, niter, pr, thin, ntuning, stop_tuning) {
  .Call(`_DSDprior_sample_logitm`, y, X, Z_list, K_list, rank_K_g, A_list, e_list, S_beta, pri_s2e, pri_s2b, pri_s2g, beta_init, g_init_list, S2g_init, s2e_init, FFe, FFb, FFg, niter, pr, thin, ntuning, stop_tuning)
}


#' @title Draw posterior samples from a Bayesian linear mixed model
#'
#' @description Allow the user to specify a Bayesian linear mixed model, with flexibility in the prior choices.
#'
#' @param y Vector of observed responses.
#' @param X Matrix with fixed effects. The intercept is added by the function.
#' @param pi_0 Prior elicitation parameter.
#' @param reff_list list with the design matrices of the random effects.
#' @param p_DSD Hyperparameter determining the prior behavior near by 0. Default 0.5.
#' @param q_DSD Hyperparameter determining the tail decay of the prior. Default 1.5.
#' @param model Model to be fitted on input data. Available 'Gaussian', 'Bernoulli_logit' and 'Bernoulli_probit'.
#' @param beta_init Initial values for the vector of regression coefficients. Default all 0s.
#' @param s2g_init Vector with the initial values for the random effects scales. Default all 1s.
#' @param s2e_init Initial value for the data variance. Default 1.
#' @param niter Number of MCMC iterations.
#' @param pr Frequency of iterations for printing the progress.
#' @param thin Number of thinned observations.
#'
#'
#' @return
#'
#' The function returns a list with the drowns from the posterior distribution of the model parameters.
#'
#'
#'
#'
#'
#' @export
#'





sample_modelDSD <- function(y, X, reff_list,
                            pi_0 = 0.5, p_DSD = 0.5, q_DSD = 1.5,
                            model = c("Gaussian", "Bernoulli_logit", "Bernoulli_probit"),
                            beta_init = rep(0, ncol(X)),
                            s2g_init = rep(1, length(reff_list)),
                            s2e_init = 1,
                            niter = 25000, pr = niter / 10, thin = 10){


  model <- match.arg(model)

  # addding intercept and then standardising
  if(!all(X[,1]==1)){
    X <- cbind(1, X)
  }
  X[,-1] <- scale(X[,-1], scale = F)[,]

  # Basic quantities
  n <- length(y)
  n_eff <- length(reff_list)
  if(is.null(names(reff_list))){
    names(reff_list) <- paste0("reff",1:n_eff)
  }


  # Tuning parameters RW proposal
  FFe = 1.5
  FFb = rep(1.5, ncol(X))
  FFg = rep(1.5, n_eff)
  ntuning = 100
  stop_tuning = 5000

  # Setting Random effects
  K_list <- purrr::map(reff_list, ~.$K)
  names(K_list)<-NULL
  Z_list <- purrr::map(reff_list, ~.$Z)
  rank_K <- purrr::map_dbl(reff_list, ~.$rank)
  names(rank_K)<-NULL
  A_list <- purrr::map(reff_list, ~.$A)
  names(A_list)<-NULL
  e_list <- purrr::map(reff_list, ~.$e)
  names(e_list)<-NULL
  pars_list <- purrr::map(reff_list, ~DSDprior::get_DSD_par(.))
  alpha_til_vector <- purrr::map_dbl(pars_list, ~.$alpha_til)
  beta_til_vector <- purrr::map_dbl(pars_list, ~.$beta_til)



  # Setting b_pri
  if(model == "Gaussian"){
    c <- var(y)
  }else if(model == "Bernoulli_logit"){
    c <- 1 / (mean(y) * (1 - mean(y)))
  }else if(model == "Bernoulli_probit"){
    c <- (mean(y) * (1 - mean(y))) / dnorm(qnorm(mean(y)))^2
  }


  b_pri <- DSDprior::get_b_DSD(pi_0 = pi_0, c = c, n = n, p_DSD = p_DSD, q_DSD = q_DSD)


  # Setting priors
  pri_s2e <-  DSDprior::W_pri(var(y), .5)
  v_beta <- 100
  S_beta <- diag(ncol(X)) * v_beta
  v_alpha <- 10000*c
  S_beta[1, 1] <- v_alpha

  # scales for random effects
  pri_s2g <- list()
  for(i in 1:n_eff){
    pri_s2g[[i]] <- DSDprior::DSD_pri((n - 1) / 2, (n - 1) / 2,
                                     alpha_til_vector[i], beta_til_vector[i],
                                     b_pri, p_DSD, q_DSD)
  }
  # on random effects
  pri_s2b <- vector("list", ncol(X))
  pri_s2b[[1]] <- list("Mixture" = "No")
  if(ncol(X)>1){
    for(k in 2:ncol(X)) {
      pri_s2b[[k]] <- DSDprior::DSD_pri_beta((n - 1) / 2, X[,k], b_pri, q_DSD)
    }
  }
  g_init_list = purrr::map2(rep(0, length(Z_list)), sapply(Z_list, ncol), rep)


  # sampling
  if(model == "Gaussian"){
    out_mcmc <- sample_blmm(y = y, X = X,
                            Z_list = Z_list, K_list = K_list,
                            rank_K_g = rank_K,
                            A_list = A_list, e_list = e_list,
                            S_beta =  S_beta, pri_s2e = pri_s2e,
                            pri_s2b = pri_s2b,
                            pri_s2g = pri_s2g,
                            beta_init = beta_init,
                            g_init_list = g_init_list,
                            S2g_init = s2g_init, s2e_init = s2e_init,
                            FFe = FFe, FFb = FFb, FFg = FFg,
                            niter = niter,pr = pr,thin = thin,
                            ntuning = ntuning, stop_tuning = stop_tuning)

  }else if(model == "Bernoulli_logit"){
    out_mcmc <- sample_logitm(y = y, X = X,
                              Z_list = Z_list, K_list = K_list,
                              rank_K_g = rank_K,
                              A_list = A_list, e_list = e_list,
                              S_beta =  S_beta, pri_s2e = pri_s2e,
                              pri_s2b = pri_s2b,
                              pri_s2g = pri_s2g,
                              beta_init = beta_init,
                              g_init_list = g_init_list,
                              S2g_init = s2g_init, s2e_init = s2e_init,
                              FFe = FFe, FFb = FFb, FFg = FFg,
                              niter = niter,pr = pr,thin = thin,
                              ntuning = ntuning,stop_tuning = stop_tuning)
  }else if(model == "Bernoulli_probit"){
    out_mcmc <- sample_bprobm(y = y, X = X,
                              Z_list = Z_list, K_list = K_list,
                              rank_K_g = rank_K,
                              A_list = A_list, e_list = e_list,
                              S_beta =  S_beta, pri_s2e = pri_s2e,
                              pri_s2b = pri_s2b,
                              pri_s2g = pri_s2g,
                              beta_init = beta_init,
                              g_init_list = g_init_list,
                              S2g_init = s2g_init, s2e_init = s2e_init,
                              FFe = FFe, FFb = FFb, FFg = FFg,
                              niter = niter,pr = pr,thin = thin,
                              ntuning = ntuning,stop_tuning = stop_tuning)
  }




  out_mcmc$s2e <- matrix(out_mcmc$s2e, ncol = 1)
  colnames(out_mcmc$s2e) <- "s2_epsilon"

   # Fixed effects
  if(is.null(colnames(X))){
    colnames(out_mcmc$beta) <- paste0("beta_", 0:(ncol(X)-1))}
  if(!is.null(colnames(X))){
    colnames(out_mcmc$beta) <- colnames(X)}

  # random effects
  if(is.null(names(Z_list))){
    names(Z_list) <- paste0("nu_",1:length(Z_list))
  }
  names(out_mcmc$gamma) <- names(Z_list)

  # reff's scalers
  colnames(out_mcmc$s2g) <- paste0("s2_", names(Z_list))

  out <- list()
  out[["beta"]] <- out_mcmc[["beta"]]
  for(k in 1:length(Z_list)){
    out[[names(Z_list)[k]]] <- out_mcmc$gamma[[names(Z_list)[k]]]
    colnames(out[[names(Z_list)[k]]]) <- paste0(names(Z_list)[k], "_", 1:ncol(Z_list[[k]]))
  }
  if(model == "Gaussian"){
    out[["s2"]] <- cbind(out_mcmc[["s2e"]], out_mcmc[["s2g"]])
  }else{
    out[["s2"]] <- out_mcmc[["s2g"]]
  }

  out <- purrr::map(out, dplyr::as_tibble)


  return(out)

}



#' @title Draw posterior samples from a Bayesian linear mixed model
#'
#' @description Allow the user to specify a Bayesian linear mixed model, with flexibility in the prior choices.
#'
#' @param y Vector of observed responses.
#' @param X Matrix with fixed effects. The intercept is added by the function.
#' @param reff_list list with the design matrices of the random effects.
#' @param pri_s2e Prior for the variance parameter in case of Gaussian model.
#' @param pri_s2g List of priors for the scale hyperparameters ruling the random effects.
#' @param pri_s2b List of priors for the regression coefficients related to fixed effects.
#' @param model Model to be fitted on input data. Available 'Gaussian', 'Bernoulli_logit' and 'Bernoulli_probit'.
#' @param beta_init Initial values for the vector of regression coefficients. Default all 0s.
#' @param s2g_init Vector with the initial values for the random effects scales. Default all 1s.
#' @param s2e_init Initial value for the data variance. Default 1.
#' @param niter Number of MCMC iterations.
#' @param pr Frequency of iterations for printing the progress.
#' @param thin Number of thinned observations.
#'
#' @return
#'
#' The function returns a list with the drowns from the posterior distribution of the model parameters.
#'
#'
#'
#'
#'
#'
#' @export
#'





sample_model_manual <- function(y, X, reff_list,
                            model = c("Gaussian", "Bernoulli_logit", "Bernoulli_probit"),
                            pri_s2e = NULL, pri_s2g, pri_s2b,
                            beta_init = rep(0, ncol(X)),
                            s2g_init = rep(1, length(reff_list)),
                            s2e_init = 1,
                            niter = 25000, pr = niter / 10, thin = 10){


  model <- match.arg(model)

  # addding intercept and then standardising
  if(!all(X[,1]==1)){
    X <- cbind(1, X)
  }
  X[,-1] <- scale(X[,-1], scale = F)[,]

  # Basic quantities
  n <- length(y)
  n_eff <- length(reff_list)
  if(is.null(names(reff_list))){
    names(reff_list) <- paste0("reff",1:n_eff)
  }


  # Tuning parameters RW proposal
  FFe = 1.5
  FFb = rep(1.5, ncol(X))
  FFg = rep(1.5, n_eff)
  ntuning = 100
  stop_tuning = 5000

  # Setting Random effects
  K_list <- purrr::map(reff_list, ~.$K)
  names(K_list)<-NULL
  Z_list <- purrr::map(reff_list, ~.$Z)
  rank_K <- purrr::map_dbl(reff_list, ~.$rank)
  names(rank_K)<-NULL
  A_list <- purrr::map(reff_list, ~.$A)
  names(A_list)<-NULL
  e_list <- purrr::map(reff_list, ~.$e)
  names(e_list)<-NULL



  # Setting priors
  v_beta <- 100
  S_beta <- diag(ncol(X)) * v_beta
  v_alpha <- 10000
  S_beta[1, 1] <- v_alpha

  g_init_list = purrr::map2(rep(0, length(Z_list)), sapply(Z_list, ncol), rep)

  if (is.null(pri_s2e)){
    pri_s2e <-  DSDprior::W_pri(var(y), .5)
  }


  # sampling
  if(model == "Gaussian"){
    out_mcmc <- sample_blmm(y = y, X = X,
                            Z_list = Z_list, K_list = K_list,
                            rank_K_g = rank_K,
                            A_list = A_list, e_list = e_list,
                            S_beta =  S_beta, pri_s2e = pri_s2e,
                            pri_s2b = pri_s2b,
                            pri_s2g = pri_s2g,
                            beta_init = beta_init,
                            g_init_list = g_init_list,
                            S2g_init = s2g_init, s2e_init = s2e_init,
                            FFe = FFe, FFb = FFb, FFg = FFg,
                            niter = niter,pr = pr,thin = thin,
                            ntuning = ntuning, stop_tuning = stop_tuning)

  }else if(model == "Bernoulli_logit"){
    out_mcmc <- sample_logitm(y = y, X = X,
                              Z_list = Z_list, K_list = K_list,
                              rank_K_g = rank_K,
                              A_list = A_list, e_list = e_list,
                              S_beta =  S_beta, pri_s2e = pri_s2e,
                              pri_s2b = pri_s2b,
                              pri_s2g = pri_s2g,
                              beta_init = beta_init,
                              g_init_list = g_init_list,
                              S2g_init = s2g_init, s2e_init = s2e_init,
                              FFe = FFe, FFb = FFb, FFg = FFg,
                              niter = niter,pr = pr,thin = thin,
                              ntuning = ntuning,stop_tuning = stop_tuning)
  }else if(model == "Bernoulli_probit"){
    out_mcmc <- sample_bprobm(y = y, X = X,
                              Z_list = Z_list, K_list = K_list,
                              rank_K_g = rank_K,
                              A_list = A_list, e_list = e_list,
                              S_beta =  S_beta, pri_s2e = pri_s2e,
                              pri_s2b = pri_s2b,
                              pri_s2g = pri_s2g,
                              beta_init = beta_init,
                              g_init_list = g_init_list,
                              S2g_init = s2g_init, s2e_init = s2e_init,
                              FFe = FFe, FFb = FFb, FFg = FFg,
                              niter = niter,pr = pr,thin = thin,
                              ntuning = ntuning,stop_tuning = stop_tuning)
  }




  out_mcmc$s2e <- matrix(out_mcmc$s2e, ncol = 1)
  colnames(out_mcmc$s2e) <- "s2_epsilon"

  # Fixed effects
  if(is.null(colnames(X))){
    colnames(out_mcmc$beta) <- paste0("beta_", 0:(ncol(X)-1))}
  if(!is.null(colnames(X))){
    colnames(out_mcmc$beta) <- colnames(X)}

  # random effects
  if(is.null(names(Z_list))){
    names(Z_list) <- paste0("nu_",1:length(Z_list))
  }
  names(out_mcmc$gamma) <- names(Z_list)

  # reff's scalers
  colnames(out_mcmc$s2g) <- paste0("s2_", names(Z_list))

  out <- list()
  out[["beta"]] <- out_mcmc[["beta"]]
  for(k in 1:length(Z_list)){
    out[[names(Z_list)[k]]] <- out_mcmc$gamma[[names(Z_list)[k]]]
    colnames(out[[names(Z_list)[k]]]) <- paste0(names(Z_list)[k], "_", 1:ncol(Z_list[[k]]))
  }
  if(model == "Gaussian"){
    out[["s2"]] <- cbind(out_mcmc[["s2e"]], out_mcmc[["s2g"]])
  }else{
    out[["s2"]] <- out_mcmc[["s2g"]]
  }

  out <- purrr::map(out, dplyr::as_tibble)


  return(out)

}

