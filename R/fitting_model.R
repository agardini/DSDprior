sample_blmm <- function(y, X, Z_list, K_list, rank_K_g, A_list, e_list, S_beta, pri_s2e, pri_s2b, pri_s2g, beta_init, g_init_list, S2g_init, s2e_init, FFe, FFb, FFg, niter, pr, thin, ntuning, stop_tuning) {
  .Call(`_DSDprior_sample_blmm`, y, X, Z_list, K_list, rank_K_g, A_list, e_list, S_beta, pri_s2e, pri_s2b, pri_s2g, beta_init, g_init_list, S2g_init, s2e_init, FFe, FFb, FFg, niter, pr, thin, ntuning, stop_tuning)
}
sample_bprobm <- function(y, X, Z_list, K_list, rank_K_g, A_list, e_list, S_beta, pri_s2e, pri_s2b, pri_s2g, beta_init, g_init_list, S2g_init, s2e_init, FFe, FFb, FFg, niter, pr, thin, ntuning, stop_tuning) {
  .Call(`_DSDprior_sample_bprobm`, y, X, Z_list, K_list, rank_K_g, A_list, e_list, S_beta, pri_s2e, pri_s2b, pri_s2g, beta_init, g_init_list, S2g_init, s2e_init, FFe, FFb, FFg, niter, pr, thin, ntuning, stop_tuning)
}

sample_logitm <- function(y, X, Z_list, K_list, rank_K_g, A_list, e_list, S_beta, pri_s2e, pri_s2b, pri_s2g, beta_init, g_init_list, S2g_init, s2e_init, FFe, FFb, FFg, niter, pr, thin, ntuning, stop_tuning) {
  .Call(`_DSDprior_sample_logitm`, y, X, Z_list, K_list, rank_K_g, A_list, e_list, S_beta, pri_s2e, pri_s2b, pri_s2g, beta_init, g_init_list, S2g_init, s2e_init, FFe, FFb, FFg, niter, pr, thin, ntuning, stop_tuning)
}

sample_poissonm_nob_auto <- function(y, X, Z_list, K_list, rank_K_g, A_list, e_list, S_beta, pri_s2b, pri_s2g, beta_init, g_init_list, S2g_init, FFb, FFg, n0, log_offset, m_mix_orig, v_mix_orig, w_mix_orig, m_mix_adj, v_mix_adj, w_mix_adj, check_mix_min, check_mix_max, b_gibbs_start, b_check, threshold_MH, threshold_adj, niter, pr, thin, ntuning, stop_tuning) {
  .Call(`_DSDprior_sample_poissonm_nob_auto`, y, X, Z_list, K_list, rank_K_g, A_list, e_list, S_beta, pri_s2b, pri_s2g, beta_init, g_init_list, S2g_init, FFb, FFg, n0, log_offset, m_mix_orig, v_mix_orig, w_mix_orig, m_mix_adj, v_mix_adj, w_mix_adj, check_mix_min, check_mix_max, b_gibbs_start, b_check, threshold_MH, threshold_adj, niter, pr, thin, ntuning, stop_tuning)
}


#' @title Draw posterior samples from a Bayesian linear mixed model
#'
#' @description Allow the user to specify a Bayesian linear mixed model, with flexibility in the prior choices.
#'
#' @param y Vector of observed responses.
#' @param X Matrix with fixed effects. The intercept is added by the function.
#' @param offset Vector with offset values. Default NULL, possibly useful for Poisson regression models.
#' @param pi_0 Prior elicitation parameter.
#' @param reff_list list with the design matrices of the random effects.
#' @param p_DSD Hyperparameter determining the prior behavior near by 0. Default 0.5.
#' @param q_DSD Hyperparameter determining the tail decay of the prior. Default 1.5.
#' @param model Model to be fitted on input data. Available 'Gaussian', 'Bernoulli_logit', 'Bernoulli_probit', "Poisson".
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





sample_modelDSD <- function(y, X, offset = NULL, reff_list,
                            pi_0 = 0.5, p_DSD = 0.5, q_DSD = 1.5,
                            model = c("Gaussian", "Bernoulli_logit", "Bernoulli_probit", "Poisson"),
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


print("bb")
  # Setting b_pri
  if(model == "Gaussian"){
    c <- var(y)
  }else if(model == "Bernoulli_logit"){
    c <- 1 / (mean(y) * (1 - mean(y)))
  }else if(model == "Bernoulli_probit"){
    c <- (mean(y) * (1 - mean(y))) / dnorm(qnorm(mean(y)))^2
  }else if(model == "Poisson"){
    print("aa")
    c <- 1 / mean(y)

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
  }else if(model == "Poisson"){
    # data with mixture dimension
    ord <- order(y)
    y_s <- y[ord]
    X_s <- X[ord, ]
    n_0 <- sum(y_s == 0)
    n_aux <- 2 * n - n_0
    if(ncol(X)>1){
      X_mix <- rbind(X_s, X_s[y_s>0,])
    }else{
      X_mix <- matrix(c(X_s, X_s[y_s>0]), ncol=1)
    }
    E_mix <- offset[ord]
    E_mix <- c(E_mix, E_mix[y_s>0])
    for(i in 1:length(Z_list)){
      Z_list[[i]] <- rbind(Z_list[[i]], Z_list[[i]][y_s>0,])
    }




    # setting auxiliary variables
    w_list_orig <- vector("list", n_aux)
    m_list_orig <- vector("list", n_aux)
    v_list_orig <- vector("list", n_aux)
    for(k in 1:n){
      w_list_orig[[k]] <- DSDprior::original_mixtures[[1]]$w
      m_list_orig[[k]] <- DSDprior::original_mixtures[[1]]$m
      v_list_orig[[k]] <- DSDprior::original_mixtures[[1]]$v
    }
    for(k in (n + 1):n_aux){
      index <- y_s[n_0 + k - n]
      if(index > 30000){
        index<-30000
      }
      w_list_orig[[k]] <- DSDprior::original_mixtures[[index]]$w
      m_list_orig[[k]] <- DSDprior::original_mixtures[[index]]$m
      v_list_orig[[k]] <- DSDprior::original_mixtures[[index]]$v
    }

    ################################
    w_list_adj <- vector("list", n_aux)
    m_list_adj <- vector("list", n_aux)
    v_list_adj <- vector("list", n_aux)
    for(k in 1:n){
      w_list_adj[[k]] <- DSDprior::adjusted_mixtures[[1]]$w
      m_list_adj[[k]] <- DSDprior::adjusted_mixtures[[1]]$m
      v_list_adj[[k]] <- DSDprior::adjusted_mixtures[[1]]$v
    }
    for(k in (n + 1):n_aux){
      index <- y_s[n_0 + k - n]
      if(index > 30000){
        index<-30000
      }
      w_list_adj[[k]] <- DSDprior::adjusted_mixtures[[index]]$w
      m_list_adj[[k]] <- DSDprior::adjusted_mixtures[[index]]$m
      v_list_adj[[k]] <- DSDprior::adjusted_mixtures[[index]]$v
    }
    appo <- y_s[y_s>0]
    index_upper <- ifelse(appo>30000, yes = 30000, no = appo)
    # print(index_upper)
    check_mix <- c(rep(DSDprior::eps_max[1], n), DSDprior::eps_max[index_upper])
    ind_check <- 1000
    check_mix_max <- c(rep(DSDprior::eps_max[1], n), DSDprior::eps_max[index_upper])
    check_mix_min <- c(rep(DSDprior::eps_min[1], n), DSDprior::eps_min[index_upper])

    # tuning parameters MCMC algorithm
    T1 = 500
    T2 = 250
    pL = 0.05
    pU = 0.05

    out_mcmc <- sample_poissonm_nob_auto(y = y_s, X = X_mix,
                                         Z_list = Z_list, K_list = K_list,rank_K_g = rank_K,
                                         A_list = A_list, e_list = e_list, S_beta = S_beta,
                                         pri_s2b = pri_s2b,
                                         pri_s2g = pri_s2g,
                                         beta_init = beta_init,
                                         g_init_list = g_init_list,
                                         S2g_init = s2g_init,
                                         FFg = FFg, FFb = FFb,
                                         n0 = n_0, log_offset = log(E_mix),
                                         m_mix_orig = m_list_orig, v_mix_orig = v_list_orig, w_mix_orig = w_list_orig,
                                         m_mix_adj = m_list_adj, v_mix_adj = v_list_adj, w_mix_adj = w_list_adj,
                                         check_mix_max = check_mix_max, check_mix_min = check_mix_min,
                                         b_gibbs_start = T1, b_check = T2, threshold_MH = pL, threshold_adj = pU,
                                         niter = n_iter, pr = pr, thin = thin,
                                         ntuning = ntuning, stop_tuning = stop_tuning)
  }

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
    out_mcmc$s2e <- matrix(out_mcmc$s2e, ncol = 1)
    colnames(out_mcmc$s2e) <- "s2_epsilon"
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
#' @param offset Vector with offset values. Default NULL, possibly useful for Poisson regression models.
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





sample_model_manual <- function(y, X, offset = NULL, reff_list,
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
  }else if(model == "Poisson"){
    # data with mixture dimension
    ord <- order(y)
    y_s <- y[ord]
    X_s <- X[ord, ]
    n_0 <- sum(y_s == 0)
    n_aux <- 2 * n - n_0
    if(ncol(X)>1){
      X_mix <- rbind(X_s, X_s[y_s>0,])
    }else{
      X_mix <- matrix(c(X_s, X_s[y_s>0]), ncol=1)
    }
    E_mix <- offset[ord]
    E_mix <- c(E_mix, E_mix[y_s>0])
    for(i in 1:length(Z_list)){
      Z_list[[i]] <- rbind(Z_list[[i]], Z_list[[i]][y_s>0,])
    }




    # setting auxiliary variables
    w_list_orig <- vector("list", n_aux)
    m_list_orig <- vector("list", n_aux)
    v_list_orig <- vector("list", n_aux)
    for(k in 1:n){
      w_list_orig[[k]] <- DSDprior::original_mixtures[[1]]$w
      m_list_orig[[k]] <- DSDprior::original_mixtures[[1]]$m
      v_list_orig[[k]] <- DSDprior::original_mixtures[[1]]$v
    }
    for(k in (n + 1):n_aux){
      index <- y_s[n_0 + k - n]
      if(index > 30000){
        index<-30000
      }
      w_list_orig[[k]] <- DSDprior::original_mixtures[[index]]$w
      m_list_orig[[k]] <- DSDprior::original_mixtures[[index]]$m
      v_list_orig[[k]] <- DSDprior::original_mixtures[[index]]$v
    }

    ################################
    w_list_adj <- vector("list", n_aux)
    m_list_adj <- vector("list", n_aux)
    v_list_adj <- vector("list", n_aux)
    for(k in 1:n){
      w_list_adj[[k]] <- DSDprior::adjusted_mixtures[[1]]$w
      m_list_adj[[k]] <- DSDprior::adjusted_mixtures[[1]]$m
      v_list_adj[[k]] <- DSDprior::adjusted_mixtures[[1]]$v
    }
    for(k in (n + 1):n_aux){
      index <- y_s[n_0 + k - n]
      if(index > 30000){
        index<-30000
      }
      w_list_adj[[k]] <- DSDprior::adjusted_mixtures[[index]]$w
      m_list_adj[[k]] <- DSDprior::adjusted_mixtures[[index]]$m
      v_list_adj[[k]] <- DSDprior::adjusted_mixtures[[index]]$v
    }
    appo <- y_s[y_s>0]
    index_upper <- ifelse(appo>30000, yes = 30000, no = appo)
    # print(index_upper)
    check_mix <- c(rep(DSDprior::eps_max[1], n), DSDprior::eps_max[index_upper])
    ind_check <- 1000
    check_mix_max <- c(rep(DSDprior::eps_max[1], n), DSDprior::eps_max[index_upper])
    check_mix_min <- c(rep(DSDprior::eps_min[1], n), DSDprior::eps_min[index_upper])

    # tuning parameters MCMC algorithm
    T1 = 500
    T2 = 250
    pL = 0.05
    pU = 0.05

    out_mcmc <- sample_poissonm_nob_auto(y = y_s, X = X_mix,
                                         Z_list = Z_list, K_list = K_list,rank_K_g = rank_K,
                                         A_list = A_list, e_list = e_list, S_beta = S_beta,
                                         pri_s2b = pri_s2b,
                                         pri_s2g = pri_s2g,
                                         beta_init = beta_init,
                                         g_init_list = g_init_list,
                                         S2g_init = s2g_init,
                                         FFg = FFg, FFb = FFb,
                                         n0 = n_0, log_offset = log(E_mix),
                                         m_mix_orig = m_list_orig, v_mix_orig = v_list_orig, w_mix_orig = w_list_orig,
                                         m_mix_adj = m_list_adj, v_mix_adj = v_list_adj, w_mix_adj = w_list_adj,
                                         check_mix_max = check_mix_max, check_mix_min = check_mix_min,
                                         b_gibbs_start = T1, b_check = T2, threshold_MH = pL, threshold_adj = pU,
                                         niter = n_iter, pr = pr, thin = thin,
                                         ntuning = ntuning, stop_tuning = stop_tuning)
  }






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
    out_mcmc$s2e <- matrix(out_mcmc$s2e, ncol = 1)
    colnames(out_mcmc$s2e) <- "s2_epsilon"
    out[["s2"]] <- cbind(out_mcmc[["s2e"]], out_mcmc[["s2g"]])
  }else{
    out[["s2"]] <- out_mcmc[["s2g"]]
  }

  out <- purrr::map(out, dplyr::as_tibble)


  return(out)

}

