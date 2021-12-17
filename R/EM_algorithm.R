#' Actual EM algorithm for CENTIPEDE
#' @param k number of mixed distributions
#' @param l number of rows in X
#' @param s number of columns in X
#' @param X observed experimental data matrix
#' @param R number of reads in each motif, sum of each row in X
#' @param G prior information matrix (PWM+Cons.Score+TSS)
#' @param Z binary variable used to initialize beta's
#' @param pi_mat prior probability matrix
#' @param B beta's for the logistic regression used to calculate pi_l
#' @param a0 parameters for negative binomial distribution unbound
#' @param a1 parameters for negative binomial distribution bound
#' @param tau0 parameters for negative binomial distribution unbound
#' @param tau1 parameters for negative binomial distribution bound
#' @param lambda parameters for the multinomial distribution
#' @param prob_mat posterior probability matrix
#' @param log_prob_vec log of posterior probability of bound
#' @param max_log maximum of log posterior bound probability
#' @param loglik log likelihood
#' @param tol threshold for testing convergence
#' @param logittau0 logit scale of tau0
#' @param logittau1 logit scale of tau1
#' @param logitbeta logit scale of beta
#' @return A list of convergence, updated parameters, posterior probabilities, iterations, log likelihood list, bound or unbound status
#' #'@export



CentipedeMixEm <- setRefClass("CentipedeMixEm",
                              fields=list(k = "integer",  ## number of mixed distributions
                                          l = "integer",  ## number of rows in X
                                          s = "integer",  ## number of columns in X
                                          X = "matrix",  ## observed experimental data matrix
                                          R = "vector",  ## number of reads in each motif, sum of each row in X
                                          G = "matrix",  ## prior information matrix (PWM+Cons.Score+TSS)
                                          Z = "vector",  ## binary variable used to initialize beta's
                                          pi_mat = "matrix",  ## pi_l = Pr(Z_l=1) = pi_vec[,1]; 1-pi_l = Pr(Z_l=0) = pi_vec[,2]
                                          B = "vector",  ## beta's for the logistic regression used to calculate pi_l
                                          a0 = "numeric",  ## parameters for the negative binomial distribution
                                          a1 = "numeric",
                                          tau0 = "numeric",
                                          tau1 = "numeric",
                                          lambda = "vector",  ## parameters for the multinomial distribution
                                          prob_mat = "matrix",  ## l*k matrix, prob_mat[,1]=posterior probability of Z_l=1|X,R,G (first column)
                                          ## prob_mat[,2]=posterior probability of Z_l=0|X,R,G (second column)
                                          log_prob_vec = "vector", ##for more stable storage of probability matrix
                                          max_log = "numeric", #save current scaling factor
                                          loglik = "numeric",
                                          tol = "numeric",
                                          logittau0 = "numeric",
                                          logittau1 = "numeric",
                                          logit_beta = "vector"
                              ))


#' method to initialize CentipedeMixEm class (called with CentipedeMixEm$new)
#' @param input_data - data to initailize
#' @param num_components - number of components (k)
CentipedeMixEm$methods(initialize = function(experiment_dat, num_components, prior_info, init_binary){
  X <<- experiment_dat ## Use <<- to assign fields
  k <<- num_components
  l <<- nrow(X)
  s <<- ncol(X)
  G <<- prior_info
  R <<- rowSums(X)
  Z <<- init_binary
})

#' method to initialize parameters for E-M
CentipedeMixEm$methods(init.paras = function(){
  ## set up tol for checking convergence
  tol <<- 1e-100

  # initialize beta in a logit scale
  flogit <- function(logit_beta){
    -sum(Z*(cbind(1,G)%*%logit_beta) - log(1+exp(cbind(1,G)%*%logit_beta)))
  }

  glogit <- function(logit_beta){
    -t(Z - plogis(cbind(1,G) %*% logit_beta)) %*% cbind(1,G)
  }

  logit_beta <<- optim(rep(0,ncol(G)+1), flogit, glogit, method = "BFGS")$par

  pi_mat <<- matrix(0, nrow = l, ncol = k)
  ## stable sigmoid function
  sigmoid <- function(input){
    return(ifelse(input>=0, 1/(1+exp(-input)), exp(input)/(1+exp(input))))
  }
  pi_mat[,1] <<- plogis(cbind(1,G) %*% logit_beta)
  pi_mat[,2] <<- 1 - pi_mat[,1] ##define priors later
  print(c("initial pi:",range(pi_mat)))
  ## initialize lambda:

  ## initialize lambda
  lambda <<- as.vector(t(t(Z) %*% X))
  lambda <<- lambda/sum(lambda)
  print(c("checking initial lambda:", sum(lambda == 0)))


  tau0 <<- 0.004895833

  tau1 <<- 0.002229727
  logittau0 <<- qlogis(tau0)
  logittau1 <<- qlogis(tau1)

  a0 <<- mean(R)*tau0/(1-tau0)  ## initial a0 based on mean of negbinom
  a1 <<- mean(R)*tau1/(1-tau1)  ## initial a0 based on mean of negbinom

  prob_mat <<- matrix(NA, nrow = l, ncol = k)
  #prob_mat[,1] <<- plogis(cbind(1,G) %*% B)
  #prob_mat[,2] <<- 1 - prob_mat[,1]
  log_prob_vec <<- rep(NA,l)
  likelihood <- -Inf
  loglik <<- -Inf
  max_log <<- -Inf
})


#' E-step for E-M algorithm
CentipedeMixEm$methods(update.prob = function(){
  ##Update pl which is the prob_mat
  #function to prevent extreme values
  clipExtremes <- function(x,tol){
    Q <- quantile(x,c(tol,1-tol))
    x[x<Q[1]] <- Q[1]
    x[x>Q[2]] <- Q[2]
    x
  }

  #Here we calculate the posterior odds instead of the inverse of the posterior odds, all ratios should be flipped with respect to previous implementation

  prior_log_pi_ratio <- cbind(1,G)%*%logit_beta #calculate the log(pi/(1-pi)) directly from priors and beta
  print(c("prior_log_pi_ratio",range(prior_log_pi_ratio)))

  negbinum <- (-1)*(R+a0) * log(1+exp(logittau0)) + a0 * logittau0 + lgamma(R+a0) - lgamma(a0)

  negbidenom <- (-1)*(R+a1) * log(1+exp(logittau1)) + a1 * logittau1 + lgamma(R+a1) - lgamma(a1)

  negbi <- negbidenom - negbinum

  lambda_log_1 <- X%*%(log(s*lambda))
  lambda_log_1 <- clipExtremes(lambda_log_1,0.0001) #make more stable
  lambda_log <- lambda_log_1

  ## stable sigmoid function
  sigmoid <- function(input){
    return(ifelse(input>=0, 1/(1+exp(-input)), exp(input)/(1+exp(input))))
  }

  stable_log_sigmoid <- function(input){
    ifelse(input>0, -log(1+exp(-input)), input-log(exp(input)+1))
  }

  ## calculate log odds ratio
  log_or <- as.vector(prior_log_pi_ratio + negbi + lambda_log)
  log_prob_vec <<- stable_log_sigmoid(log_or)
  max_log <<- max(log_prob_vec)
  log_prob_vec <<- log_prob_vec + max_log


  prob_mat[,1] <<- plogis(log_or) #stable calculation from log odds
  prob_mat[,2] <<- 1 - prob_mat[,1] ## 11/28

  #complete and log-likelihood function, this is new
  stableloglike0_1 <- lgamma(R + a0) - lgamma(a0) + a0*logittau0 - (R + a0)*log(1+exp(logittau0)) - R*log(s)

  stableloglike0_1 <- clipExtremes(stableloglike0_1,0.0001)

  stableloglike0_1 <- prob_mat[,2]*stableloglike0_1

  ## add bound terms
  stableloglike1_1 <- X %*% log(lambda)
  stableloglike1_2 <- lgamma(R + a1) - lgamma(a1) + a1*logittau1 - (R + a1)*log(1+exp(logittau1))
  stableloglike1 <- stableloglike1_1 + stableloglike1_2

  stableloglike1 <- clipExtremes(stableloglike1, 0.001)
  stableloglike1 <- prob_mat[,1]*stableloglike1

  loglik0 <- sum(stableloglike0_1) + sum(stableloglike1) - sum(log(1+exp(cbind(1,G)%*%logit_beta)))

  Aux <- log(1+exp(log_or))
  Aux[log_or >18] <- log_or[log_or>18] #to make more stable for large log_or
  loglik <<- sum(Aux) + loglik0 #complete log likelihood
})

#' M-step for E-M algorithm to update B
CentipedeMixEm$methods(update.B= function(){
  stable_log <- function(input){
    stable_log = ifelse(input >= 0, input+log(1+exp(-input)), log(1+exp(input)))
    return(stable_log)
  }

  ## update beta on a logit scale
  flogit <- function(logit_beta){
    -sum(prob_mat[,1]*(cbind(1,G)%*%logit_beta) - log(1+exp(cbind(1,G)%*%logit_beta)))
  }

  glogit <- function(logit_beta){
    -t(prob_mat[,1] - plogis(cbind(1,G) %*% logit_beta)) %*% cbind(1,G)
  }

  logit_beta <<- optim(rep(0,ncol(G)+1), flogit, glogit, method = "BFGS")$par


  object_fun <- function(B){
    -mean(exp( log_prob_vec - max_log)*(cbind(1,G)%*%B)-stable_log(cbind(1,G)%*%B))
  }

  ## define a numerically stable sigmoid function
  sigmoid <- function(input){
    return(ifelse(input>=0, 1/(1+exp(-input)), exp(input)/(1+exp(input))))
  }

  gradient1 <- function(B){
    -colMeans(as.vector(exp( log_prob_vec - max_log)-sigmoid(cbind(1,G) %*% B))*cbind(1,G))
  }
  gradient <- function(B){
    -t(exp( log_prob_vec - max_log)-plogis(cbind(1,G)%*%B))%*%cbind(1,G)
  }
})

#' M-step for E-M algorithm to update pi_vec
CentipedeMixEm$methods(update.pi = function(){
  sigmoid <- function(input){
    return(ifelse(input>=0, 1/(1+exp(-input)), exp(input)/(1+exp(input))))
  }
  pi_mat[,1] <<- plogis(cbind(1,G)%*%logit_beta)
  pi_mat[,2] <<- 1 - pi_mat[,1]
})


#' M-step for E-M algorithm to update a0
CentipedeMixEm$methods(update.a0 = function(){
  object_a0 <- function(a0tau0){
    a0 <<- exp(a0tau0[1])
    logittau0 <<- a0tau0[2]
    if((a0<1E200) && (a0>1E-200)){
      LogLik <- -(mean(prob_mat[,2] * ( lgamma(R+a0) - lgamma(a0) + a0*logittau0 - (R+a0)*log(1+exp(logittau0)) ) ))
    }else{
      LogLik <- 1E300
    }
    LogLik
  }

  grad_a0 <- function(a0tau0){
    a0 <<- exp(a0tau0[1])
    logittau0 <<- a0tau0[2]
    tau0 <<- plogis(logittau0)
    gloga0 = -a0*mean(prob_mat[,2]*(digamma(a0+R)-digamma(a0)+logittau0 -log(1+exp(logittau0))))

    glogtau0 = -mean(prob_mat[,2]*(a0 - (a0+R)*(tau0)))
    return(c(gloga0,glogtau0))
  }

  loga0tau0 <- c(log(a0),logittau0)
  optim_out <- optim(loga0tau0, object_a0, grad_a0,method = "L-BFGS-B")
  a0 <<- exp(optim_out$par[1]) # assume log
  logittau0 <<- optim_out$par[2] #assume log
})


#' M-step for E-M algorithm to update a1
CentipedeMixEm$methods(update.a1 = function(){

  object_a1 <- function(a1tau1){
    a1 <<- exp(a1tau1[1])
    logittau1 <<- a1tau1[2]
    if((a1<1E200) && (a1>1E-200)){
      LogLik <- -(mean(prob_mat[,1]* ( lgamma(R+a1) - lgamma(a1) + a1*logittau1 - (R+a1)*log(1+exp(logittau1)) ) ))
    }else{
      LogLik <- 1E300
    }
    LogLik
  }

  grad_a1 <- function(a1tau1){
    a1 <<- exp(a1tau1[1])
    logittau1 <<- a1tau1[2]
    tau1 <<- plogis(logittau1)

    gloga1 = -a1*mean(prob_mat[,1]*(digamma(a1+R)-digamma(a1)+logittau1 -log(1+exp(logittau1))))

    glogtau1 = -mean(prob_mat[,1]*(a1 - (a1+R)*(tau1)))
    return(c(gloga1,glogtau1))
  }
  loga1tau1 <- c(log(a1),logittau1)
  optim_out <- optim(loga1tau1, object_a1,grad_a1,method = "L-BFGS-B")
  a1 <<- exp(optim_out$par[1])
  logittau1 <<- optim_out$par[2]
})


#' M-step for E-M algorithm to update lambda
CentipedeMixEm$methods(update.lambda = function(){
  lambda <<- colSums(exp( log_prob_vec - max_log)*X)/sum(exp( log_prob_vec - max_log)*X)
})

#' M-step for E-M algorithm to update tau0
CentipedeMixEm$methods(update.tau0 = function(){
  tau0 <<- exp(log(a0) + log(sum(prob_mat[,2])-log(sum(R*prob_mat[,2]))))
})


#' M-step for E-M algorithm to update tau1
CentipedeMixEm$methods(update.tau1 = function(){
  tau1 <<- exp(log(a1) + log(sum(exp( log_prob_vec - max_log)))-log(sum(R*exp( log_prob_vec - max_log))))
})

#' check convergence of mixture normal EM
CentipedeMixEm$methods(check.tol = function(fmax,fmin,ftol){
  delta = abs(fmax - fmin)
  accuracy = (abs(fmax) + abs(fmin))*ftol
  return(delta < (accuracy + tol))
})

#' main function for E-M algorithm
CentipedeMixEm$methods(run.EM = function(max_iter=1000L,loglik_tol=1e-5){
  convergence = 1L
  init.paras()  ## initialize parameter
  loglik_list = NULL
  for(iter in 1:max_iter){
    cat("iter:", iter)
    loglik0 <- loglik ## log-likelihood of previous steps
    update.prob() # E-step
    update.B()  ## M-step for beta's
    print(c("updated beta:",B))
    update.pi()   # M-step for pi_vec
    update.tau0()   # M-step for tau0
    update.tau1()   # M-step for tau1
    update.a0()   # M-step for a0
    update.a1()   # M-step for a1
    update.lambda()  # M-step for lambda
    loglik_list = c(loglik_list,loglik) # append log-likelihood
    if(check.tol(loglik0,loglik,loglik_tol)){
      convergence = 0 # converged
      break
    }
  }
  return(list(convergence=convergence,beta = B, a0 = a0, a1 = a1, tau0 = tau0, tau1 = tau1, lambda = lambda,
              pi_mat=pi_mat,iter=iter,loglik_list=loglik_list, bound_status = apply(prob_mat, 1, which.max)))
})








































