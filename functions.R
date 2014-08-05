##########################################################################################
## McGLM package version 0.1 - Multivariate Covariance Generalized Linear Models #########
## Author: Wagner Hugo Bonat LEG/IMADA ###################################################
## Date: 28/07/2014 ######################################################################
## Modified: 29/07/2014 ##################################################################
##########################################################################################

##########################################################################################
## Link functions ########################################################################
##########################################################################################

## Reciprocal for Gamma model
my.inverse <- function(beta, X, extra, derivative = TRUE){
  saida <- list()
  Xbeta = X%*%beta
  saida[[1]] <-  1/Xbeta
  Xbeta <- as.numeric(Xbeta)
  mu <- as.numeric(1/Xbeta)
  #mu <- as.numeric(saida[[1]])
  n <- length(mu)
  if(derivative == TRUE){saida[[2]] = t(Diagonal(n, -Xbeta^2)%*%X)}
  saida[[3]] <-  as.numeric(saida[[1]])
   return(saida)
}

## Inverse Gaussian
my.inv.gaussian <- function(beta, X, extra, derivative = TRUE){
  saida <- list()
  Xbeta = X%*%beta
  saida[[1]] <-  sqrt(1/Xbeta)
  mu <- as.numeric(saida[[1]])
  n <- length(mu)
  if(derivative == TRUE){saida[[2]] = t(Diagonal(n, -as.numeric(Xbeta)^3/2)%*%X)}
  saida[[3]] <-  as.numeric(saida[[1]])
   return(saida)
}
  
## Identity function for Gaussian data
my.identity <- function(beta,X,extra,derivative = TRUE){
  saida <- list()
  Xbeta = X%*%beta
  saida[[1]] <- Xbeta
  UM <- rep(1,length(Xbeta))
  n <- length(UM)
  if(derivative == TRUE){saida[[2]] = t(Diagonal(n, UM)%*%X)}
  saida[[3]] <- as.numeric(Xbeta)
   return(saida)
}
  
## Logit for Binomial data
my.logit <- function(beta, X, extra, derivative = TRUE){
  saida <- list()
  Xbeta = X%*%beta
  p <- exp(Xbeta)/(1+exp(Xbeta))
  p <- as.numeric(p)
  saida[[1]] <- extra*p
  n <- length(p)
  p.int <- p*(1-p)
  if(derivative == TRUE){saida[[2]] = t(extra*Diagonal(n,p.int)%*%X)}
  saida[[3]] <- p
   return(saida)
}

## Cauchit
my.cauchit <- function(beta,X,extra,derivative=TRUE){
  saida <- list()
  Xbeta = X%*%beta
  p <- as.numeric(cauchit(Xbeta,inverse=TRUE))
  saida[[1]] <- extra*p
  derivada <- as.numeric(cauchit(p, deriv=1))
  n <- length(p)
  if(derivative == TRUE){saida[[2]] = t(extra*Diagonal(n, derivada)%*%X)}
  saida[[3]] <- p
   return(saida)
}

## Probit
my.probit <- function(beta,X,extra,derivative=TRUE){
  saida <- list()
  Xbeta = X%*%beta
  p <- probit(as.numeric(Xbeta),inverse=TRUE)
  saida[[1]] <- extra*p
  derivada <- as.numeric(probit(p, deriv=1))
  n <- length(p)
  if(derivative == TRUE){saida[[2]] = t(extra*Diagonal(n, derivada)%*%X)}
  saida[[3]] <- p
   return(saida)
}

## Complementary log-log
my.cloglog <- function(beta,X,extra,derivative=TRUE){
  saida <- list()
  Xbeta = X%*%beta
  p <- as.numeric(cloglog(Xbeta,inverse=TRUE))
  saida[[1]] <- extra*p
  derivada <- as.numeric(cloglog(p, deriv=1))
  n <- length(p)
  if(derivative == TRUE){saida[[2]] = t(extra*Diagonal(n, derivada)%*%X)}
  saida[[3]] <- p
   return(saida)
}

## Log for Poisson data
my.log <- function(beta, X, extra, derivative = TRUE){
  saida <- list()
  Xbeta = X%*%beta
  if(is.null(extra)){saida[[1]] <- exp(Xbeta)}
  if(!is.null(extra)){saida[[1]] <- exp(Xbeta + log(extra))}
  n <- length(saida[[1]])
  mu <- as.numeric(saida[[1]])
  if(derivative == TRUE){saida[[2]] <- t(Diagonal(n,mu)%*%X)}
  saida[[3]] <- as.numeric(saida[[1]])
  return(saida)
}

## Square root link function
my.sqrt <- function(beta, X, extra, derivative = TRUE){
  saida <- list()
  Xbeta = X%*%beta
  if(is.null(extra)){saida[[1]] <- Xbeta^2}
  if(!is.null(extra)){saida[[1]] <- extra*Xbeta^2}
  mu <- as.numeric(saida[[1]])
  n <- length(mu)
  if(derivative == TRUE){saida[[2]] <- t(Diagonal(n, mu)%*%X)}
  saida[[3]] <- saida[[1]]
  return(saida)
}


##########################################################################################
## Generic functions to call the link functions ##########################################
##########################################################################################
link.function <- function(beta, X, link="log", derivative = FALSE, extra){
  if(link == "log"){saida = my.log(beta = beta, X=X, derivative = derivative, extra = extra)}
  if(link == "sqrt"){saida = my.sqrt(beta = beta, X=X, derivative = derivative, extra = extra)}
  if(link == "logit"){saida = my.logit(beta = beta, X=X, derivative = derivative, extra = extra)}
  if(link == "identity"){saida = my.identity(beta = beta, X=X, derivative = derivative, extra = extra)}
  if(link == "inverse"){saida = my.inverse(beta = beta, X=X, derivative = derivative, extra = extra)}
  if(link == "1/mu^2"){saida = my.inv.gaussian(beta = beta, X=X, derivative = derivative, extra = extra)}
  if(link == "cauchit"){saida = my.cauchit(beta = beta, X=X, derivative = derivative, extra = extra)}
  if(link == "probit"){saida = my.probit(beta = beta, X=X, derivative = derivative, extra = extra)}
  if(link == "cloglog"){saida = my.cloglog(beta = beta, X=X, derivative = derivative, extra = extra)}
  return(saida)
}

##########################################################################################
## Correlation matrix between response variables #########################################
##########################################################################################
monta.Sigmab <- function(n.resp, rho.vector){
    n.parameter <- n.resp*(n.resp - 1)/2
    Sigmab <- diag(n.resp)
    Sigmab[lower.tri(Sigmab)] <- rho.vector
    Sigmab <- forceSymmetric(t(Sigmab))
    return(Sigmab)
}

##########################################################################################
## Derivative of correlation matrix between response variables ###########################
##########################################################################################
D.Sigmab <- function(n.resp){
    position <- combn(n.resp,2)
    list.Derivative <- list()
    n.par <- n.resp*(n.resp-1)/2
    for(i in 1:n.par){
        Derivative <- matrix(0, ncol = n.resp, nrow = n.resp)
        Derivative[position[1,i],position[2,i]] <- Derivative[position[2,i],position[1,i]] <- 1
        list.Derivative[i][[1]] <- Derivative}
    return(list.Derivative)
}

##########################################################################################
## Complete derivative of C with respect rho #############################################
##########################################################################################
D.C.rho <- function(derivative, Id, c.Sw){
    saida <- c.Sw%*%kronecker(derivative,Id)%*%t(c.Sw)
    return(saida)
}

##########################################################################################
## Square root of V and its derivative with respect power parameter ######################
##########################################################################################
monta.V <- function(mu, power, type, n.obs, derivative = TRUE, power.fixed = FALSE){
    if(type == "identity"){
        V.sqrt <- Diagonal(n.obs, sqrt(1))
        saida <- list("V.sqrt" = V.sqrt)
    }

    if(type == "tweedie" & derivative == FALSE | power.fixed == TRUE){
        V.sqrt <- Diagonal(n.obs, sqrt(mu^power))
        saida <- list("V.sqrt" = V.sqrt)
    }
      
    if(type == "tweedie" & derivative == TRUE){
        mu.power <- mu^power
        V.sqrt <- Diagonal(n.obs, sqrt(mu.power))
        derivada <- (mu.power*log(mu))/(2*sqrt(mu.power))
        D.sqrt.V <- Diagonal(n.obs, derivada)
        saida <- list("V.sqrt" = V.sqrt, "D.sqrt.V" = D.sqrt.V)
    }

    if(type == "poisson.tweedie" & derivative == FALSE | power.fixed == TRUE){
        V.sqrt <- Diagonal(n.obs, sqrt(mu^power))
        saida <- list("V.sqrt" = V.sqrt)
    }
       
    if(type == "poisson.tweedie" & derivative == TRUE){
        mu.power <- mu^power
        V.sqrt <- Diagonal(n.obs, sqrt(mu.power))
        derivada <- (mu.power*log(mu))/(2*sqrt(mu.power))
        D.sqrt.V <- Diagonal(n.obs, derivada)
        saida <- list("V.sqrt" = V.sqrt, "D.sqrt.V" = D.sqrt.V)
    }

    if(type == "mu(1-mu)"){
       V.sqrt <- Diagonal(n.obs, sqrt(mu*(1-mu)))
       saida <- list("V.sqrt" = V.sqrt)}
   
    return(saida)
}

##########################################################################################
## Derivative of sqrt of V with respect to beta's ########################################
##########################################################################################
D.sqrt.V.beta <- function(beta, X, power, link, variance, extra){
    mu = link.function(beta = beta, X = X, link = link, derivative = TRUE, extra = extra)
    n.obs <- dim(X)[1]
    if(variance == "tweedie"){
        mu.p <- mu[[1]]^power
        mu.p1 <- mu[[1]]^(power-1)
        derivada <- (mu.p1*power)/(2*sqrt(mu.p))}
    if(variance == "poisson.tweedie"){
        mu.p <- mu[[1]]^power
        mu.p1 <- mu[[1]]^(power-1)
       derivada <- (mu.p1*power)/(2*sqrt(mu.p))}
    if(variance == "mu(1-mu)"){
        derivada <- (1 - 2*mu[[1]])/(2*sqrt(mu[[1]] - mu[[1]]^2))}
    if(variance == "identity"){
        derivada <- 0}
    if(variance == "mu"){
        derivada <- 1}
    saida <- list()
    n.beta <- length(beta)
    for(i in 1:n.beta){saida[[i]] <- Diagonal(n.obs, derivada*mu[[2]][i,])}
    return(saida)
}

##########################################################################################
## Derivative of Sigma Wi with respect to beta's #########################################
##########################################################################################
D.SigmaWi.beta <- function(derivative, Omega, V.sqrt){
    saida <- derivative%*%Omega%*%V.sqrt + V.sqrt%*%Omega%*%derivative
    return(saida)
}

##########################################################################################
## Derivative of cholesky Sigma  with respect to beta's ##################################
##########################################################################################
D.Sigma.beta <- function(mean.struc.hat, link, list.hat, variance, X.list, extra.list, Omega.args, power.fixed = FALSE){
    n.resp <- length(mean.struc.hat$list.mean)
    n.obs <- length(mean.struc.hat$list.mean[[1]][[1]])
    D.chol.beta <- list()
    for(i in 1:n.resp){
        #Omegai <- monta.Omega(tau.vector = list.hat$tau.list[[i]][[1]], n.obs = n.obs, derivative = TRUE)
        Omegai <- monta.Omega(tau.vector = list.hat$tau.list[[i]], Omega.args = Omega.args[[i]])
        sqrt.Vi <- monta.V(mu = mean.struc.hat$list.mean[[i]][[1]],
                           power = list.hat$power.vector[i], type = variance[i], n.obs = n.obs, derivative = TRUE, power.fixed = power.fixed)
        derivadas = D.sqrt.V.beta(list.hat$Regression[[i]], X = X.list[[i]], power = list.hat$power.vector[i], link = link[i],
                      variance = variance[i], extra = extra.list[[i]])
        S.Wi <- Sigma.Wi(V.sqrt = sqrt.Vi$V.sqrt, Omega = Omegai$Omega)
        if(variance[i] == "poisson.tweedie"){
            S.Wi <- Diagonal(n.obs, mean.struc.hat$list.mean[[i]][[1]]) + S.Wi
            derivadas2 <- D.sqrt.V.beta(list.hat$Regression[[i]], X = X.list[[i]], power = list.hat$power.vector[i], link = link[i],
                                        variance = "mu", extra = extra.list[[i]])
            for(j in 1:length(derivadas)){
            derivadas[[j]] <- derivadas[[j]] + derivadas2[[j]]}
        }
        if(class(S.Wi) == "dtCMatrix"){c.Swi <- Diagonal(n.obs, sqrt(diag(S.Wi)))}
        if(class(S.Wi) != "dtCMatrix"){c.Swi <- t(chol(S.Wi))}
        inv.c.Swi <- solve(c.Swi)
        ## Computing the derivatives
        D.Swi.beta <- lapply(derivadas, D.SigmaWi.beta, Omega = Omegai$Omega, V.sqrt = sqrt.Vi$V.sqrt)
        D.chol.beta[[i]] <- lapply(derivadas, D.chol, chol.S = c.Swi, inv.chol.S = inv.c.Swi)
    }
    D.beta <- monta.block.beta(n.resp = n.resp, n.obs = n.obs, list.derivatives = D.chol.beta)
    saida <- list("D.beta" = D.beta)
    return(saida)
}

##########################################################################################
## Derivatives of C with respect beta ####################################################
##########################################################################################
D.C.beta <- function(mean.struc.hat, link, list.hat, variance, X.list, extra.list, Omega.args, power.fixed = FALSE){
    n.resp <- length(mean.struc.hat$list.mean)
    n.obs <- length(mean.struc.hat$list.mean[[1]][[1]])
    Sb <- monta.Sigmab(n.resp = n.resp, rho.vector = list.hat$rho)
    I <- Diagonal(n.obs,1)
    kron.Sb <- kronecker(Sb,I)
    c.Sw <- monta.Sigma.block(list.mean = mean.struc.hat$list.mean,
                              power.vector = list.hat$power.vector, tau.list = list.hat$tau.list,
                              type = variance, Omega.args = Omega.args, power.fixed = FALSE)
    C <- c.Sw$Sigma$Sigma%*%kron.Sb%*%t(c.Sw$Sigma$Sigma)
    invC <- t(c.Sw$Sigma$inv.Sigma)%*%kronecker(solve(Sb),I)%*%c.Sw$Sigma$inv.Sigma
    D.beta <- D.Sigma.beta(mean.struc.hat = mean.struc.hat, link = link, list.hat = list.hat, variance = variance, X.list = X.list, extra.list = extra.list,
                           power.fixed = FALSE, Omega.args = Omega.args)
    D.C.b <- lapply(D.beta[[1]], deriv.prod, Sb = kron.Sb, c.Sw = c.Sw$Sigma$Sigma)
    return(D.C.b)
}

##########################################################################################
## Function to build block diagonal matrix from a list of matrices #######################
##########################################################################################
monta.block.beta <- function(n.resp, n.obs, list.derivatives){
    n.param <- sum(as.numeric(lapply(list.derivatives, length)))
    list.derivatives2 <- do.call(c, list.derivatives)
    saida <- list()
    mat0 <- Matrix(0, ncol = n.obs, nrow = n.obs, sparse = TRUE)
    for(i in 1:n.param){
        saida[[i]] <- list()}
    for(i in 1:n.param){
        for(k in 1:n.resp){
            saida[[i]][k][1] <- mat0
        }
    }
    posicao <- as.numeric(lapply(list.derivatives,length))
    posicao <- rep(1:n.resp,posicao)
    for(i in 1:n.param){
        saida[[i]][posicao[i]][[1]] <- list.derivatives2[[i]]}
    derivadas <- list()
    for(i in 1:n.param){
        derivadas[[i]] <- bdiag(saida[[i]])}
    return(derivadas)
}

##########################################################################################
## Omega matrix and its derivatives ######################################################
##########################################################################################
#monta.Omega <- function(tau.vector, n.obs, derivative = TRUE){
#    R1 <- Diagonal(n.obs, 1)
#    Omega <- tau.vector[1]*R1
#    if(derivative == TRUE){D.Omega <- R1}
#    if(derivative == FALSE){D.Omega <- NULL}
#    saida <- list("Omega" = Omega, "D.Omega" = list("tau1" = D.Omega))
#    return(saida)
#}

##########################################################################################
## Build SigmaWi #########################################################################
##########################################################################################
Sigma.Wi <- function(V.sqrt, Omega){
    S.Wi <- V.sqrt%*%Omega%*%V.sqrt
    return(S.Wi)
}

##########################################################################################
## Derivative of Sigma.wi with respect power parameter ###################################
##########################################################################################
D.SigmaWi.power <- function(derivative, Omega, V.sqrt){
    saida <- derivative%*%Omega%*%V.sqrt + V.sqrt%*%Omega%*%derivative
    return(saida)
}

##########################################################################################
## Derivative of Sigma.wi with respect phi and tau parameters ############################
##########################################################################################
D.SigmaWi.tau <- function(derivative, V.sqrt){
    saida <- V.sqrt%*%derivative%*%V.sqrt
    return(saida)
}

##########################################################################################
## Build block diagonal matrix and its derivative ########################################
##########################################################################################
monta.Sigma.block <- function(list.mean, power.vector, tau.list, type, Omega.args, power.fixed = FALSE){
    saida <- list()
    saida.inverse <- list()
    n.resp <- length(list.mean)
    n.obs <- length(list.mean[[1]][[1]])
    D.chol.power <- list()
    D.chol.tau <- list()
    cholesky <- list()
    inv.cholesky <- list()
    for(i in 1:n.resp){
        #Omegai <- monta.Omega(tau.vector = tau.list[[i]], n.obs = n.obs, derivative = TRUE)
        Omegai <- monta.Omega(tau.vector = tau.list[[i]], Omega.args = Omega.args[[i]])
        sqrt.Vi <- monta.V(mu = list.mean[[i]][[1]], power = power.vector[i], type = type[i], n.obs = n.obs, derivative = TRUE, power.fixed = power.fixed)
        S.Wi <- Sigma.Wi(V.sqrt = sqrt.Vi$V.sqrt, Omega = Omegai$Omega)
        if(type[i] == "poisson.tweedie"){S.Wi <- Diagonal(n.obs, list.mean[[i]][[1]]) + S.Wi}
        if(class(S.Wi) == "dtCMatrix"){c.Swi <- Diagonal(n.obs, sqrt(diag(S.Wi)))}
        if(class(S.Wi) != "dtCMatrix"){c.Swi <- t(chol(S.Wi))}
        inv.c.Swi <- solve(c.Swi)
        ## Computing the derivatives
        if(power.fixed == FALSE && type[i] == "tweedie" | type[i] == "poisson.tweedie"){
            D.Swi.power <- D.SigmaWi.power(derivative = sqrt.Vi$D.sqrt.V, Omega = Omegai$Omega, V.sqrt = sqrt.Vi$V.sqrt)
            D.chol.power[[i]] <- D.chol(derivative = D.Swi.power, chol.S = c.Swi, inv.chol.S = inv.c.Swi)}
        D.Swi.tau <- lapply(Omegai$D.Omega, D.SigmaWi.tau, V.sqrt = sqrt.Vi$V.sqrt)
        D.chol.tau[[i]] <- lapply(D.Swi.tau, D.chol, chol.S = c.Swi, inv.chol.S = inv.c.Swi)
        cholesky[[i]] <- c.Swi
        inv.cholesky[[i]] <- inv.c.Swi
    }
    Sigma.block <- bdiag(cholesky)
    inv.Sigma.block <- bdiag(inv.cholesky)
    D.tau <- monta.block.beta(n.resp = n.resp, n.obs = n.obs, list.derivatives = D.chol.tau)
    if(length(D.chol.power) != 0){D.power <- monta.block(n.resp = n.resp,n.obs=n.obs, list.derivatives =  D.chol.power)}
    if(length(D.chol.power) != 0){D.tau <- c(D.power,D.tau)}
    saida <- list("Sigma" = list("Sigma" = Sigma.block, "inv.Sigma.block" = inv.Sigma.block), "Derivatives" = list(D.tau))
    return(saida)
}

##########################################################################################
## Function to build block diagonal matrix from a list of matrices #######################
##########################################################################################
monta.block <- function(n.resp, n.obs, list.derivatives){
    idx = lapply(list.derivatives, is.null)
    n.param <- length(list.derivatives)
    saida <- list()
    mat0 <- Matrix(0, ncol = n.obs, nrow = n.obs, sparse = TRUE)
    for(i in 1:n.param){
        saida[[i]] <- list()}
    for(i in 1:n.param){
        for(k in 1:n.resp){
            saida[[i]][k][1] <- mat0
        }
    }
    for(i in 1:n.param){
        if(idx[i][[1]] == TRUE){saida[[i]] <- NA}
        if(idx[i][[1]] == FALSE){saida[[i]][i][1] <- list.derivatives[[i]]}
    }
    derivadas <- list()
    for(i in 1:n.param){
        if(idx[i][[1]] == TRUE){derivadas[[i]] <- NA}
        if(idx[i][[1]] == FALSE){derivadas[[i]] <- bdiag(saida[[i]])}
    }
    derivadas = derivadas[!is.na(derivadas)]
    return(derivadas)
}

##########################################################################################
## Derivative of Cholesky ################################################################
##########################################################################################
D.chol <- function(derivative, chol.S, inv.chol.S){
    mat1 <- inv.chol.S%*%derivative%*%t(inv.chol.S)
    diagonal.term <- diag(mat1)/2
    dimensao <- dim(inv.chol.S)
    mat2 <- Matrix(0,nrow = dimensao[1], ncol = dimensao[2])
    mat2[lower.tri(mat2)] <- mat1[lower.tri(mat1)]
    diag(mat2) <- diagonal.term
    saida <- chol.S%*%mat2
    return(saida)
}

##########################################################################################
## Generic derivative of product of matrices #############################################
##########################################################################################
deriv.prod <- function(derivative, Sb, c.Sw){
    saida <- derivative%*%Sb%*%t(c.Sw) + c.Sw%*%Sb%*%t(derivative)
    return(saida)
}

##########################################################################################
## Build C matrix ########################################################################
##########################################################################################
monta.C <- function(list.mean, power.vector, tau.list, rho.vector, type, Omega.args, power.fixed = FALSE){
    n.resp <- length(list.mean)
    n.obs <- length(list.mean[[1]][[1]])
    Sb <- monta.Sigmab(n.resp = n.resp, rho.vector = rho.vector)
    I <- Diagonal(n.obs,1)
    kron.Sb <- kronecker(Sb,I)
    c.Sw <- monta.Sigma.block(list.mean = list.mean, power.vector = power.vector, tau.list = tau.list, type = type,
                              Omega.args = Omega.args, power.fixed = FALSE)
    C <- c.Sw$Sigma$Sigma%*%kron.Sb%*%t(c.Sw$Sigma$Sigma)
    invC <- t(c.Sw$Sigma$inv.Sigma)%*%kronecker(solve(Sb),I)%*%c.Sw$Sigma$inv.Sigma
    D.C <- lapply(c.Sw$Derivatives[[1]], deriv.prod, Sb = kron.Sb, c.Sw = c.Sw$Sigma$Sigma)
    D.Sb.rho <- D.Sigmab(n.resp = n.resp)
    D.C.correlation <- lapply(D.Sb.rho, D.C.rho, Id = I, c.Sw = c.Sw$Sigma$Sigma)
    saida1 <- list("C" = C, "inv.C" = invC)
    saida <- list("Covariance" = saida1, "Derivatives" = c(D.C,D.C.correlation))
    return(saida)
}

##########################################################################################
## Generic score estimating function - Regression parameters #############################
##########################################################################################
quasi.score <- function(E.Y, Y, t.D, inv.C){
    e1 = t.D%*%inv.C
    esc <- e1%*%(Y-E.Y)
    Sensi <- e1%*%t(t.D)
    saida <- list("Score" = esc, "S.beta" = -Sensi)
    return(saida)
}
    
##########################################################################################
## Core of Pearson estimating function ###################################################
##########################################################################################
pearson.function <- function(d.Sigma, ## Derivative of covariance structure
                             res, ## Residuals vector
                             inv.Sigma ## inverse of covariance matrix
                             ){
    traco <- inv.Sigma%*%d.Sigma
    p1 <- t(res)%*%traco%*%inv.Sigma%*%res - sum(diag(traco))
    return(as.numeric(p1))
}
##########################################################################################
## Correct Pearson estimating function ###################################################
##########################################################################################
pearson.function.correct <- function(d.Sigma, ## Derivative of covariance structure
                                     res, ## Residuals vector
                                     inv.Sigma, ## inverse of covariance matrix
                                     Sigma, ## covariance matrix
                                     D ## derivative of mean
                             ){
    traco <- inv.Sigma%*%d.Sigma
    W.lambda <- traco%*%inv.Sigma
    J.beta <- t(D)%*%inv.Sigma%*%D
    term1 <- t(D)%*%W.lambda%*%D
    #correction <- -sum(diag(solve(term1,J.beta)))
    correction <- -sum(diag(solve(J.beta,term1)))
    p1 <- t(res)%*%W.lambda%*%res - sum(diag(traco)) + correction
    return(as.numeric(p1))
}

##########################################################################################
## Auxiliar function to multiply matrix ##################################################
##########################################################################################
multiplica <- function(d.Sigma, inv.Sigma){
    saida <- inv.Sigma%*%d.Sigma
    return(saida)
}

## Compute the weight matrix
multiplica.J.beta.parameter <- function(peso,D){
    saida <- t(D)%*%peso%*%D
    return(saida)
}

multiplica.weight <- function(d.Sigma, inv.Sigma){
    peso <- inv.Sigma%*%d.Sigma%*%inv.Sigma
    saida <- peso
    return(saida)
}

##########################################################################################
## Sensitive matrix ######################################################################
##########################################################################################
sensitive.covariance <- function(inv.Sigma, ## Inverse of covariance matrix
                                 d.Sigma.list ## list with all derivatives of var-cov matrix
                                 ){
      produto <- lapply(d.Sigma.list, multiplica, inv.Sigma = inv.Sigma)
      n.par <- length(d.Sigma.list)
      Sensitive <- matrix(NA, ncol = n.par, nrow = n.par)
      for(j in 1:n.par){
          for(i in 1:n.par){
              Sensitive[i,j] <- -sum(diag(produto[[i]]%*%produto[[j]]))
          }
      }
      return(Sensitive)
  }


##########################################################################################
## Pearson estimating function evaluating for whole sample ###############################
##########################################################################################
pearson.est.fun <- function(list.mean, type, power.vector, tau.list, rho.vector, res, Omega.args, correct = FALSE, D.block){
    MAT <- monta.C(list.mean = list.mean, type = type, power.vector = power.vector, tau.list = tau.list, rho.vector = rho.vector, Omega.args = Omega.args)
    if(correct == FALSE){
        escore <- do.call(c,lapply(MAT$Derivatives, pearson.function, res = res, inv.Sigma = MAT$Covariance$inv.C))
        Sensitive = sensitive.covariance(inv.Sigma = MAT$Covariance$inv.C, d.Sigma.list = MAT$Derivatives)
    }
    if(correct == TRUE){
        escore <- do.call(c,lapply(MAT$Derivatives, pearson.function.correct, res = res,
                                   Sigma = MAT$Covariance$C, inv.Sigma = MAT$Covariance$inv.C, D = t(D.block)))
        Sensitive = sensitive.covariance(inv.Sigma = MAT$Covariance$inv.C, d.Sigma.list = MAT$Derivatives)
    }
    #rm(MAT)
    #gc(reset = TRUE)
    saida <- list("escore" = escore, "Sensitive" = Sensitive)
    return(saida)
}

##########################################################################################
## Transform list.initial to vector of beta.ini and cov.ini ##############################
##########################################################################################
list2vec <- function(list.initial){
    beta.ini <- do.call(c,list.initial$Regression)
    cov.ini <- c(na.exclude(list.initial$power.vector), do.call(c,list.initial$tau.list), list.initial$rho)
    return(list("beta.ini" = as.numeric(beta.ini), "cov.ini" = as.numeric(cov.ini)))
}

getInformation <- function(list.initial){
    n.betas <- lapply(list.initial$Regression, length)
    n.taus <- lapply(list.initial$tau.list, length)
    n.power <- length(na.exclude(list.initial$power.vector))
    n.rho <- length(list.initial$rho)
    n.cov <- n.power+n.rho+sum(do.call(c,n.taus))
    saida <- list("n.betas" = n.betas, "n.taus" = n.taus, "n.power" = n.power, "n.rho" = n.rho, "n.cov" = n.cov)
    return(saida)
}


updatedBeta <- function(list.initial, betas, information){
    n.resp <- length(list.initial$Regression)
    cod <- rep(1:n.resp,information$n.betas)
    temp <- data.frame("beta" = betas,cod)
    for(k in 1:n.resp){list.initial$Regression[[k]] <- temp[which(temp$cod == k),]$beta}
    return(list.initial)
}

updatedCov <- function(list.initial, covariance, information){
    n.resp <- length(list.initial$Regression)
    #posicao.na <- which(is.na(list.initial$power.vector))
    posicao.true <- which(!is.na(list.initial$power.vector))
    list.initial$power.vector <- rep(NA, n.resp)
    if(information$n.power !=0){power.vector <- covariance[1:information$n.power]}
    #power.vector[variance == "tweedie"]
    cod.cov <- rep(1:n.resp,information$n.taus)
    n.cov <- information$n.cov - information$n.rho
    temp.cov <- data.frame("tau" = covariance[c(information$n.power+1):n.cov], cod.cov)
    for(k in 1:n.resp){list.initial$tau.list[[k]] <- temp.cov[which(temp.cov$cod == k),]$tau}
    list.initial$rho <- covariance[c(n.cov+1):information$n.cov]
    if(information$n.power !=0){list.initial$power.vector[posicao.true] <- power.vector}
    return(list.initial)
}
    

##########################################################################################
## Chaser algorithm ######################################################################
##########################################################################################
chaser <- function(list.initial, link, type, y.vector, X.list, extra.list, Omega.args, tol=0.0001, max.iter = 20, correct = FALSE){
   parametros <- list2vec(list.initial)
   inf <- getInformation(list.initial)
   solucao.beta <- matrix(NA, max.iter,length(parametros$beta.ini))
   solucao.cov <- matrix(NA, max.iter, length(parametros$cov.ini))
   solucao.beta[1,] <- parametros$beta.ini
   solucao.cov[1,] <- parametros$cov.ini
   beta.ini <- parametros$beta.ini
   cov.ini <- parametros$cov.ini
   n.resp <- length(list.initial$Regression)
    for(i in 2:max.iter){
      mean.struc <- mean.matrix(list.initial = list.initial, X.list = X.list, n.resp = n.resp, link = link, extra.list = extra.list)
      MAT = monta.C(list.mean = mean.struc$list.mean, power.vector = list.initial$power.vector,
                    tau.list = list.initial$tau.list, rho.vector = list.initial$rho, type = type,
                    Omega.args = Omega.args, power.fixed = power.fixed)
      beta.temp <- quasi.score(E.Y = mean.struc$mean.vector, Y = y.vector, t.D = mean.struc$D.block, inv.C =  MAT$Covariance$inv.C)
      solucao.beta[i,] <- as.numeric(beta.ini - solve(beta.temp$S.beta, beta.temp$Score))
      list.initial <- updatedBeta(list.initial, solucao.beta[i,], information = inf)
      mean.struc <- mean.matrix(list.initial = list.initial, X.list = X.list, n.resp = n.resp, link = link, extra.list = extra.list)      
      res = y.vector - mean.struc$mean.vector
      cov.temp <- pearson.est.fun(list.mean = mean.struc$list.mean, type = type, power.vector = list.initial$power.vector,
                           tau.list = list.initial$tau.list, rho.vector = list.initial$rho, res = res, correct = correct,
                                  Omega.args = Omega.args, D.block = mean.struc$D.block)
      solucao.cov[i,] <- as.numeric(cov.ini - solve(cov.temp$Sensitive, cov.temp$escore))
      beta.ini <- solucao.beta[i,]
      cov.ini <- solucao.cov[i,]
      print(round(cov.ini,3))
      list.initial <- updatedCov(list.initial, solucao.cov[i,], information = inf)
      tolera <- abs(c(solucao.beta[i,],solucao.cov[i,])  - c(solucao.beta[i-1,],solucao.cov[i-1,]))
    if(all(tolera < tol) == TRUE)break
  }
  saida <- list("IterationRegression" = solucao.beta, "IterarionCovariance" = solucao.cov, "Regression" = beta.ini, "Covariance" = cov.ini)
  return(saida)
}

##########################################################################################
## Variance-covariance matrix regression parameters ######################################
##########################################################################################
V.beta <- function(mean.struc, y.vector, MAT){
    S.beta = quasi.score(E.Y = mean.struc$mean.vector, Y = y.vector, t.D = mean.struc$D.block, inv.C = MAT$Covariance$inv.C)$S.beta
    V.b <- -S.beta
    return(V.b)
}
##########################################################################################
## Auxiliar function to compute cumulants ################################################
##########################################################################################
compute.cumulants <- function(res, mu, type, power.hat, phi.hat){
     if(type == "tweedie"){k4 = res^4 - 3*(phi.hat*mu^power.hat)^2}
     if(type == "poisson.tweedie"){k4 = res^4 - 3*(mu + phi.hat*mu^power.hat)^2}
     if(type == "mu(1-mu)"){k4 = res^4 - 3*(phi.hat*(mu*(1-mu)))^2}
     if(type == "identity"){k4 = res^4 - 3*(phi.hat)^2}
     return(k4)
 }

##########################################################################################
## Variance-covariance matrix covariance parameters ######################################
##########################################################################################
V.cov <- function(power.hat, phi.hat, mean.struc, y.list, MAT, type){
    n.obs <- dim(mean.struc$list.mean[[1]][1][[1]])[1]
    n.cov <- length(MAT$Derivatives)
    Var.cov <- matrix(NA, ncol = n.cov, nrow = n.cov)
    n.resp <- length(type)
    K4 <- matrix(NA, ncol = n.resp, nrow = n.obs)
    for(i in 1:n.resp){
        res = as.numeric(mean.struc$list.mean[[i]][1][[1]]) - y.list[[i]]
        K4[,i] <- compute.cumulants(res = res, mu = as.numeric(mean.struc$list.mean[[i]][1][[1]]), type = type[i],
                                    power.hat = power.hat[i], phi.hat = phi.hat[i])
    }
    k4 <- c(K4)
for(i in 1:n.cov){
 for(j in 1:n.cov){
    D.C1 <- MAT$Derivative[[i]]
    D.C2 <- MAT$Derivative[[j]]
    W1 <- MAT$Covariance$inv.C%*%D.C1%*%MAT$Covariance$inv.C
    W2 <- MAT$Covariance$inv.C%*%D.C2%*%MAT$Covariance$inv.C
    Var.cov[i,j] <- 2*sum(diag(MAT$Covariance$inv.C%*%D.C1%*%MAT$Covariance$inv.C%*%D.C2)) + sum(k4*diag(W1)*diag(W2))
}
}
    return(Var.cov)
}

##########################################################################################
## Auxiliary function to compute thirth moment ###########################################
##########################################################################################
covprod <- function(A, res, W){
    res =as.numeric(res)
    saida <- (res%*%W%*%res)%*%(t(res)%*%A)
    return(as.numeric(saida))
}

#########################################################################################
## Cross covariance matrix between regression and covariance parameters #################
##########################################################################################
V.c.b <- function(mean.struc, MAT, res){
    n.cov <- length(MAT$Derivatives)
    inv.C <- MAT$Covariance$inv.C
    Wlist <- list()
    for(i in 1:n.cov){
        Wlist[[i]] <- inv.C%*%MAT$Derivative[[i]]%*%inv.C}
    A = mean.struc$D.block%*%inv.C
    n.beta <- dim(A)[1]
    mat.V.c.b <- matrix(NA, ncol = n.cov, nrow = n.beta)
    for(j in 1:n.beta){
    for(i in 1:n.cov){
    mat.V.c.b[j,i] <- covprod(A[j,], Wlist[[i]], r = res)
}
}
    return(mat.V.c.b)
}

##########################################################################################
## Joint variance-covariance matrix ######################################################
##########################################################################################
VarCov <- function(mean.struc, list.hat, MAT, y.vector, y.list, type){
    V.b <- V.beta(mean.struc = mean.struc, y.vector = y.vector, MAT = MAT)
    power.hat <- list.hat$power.vector
    phi.hat <- sapply(list.hat$tau.list, "[[", 1)
    V.c <- V.cov(power.hat = list.hat$power.vector, phi.hat = phi.hat, mean.struc = mean.struc, y.list = y.list, MAT = MAT, type = type)
    res = y.vector - mean.struc$mean.vector
    V.cross <- V.c.b(mean.struc = mean.struc, MAT = MAT, res = res)  
    p1 <- rbind(as.matrix(V.b), as.matrix(t(V.cross)))
    p2 <- rbind(as.matrix(V.cross), as.matrix(t(V.c)))
    COV.mat <- cbind(p1,p2)
    return(COV.mat)
}
##########################################################################################
## Cross sensitivity matrix ##############################################################
##########################################################################################
S.cov.beta <- function(mean.struc.hat, link, list.hat, MAT, variance, X.list,extra.list, Omega.args, power.fixed = FALSE){
    D.beta = D.C.beta(mean.struc.hat = mean.struc.hat, link = link, list.hat = list.hat, variance = variance, X.list = X.list,
        extra.list = extra.list, Omega.args = Omega.args, power.fixed = FALSE)
    inv.C <- MAT$Covariance$inv.C
    info = getInformation(list.hat)
    n.beta = sum(as.numeric(lapply(list.hat$Regression,length)))
    S.c.b <- matrix(NA, ncol = n.beta, nrow = info$n.cov)
    for(i in 1:info$n.cov){
    for(j in 1:n.beta){
    Der.cov <- MAT$Derivative[[i]]
    Der.beta <- D.beta[[j]]
    S.c.b[i,j] <- -sum(diag(inv.C%*%Der.cov%*%inv.C%*%Der.beta))
}
}
return(S.c.b)
}

##########################################################################################
## Joint sensitivity matrix ##############################################################
##########################################################################################
inv.joint.sensitivity <- function(mean.struc.hat, link, list.hat, MAT, y.vector, variance, X.list, extra.list, Omega.args, correct){
    inv.S.beta <- solve(quasi.score(E.Y = mean.struc.hat$mean.vector , Y = y.vector, t.D = mean.struc.hat$D.block, inv.C = MAT$Covariance$inv.C)$S.beta)
    res <- y.vector - mean.struc.hat$mean.vector
    inv.S.cov = solve(pearson.est.fun(list.mean = mean.struc.hat$list.mean, type = variance, power.vector = list.hat$power.vector,
                      tau.list = list.hat$tau.list, rho.vector = list.hat$rho,res = res, correct = correct,
                      Omega.args = Omega.args, D.block = mean.struc.hat$D.block)$Sensitive)
    S.c.b = S.cov.beta(mean.struc.hat = mean.struc.hat, link = link, list.hat = list.hat, MAT = MAT, variance = variance,
        X.list = X.list, extra.list = extra.list, Omega.args = Omega.args)
    tt <- dim(S.c.b)
    mat0 <- matrix(0, ncol = tt[1], nrow = tt[2])
    cross.term <- -inv.S.cov%*%S.c.b%*%inv.S.beta
    p1 <- rbind(as.matrix(inv.S.beta),as.matrix(cross.term))
    p2 <- rbind(mat0, as.matrix(inv.S.cov))
    S.j <- cbind(p1,p2)
    return(as.matrix(S.j))
}

##########################################################################################
## Omega structured matrices #############################################################
##########################################################################################

monta.Omega <- function(tau.vector, Omega.args){
    mat <- f(par = tau.vector, formula = Omega.args$formula, models = Omega.args$models, unit.sample = Omega.args$unit.sample, data = Omega.args$data)
    Omega <- build.SXU(par = tau.vector, XU = mat)
    saida <- list("Omega" = Omega, "D.Omega" = mat)
    return(saida)
}

monta.Omega(tau.vector = c(100,0.8), Omega.args = Omega.args1)
