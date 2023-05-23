### Functions
tsmcpcox <- function(Y, Delta,X1,X2=NULL,Z, method = c("MCP", "SCAD"), c=25) {
  # tsmcpcox <- function(Y, Delta,X1,X2,Z, method = c("MCP", "SCAD"), c=25) {
  n <- length(Y)
  p <- dim(X1)[2]
  if (length(p) == 0) {
    p <- 0
  }
  p2<-dim(X2)[2]
  if (length(p2) == 0) {
    p2 <- 0
  }
  len<-seq(5,30,length=c) # kl, uniformly det point
  num.td<-NULL
  bicm<-NULL
  loc.td<-NULL
  x1_temp <-  X1[order(Z),]
  x2_temp<- X2[order(Z),]
  y_temp <- Y[order(Z)]
  delta_temp <- Delta[order(Z)]
  z_temp<-sort(Z) # ordered Z
  
  # l=10
  for(l in 1:c){
    
    m<-ceiling(len[l]*n^(1/5)) # m observation in total
    q <- floor(n/m) # starting qn
    f<- n-m*(q-1) # length of the first segment
    K_temp <- matrix(0, nrow = q, ncol = q, byrow = TRUE) # q * q matrix
    
    for (i in 1:q) K_temp[i, 1:i] <- rep(1, i)
    x1 <- NULL
    x2 <- NULL
    y <- NULL
    delta <- NULL
    z <- NULL
    # select the segment's covariate
    x1[[1]] <- as.matrix(x1_temp[1:((n - (q - 1) * m)), ]) 
    y[[1]] <- y_temp[1:((n - (q - 1) * m))]
    delta[[1]] <- delta_temp[1:((n - (q - 1) * m))]
    z[[1]] <- z_temp[1:((n - (q - 1) * m))]
    for (i in 2:q) {
      # each group
      x1[[i]] <- as.matrix(x1_temp[(n - (q - i + 1) * m + 1):((n - (q -
                                                                      i) * m)), ])
      y[[i]] <- y_temp[(n - (q - i + 1) * m + 1):((n - (q - i) * m))]
      delta[[i]] <- delta_temp[(n - (q - i + 1) * m + 1):((n - (q - i) * m))]
      z[[i]] <- z_temp[(n - (q - i + 1) * m + 1):((n - (q - i) * m))]
    }
    if(length(x2_temp)!=0)
    {
      x2[[1]] <- as.matrix(x2_temp[1:((n - (q - 1) * m)), ])
      for (i in 2:q) {
        x2[[i]] <- as.matrix(x2_temp[(n - (q - i + 1) * m + 1):((n - (q -  i) * m)), ])
      }
    }
    
    
    X_temp1 <- lapply(1:length(x1), function(j, mat, list) kronecker(mat[j,, drop = FALSE], list[[j]]), mat = K_temp, list = x1)
    Xn <- cbind(do.call("rbind", X_temp1),x2_temp) # n*2q matrix
    
    #########################################adlasso######################
    
    y_surv <- survival::Surv(y_temp,delta_temp) # survival time data structure
    
    object <- tune.fit1(Xn, y_surv, family="cox",q,m,f, penalty = method, tune = "bic") # estimate once
    # later:     object <- tune.fit(Xn, y_surv, family="cox",penalty = method, tune = "bic") # system to refine
    
    
    coeff <- rep(0,p*q+p2)
    
    coeff[object$ix] <- object$beta
    
    adcp.coef.s <- sum(abs(coeff))
    
    adcp.coef.v.m <- abs(matrix(c(coeff[1:(p*q)]), q, p , byrow = T))
    
    adcp.coef.m <- c(apply(adcp.coef.v.m, 1, max))
    
    adcp.cp <- which(adcp.coef.m != 0)
    
    if (length(adcp.cp) > 1) {
      for (i in 2:length(adcp.cp)) {
        if (adcp.cp[i] - adcp.cp[i - 1] == 1)
          adcp.cp[i] <- 0
      }
    }
    
    
    adcp.cp1 <- adcp.cp[adcp.cp > 1 & adcp.cp < q]
    
    d1 <- length(adcp.cp1)
    
    if (d1 == 0) {
      adcpcss.cp <- integer(0)
    }
    
    if (d1 >= 1) {
      # step 4 find change points. algorithm line 7
      adcpcss.cp <- NULL
      
      adcp.cp1 <- c(0, adcp.cp1, q + 1)
      for (i in 1:d1) {
        
        yt <- NULL
        x1t <- NULL
        x2t<-NULL
        deltat<-NULL
        zt<-NULL
        for (k in (adcp.cp1[i + 1] - 1):(adcp.cp1[i + 1])) 
        {
          yt <- c(yt, y[[k]])
          zt<-c(zt,z[[k]])
          deltat<-c(deltat,delta[[k]])
          x1t <- rbind(x1t, as.matrix(x1[[k]]))
          
          
        }
        if(length(x2_temp)!=0){
          for (k in (adcp.cp1[i + 1] - 1):(adcp.cp1[i + 1])){
            x2t <- rbind(x2t, as.matrix(x2[[k]]))
          }
        }
        
        # using profile-likelihood find change point.
        cp <- order(abs(z_temp-css(x1t,x2t, zt,yt,deltat,1.2*m)))[1]  # equivalent == which.min
        
        if (length(cp) == 0)
          next
        if (cp != 0) {
          if (adcp.cp1[i + 1] == 0)
            adcpcss.cp <- c(adcpcss.cp, cp)
          if (adcp.cp1[i + 1] > 0)
            adcpcss.cp <- c(adcpcss.cp, cp)
          
        }
        
        
      }
      
      if (length(adcpcss.cp) == 0)
        adcpcss.cp <- integer(0) else adcpcss.cp <- adcpcss.cp
      
    }
    
    tt <- which(abs(diff(adcpcss.cp)) < 10)
    if (length(tt) > 0)
      adcpcss.cp <- adcpcss.cp[-tt]
    
    
    #return(c(length(adcpcss.cp)))
    num.td[l]<-length(adcpcss.cp)
    adcpcss.cp1<-c(0,adcpcss.cp,n)
    loc.td[[l]]<-adcpcss.cp1
    q <- length(adcpcss.cp1)-1
    K_temp <- matrix(0, nrow = q, ncol = q, byrow = TRUE)
    
    for (i in 1:q) K_temp[i, 1:i] <- rep(1, i)
    x1 <- NULL
    for (i in 1:q) {
      x1[[i]] <- as.matrix(x1_temp[(adcpcss.cp1[i]+1):adcpcss.cp1[i+1], ])
    }
    
    
    X1_temp1 <- lapply(1:length(x1), function(j, mat, list) kronecker(mat[j,
                                                                          , drop = FALSE], list[[j]]), mat = K_temp, list = x1)
    Xn <- cbind(do.call("rbind", X1_temp1),x2_temp)
    
    ##### cox reg--> 2.4 selection of the tuning parameter and segment length
    
    coxfit<-coxph(Surv(Y,Delta)~ Xn-1) # algorithm step 9
    coeff<-coxfit$coefficients
    likvalue<-coxfit$loglik[2]
    bicm[l]<--2*likvalue+(p*(num.td[l]+1)+p2)*log(n) # calculate bic under each l
    
  }
  
  l.opt<-min(which(bicm==min(bicm))) # find the l that minize the bic
  m.opt<-ceiling(l.opt*n^(1/5))
  loc.opt<-loc.td[[l.opt]]
  num.opt<-num.td[l.opt]
  adcpcss.cp<-loc.opt[loc.opt>1&loc.opt<n]
  
  # if the estimation is correct question: what if tau unknown?????
  # if(length(loc.opt)==(length(tau)+2)){
  #   adcpcss.cp1<-loc.opt
  #   q <- length(adcpcss.cp1)-1
  #   K_temp <- matrix(0, nrow = q, ncol = q, byrow = TRUE)
  # 
  #   for (i in 1:q) K_temp[i, 1:i] <- rep(1, i)
  # 
  # 
  #   x1 <- NULL
  #   for (i in 1:q) {
  #     x1[[i]] <- as.matrix(x1_temp[(adcpcss.cp1[i]+1):adcpcss.cp1[i+1], ])
  #       }
  # 
  #   X_temp1 <- lapply(1:length(x1), function(j, mat, list) kronecker(mat[j, , drop = FALSE], list[[j]]), mat = K_temp, list = x1)
  #   Xn <- cbind(do.call("rbind", X_temp1),x2_temp)
  #   y_surv <- survival::Surv(y_temp,delta_temp)
  #   object <- tune.fit(Xn, y_surv, family="cox",penalty = method, tune = "bic")
  #   coeff<-rep(0,p*q+p2)
  # 
  #   coeff[object$ix] <-object$beta
  #   result<-list(
  #    change.points.num=length(Z[adcpcss.cp]),
  #    change.points=Z[adcpcss.cp],
  #    coeff=coeff,
  #    m=m.opt
  #    )
  #   return(result)}else{
  #     return( result<-list(
  #       change.points.num=0,
  #       change.points=NULL,
  #       coeff=coeff,
  #       m=m.opt
  #     ))
  #   }
  browser()
  
  ### new changed
  adcpcss.cp1<-loc.opt
  q <- length(adcpcss.cp1)-1
  K_temp <- matrix(0, nrow = q, ncol = q, byrow = TRUE)
  
  for (i in 1:q) K_temp[i, 1:i] <- rep(1, i)
  
  
  x1 <- NULL
  for (i in 1:q) {
    x1[[i]] <- as.matrix(x1_temp[(adcpcss.cp1[i]+1):adcpcss.cp1[i+1], ])
  }
  
  X_temp1 <- lapply(1:length(x1), function(j, mat, list) kronecker(mat[j, , drop = FALSE], list[[j]]), mat = K_temp, list = x1)
  Xn <- cbind(do.call("rbind", X_temp1),x2_temp)
  y_surv <- survival::Surv(y_temp,delta_temp)
  object <- tune.fit(Xn, y_surv, family="cox",penalty = method, tune = "bic")
  coeff<-rep(0,p*q+p2)
  
  coeff[object$ix] <-object$beta
  result<-list(
    change.points.num=length(Z[adcpcss.cp]),
    change.points=Z[adcpcss.cp],
    coeff=coeff,
    m=m.opt,
    bic.opt = bicm[l.opt]
  )
  return(result)
}


lamNames <- function(l) {
  if (length(l) > 1) {
    d <- ceiling(-log10(-max(diff(l))))
    d <- min(max(d,4), 10)
  } else {
    d <- 4
  }
  formatC(l, format="f", digits=d)
}

convexMin <- function(b, X, penalty, gamma, l2, family, penalty.factor, a, Delta=NULL) {
  n <- nrow(X)
  p <- ncol(X)
  l <- ncol(b)
  
  if (penalty=="MCP") {
    k <- 1/gamma
  } else if (penalty=="SCAD") {
    k <- 1/(gamma-1)
  } else if (penalty=="lasso") {
    return(NULL)
  }
  if (l==0) return(NULL)
  
  val <- NULL
  for (i in 1:l) {
    A1 <- if (i==1) rep(1,p) else b[,i]==0
    if (i==l) {
      L2 <- l2[i]
      U <- A1
    } else {
      A2 <- b[,i+1]==0
      U <- A1&A2
      L2 <- l2[i+1]
    }
    if (sum(!U)==0) next
    Xu <- X[,!U]
    p.. <- k*(penalty.factor[!U]!=0) - L2*penalty.factor[!U]
    if (family=="gaussian") {
      if (any(A1!=A2)) {
        eigen.min <- min(eigen(crossprod(Xu)/n - diag(p.., length(p..), length(p..)))$values)
      }
    } else if (family=="binomial") {
      if (i==l) eta <- a[i] + X%*%b[,i]
      else eta <- a[i+1] + X%*%b[,i+1]
      pi. <- exp(eta)/(1+exp(eta))
      w <- as.double(pi.*(1-pi.))
      w[eta > log(.9999/.0001)] <- .0001
      w[eta < log(.0001/.9999)] <- .0001
      Xu <- sqrt(w) * cbind(1, Xu)
      xwxn <- crossprod(Xu)/n
      eigen.min <- min(eigen(xwxn-diag(c(0, diag(xwxn)[-1]*p..)))$values)
    } else if (family=="poisson") {
      if (i==l) eta <- a[i] + X%*%b[,i]
      else eta <- a[i+1] + X%*%b[,i+1]
      mu <- exp(eta)
      w <- as.double(mu)
      Xu <- sqrt(w) * cbind(1, Xu)
      xwxn <- crossprod(Xu)/n
      eigen.min <- min(eigen(xwxn-diag(c(0, diag(xwxn)[-1]*p..)))$values)
    } else if (family=="cox") {
      eta <- if (i==l) X%*%b[,i] else X%*%b[,i+1]
      haz <- drop(exp(eta))
      rsk <- rev(cumsum(rev(haz)))
      h <- haz*cumsum(Delta/rsk)
      xwxn <- crossprod(sqrt(h) * Xu)/n
      eigen.min <- min(eigen(xwxn-diag(diag(xwxn)*p.., nrow(xwxn), ncol(xwxn)))$values)
    }
    
    if (eigen.min < 0) {
      val <- i
      break
    }
  }
  val
}

std <- function(X) {
  if (typeof(X) == 'integer') storage.mode(X) <- 'double'
  if (!inherits(X, "matrix")) {
    if (is.numeric(X)) {
      X <- matrix(as.double(X), ncol=1)
    } else {
      tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
    }
  }
  STD <- .Call("standardize", X)
  dimnames(STD[[1]]) <- dimnames(X)
  ns <- which(STD[[3]] > 1e-6)
  if (length(ns) == ncol(X)) {
    val <- STD[[1]]
  } else {
    val <- STD[[1]][, ns, drop=FALSE]
  }
  attr(val, "center") <- STD[[2]]
  attr(val, "scale") <- STD[[3]]
  attr(val, "nonsingular") <- ns
  val
}

#penalty.factor<-rep(1, ncol(X))
setupLambdaCox <- function(X, y, Delta, alpha, lambda.min, nlambda, penalty.factor) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Fit to unpenalized covariates
  ind <- which(penalty.factor!=0)
  if (length(ind)!=p) {
    nullFit <- survival::coxph(survival::Surv(y, Delta) ~ X[, -ind, drop=FALSE])
    eta <- nullFit$linear.predictors
    rsk <- rev(cumsum(rev(exp(eta))))
    s <- Delta - exp(eta)*cumsum(Delta/rsk)
  } else {
    w <- 1/(n-(1:n)+1)
    s <- Delta - cumsum(Delta*w)
  }
  
  # Determine lambda.max
  zmax <- .Call("maxprod", X, s, ind, penalty.factor) / n
  lambda.max <- zmax/alpha
  
  if (lambda.min==0) lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), len=nlambda-1)), 0)
  else lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), len=nlambda))
  lambda
}

# original
ncvsurv1 <- function(X, y, qn,mn,fn,penalty="SCAD", gamma=switch(penalty, SCAD=3.7, 3),
                     alpha=1, lambda.min=ifelse(n>p,.001,.05), nlambda=100, lambda, eps=1e-4, max.iter=10000,
                     convex=TRUE, dfmax=p, penalty.factor=rep(1, ncol(X)), warn=TRUE, returnX, ...) {
  
  # with change
  # ncvsurv1 <- function(X, y, qn,mn,fn,penalty="SCAD", gamma=switch(penalty, SCAD=3.7, 3),
  #                        alpha=1, lambda.min=0.001, nlambda=100, lambda, eps=1e-4, max.iter=10000,
  #                        convex=TRUE, dfmax=p, penalty.factor=rep(1, ncol(X)), warn=TRUE, returnX, ...){
  
  # Coersion
  
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (!inherits(y, "matrix")) {
    tmp <- try(y <- as.matrix(y), silent=TRUE)
    if (inherits(tmp, "try-error")) stop("y must be a matrix or able to be coerced to a matrix", call.=FALSE)
    if (ncol(y) != 2) stop("y must have two columns for survival data: time-on-study and a censoring indicator", call.=FALSE)
  }
  if (typeof(y) == "integer") storage.mode(y) <- "double"
  if (typeof(penalty.factor) != "double") storage.mode(penalty.factor) <- "double"
  storage.mode(qn) <- "integer"
  storage.mode(mn) <- "integer"
  storage.mode(fn) <- "integer"
  
  ## Error checking
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty", call.=FALSE)
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty", call.=FALSE)
  if (nlambda < 2) stop("nlambda must be at least 2", call.=FALSE)
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead", call.=FALSE)
  if (length(penalty.factor)!=ncol(X)) stop("penalty.factor does not match up with X", call.=FALSE)
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg", call.=FALSE)
  
  ## Set up XX, yy, lambda
  yy<-sort(y[1:fn,1])
  tOrder<-order(y[1:fn,1])
  Delta<-y[1:fn,2][tOrder]
  XX1<-X[1:fn, , drop=FALSE][tOrder,,drop=FALSE]
  for(i in 1:(qn-1)){
    yy<- c(yy,sort(y[(fn+mn*(i-1)+1):(fn+mn*i),1]))
    tOrder<-order(y[(fn+mn*(i-1)+1):(fn+mn*i),1])
    Delta<-c(Delta,y[(fn+mn*(i-1)+1):(fn+mn*i),2][tOrder])
    XX1<-rbind(XX1,X[(fn+mn*(i-1)+1):(fn+mn*i), , drop=FALSE][tOrder,,drop=FALSE])
  }
  XX<-std(XX1)
  
  if (sys.nframe() > 1 && sys.call(-1)[[1]]=="local_mfdr") return(list(X=XX, time=yy, fail=Delta))
  ns <- attr(XX, "nonsingular")
  penalty.factor <- penalty.factor[ns]
  p <- ncol(XX)
  if (missing(lambda)) {
    lambda <- setupLambdaCox(XX, yy, Delta, alpha, lambda.min, nlambda, penalty.factor)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    if (nlambda == 1) warning(gsub('ncvreg', 'ncvsurv', lambda.warning), call.=FALSE)
    user.lambda <- TRUE
  }
  
  ## Fit
  res <- .Call("strat_cox_dh", XX, Delta,qn,mn,fn,penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor,
               alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)), as.integer(warn))
  b <- matrix(res[[1]], p, nlambda)
  loss <- -1*res[[2]]
  iter <- res[[3]]
  Eta <- matrix(res[[4]], n, nlambda)
  
  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  Eta <- Eta[, ind, drop=FALSE]
  if (warn & sum(iter)==max.iter) warning("Algorithm failed to converge for some values of lambda")
  
  ## Local convexity?
  convex.min <- if (convex) convexMin(b, XX, penalty, gamma, lambda*(1-alpha), "cox", penalty.factor, Delta=Delta) else NULL
  
  ## Unstandardize
  beta <- matrix(0, nrow=ncol(X), ncol=length(lambda))
  bb <- b/attr(XX, "scale")[ns]
  beta[ns,] <- bb
  offset <- -crossprod(attr(XX, "center")[ns], bb)
  
  ## Names
  varnames <- if (is.null(colnames(X))) paste("V", 1:ncol(X), sep="") else colnames(X)
  dimnames(beta) <- list(varnames, lamNames(lambda))
  
  ## Output
  val <- structure(list(beta = beta,
                        iter = iter,
                        lambda = lambda,
                        penalty = penalty,
                        gamma = gamma,
                        alpha = alpha,
                        convex.min = convex.min,
                        loss = loss,
                        penalty.factor = penalty.factor,
                        n = n,
                        time = yy,
                        fail = Delta,
                        order = tOrder),
                   class = c("ncvsurv", "ncvreg"))
  val$Eta <- sweep(Eta, 2, offset, "-")
  if (missing(returnX)) {
    if (utils::object.size(XX) > 1e8) {
      warning("Due to the large size of X (>100 Mb), returnX has been turned off.\nTo turn this message off, explicitly specify returnX=TRUE or returnX=FALSE).")
      returnX <- FALSE
    } else {
      returnX <- TRUE
    }
  }
  if (returnX) {
    val$X <- XX
  }
  val
}

####################################################################

loglik = function(X, y, beta, family) {
  K = dim(beta)[2]
  link = cbind(1, X) %*% beta
  yrep = repmat(y, 1, K)
  if (family == "gaussian")
    return(apply((yrep - link)^2, 2, sum))
  if (family == "poisson")
    return(apply(exp(link) - yrep * link, 2, sum))
  if (family == "binomial")
    return(apply(log(1 + exp(link)) - yrep * link, 2, sum))
}

repmat = function(X, m, n) {
  ## R equivalent of repmat (matlab)
  X = as.matrix(X)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = T)
}

getdf = function(coef.beta) {
  apply(abs(coef.beta) > 1e-10, 2, sum)
}


margcoef <- function(x, y, condind = NULL, family, null.model = FALSE, iterind) {
  n = dim(x)[1]
  p = dim(x)[2]
  ones = rep(1, n)
  candind = setdiff(1:p, condind)
  if (iterind == 0) {
    if (family == "cox")
      margcoef = abs(cor(x, y[, 1])) else margcoef = abs(cor(x, y))
  } else {
    if (null.model == TRUE) {
      if (is.null(condind) == TRUE) {
        x = x[sample(1:n), ]
      }
      if (is.null(condind) == FALSE) {
        x[, candind] = x[sample(1:n), candind]
      }
    }
    margcoef = abs(sapply(candind, mg, x, y, ones, family, condind))
  }
  return(margcoef)
}

mg <- function(index, x = x, y = y, ones = ones, family = family, condind = condind) {
  margfit = switch(family, gaussian = coef(glm.fit(cbind(ones, x[, index], x[, condind]), y, family = gaussian()))[2],
                   binomial = coef(glm.fit(cbind(ones, x[, index], x[, condind]), y, family = binomial()))[2], poisson = coef(glm.fit(cbind(ones,
                                                                                                                                            x[, index], x[, condind]), y, family = poisson()))[2], cox = coef(coxph(y ~ cbind(x[, index], x[,
                                                                                                                                                                                                                                            condind])))[1])
}


obtain.ix0 <- function(x, y, s1, s2, family, nsis, iter, varISIS, perm, q, greedy, greedy.size, iterind) {
  if (iter == FALSE) {
    margcoef = margcoef(x, y, family = family, null.model = FALSE, iterind = iterind)
    rankcoef = sort(margcoef, decreasing = TRUE, index.return = TRUE)
    ix0 = rankcoef$ix[1:nsis]
  } else {
    if (varISIS == "vanilla") {
      margcoef = margcoef(x, y, family = family, null.model = FALSE, iterind = iterind)
      rankcoef = sort(margcoef, decreasing = TRUE, index.return = TRUE)
      if (perm == FALSE)
        ix0 = rankcoef$ix[1:floor((2/3) * nsis)] else {
          count = 0
          repeat {
            count = count + 1
            randcoef = margcoef(x, y, family = family, null.model = TRUE, iterind = iterind)
            if (length(which(margcoef >= quantile(randcoef, q))) > 0)
              break
            if(count > 10)
              break
          }
          if (greedy == FALSE) {
            if (length(which(margcoef >= quantile(randcoef, q))) >= 2) {
              length1 = length(which(margcoef >= quantile(randcoef, q)))
              above.thresh = rankcoef$ix[1:length1]
              ix0 = rankcoef$ix[1:floor((2/3) * nsis)]
              ix0 = sort(intersect(ix0, above.thresh))
            } else ix0 = rankcoef$ix[1:2]
          } else {
            if (greedy.size == 1)
              ix0 = rankcoef$ix[1:2] else ix0 = rankcoef$ix[1:greedy.size]
          }
        }
    } else {
      if(family == 'cox'){
        margcoef1 = margcoef(x[s1, ], y[s1,], family = family, null.model = FALSE, iterind = iterind)
        margcoef2 = margcoef(x[s2, ], y[s2,], family = family, null.model = FALSE, iterind = iterind)
      } else{
        margcoef1 = margcoef(x[s1, ], y[s1], family = family, null.model = FALSE, iterind = iterind)
        margcoef2 = margcoef(x[s2, ], y[s2], family = family, null.model = FALSE, iterind = iterind)
        
      }
      
      rankcoef1 = sort(margcoef1, decreasing = TRUE, index.return = TRUE)
      rankcoef2 = sort(margcoef2, decreasing = TRUE, index.return = TRUE)
      if (perm == FALSE) {
        if (varISIS == "aggr") {
          ix01 = rankcoef1$ix[1:floor((2/3) * nsis)]
          ix02 = rankcoef2$ix[1:floor((2/3) * nsis)]
          ix0 = sort(intersect(ix01, ix02))
          if (length(ix0) <= 1)
            ix0 = int.size.k(rankcoef1$ix, rankcoef2$ix, 2)
        }
        if (varISIS == "cons") {
          iensure = intensure(floor((2/3) * nsis), l1 = rankcoef1$ix, l2 = rankcoef2$ix, k = floor((2/3) * nsis))
          ix01 = rankcoef1$ix[1:iensure]
          ix02 = rankcoef2$ix[1:iensure]
          ix0 = sort(intersect(ix01, ix02))
        }
      } else {
        count = 0
        repeat {
          count  = count + 1
          randcoef1 = margcoef(x[s1, ], y[s1], family = family, null.model = TRUE, iterind = iterind)
          randcoef2 = margcoef(x[s2, ], y[s2], family = family, null.model = TRUE, iterind = iterind)
          if (length(which(margcoef1 >= quantile(randcoef1, q))) > 0 && length(which(margcoef2 >= quantile(randcoef2,
                                                                                                           q))) > 0)
            break
          if(count > 10) break
        }
        if (greedy == FALSE) {
          length1 = length(which(margcoef1 >= quantile(randcoef1, q)))
          length2 = length(which(margcoef2 >= quantile(randcoef2, q)))
          above.thresh.1 = rankcoef1$ix[1:length1]
          above.thresh.2 = rankcoef2$ix[1:length2]
          ix01 = rankcoef1$ix[1:floor((2/3) * nsis)]
          ix02 = rankcoef2$ix[1:floor((2/3) * nsis)]
          ix01 = sort(intersect(ix01, above.thresh.1))
          ix02 = sort(intersect(ix02, above.thresh.2))
          ix0 = sort(intersect(ix01, ix02))
          if (length(ix0) <= 1)
            ix0 = int.size.k(rankcoef1$ix, rankcoef2$ix, 2)
        } else {
          if (greedy.size == 1)
            ix0 = int.size.k(rankcoef1$ix, rankcoef2$ix, 2) else ix0 = int.size.k(rankcoef1$ix, rankcoef2$ix, greedy.size)
        }
      }
    }
  }
  return(ix0)
}

obtain.newix <- function(x, y, ix1, candind, s1, s2, family, pleft, varISIS, perm, q, greedy, greedy.size, iterind) {
  if (varISIS == "vanilla") {
    margcoef = margcoef(x, y, ix1, family = family, null.model = FALSE, iterind = iterind)
    rankcoef = sort(margcoef, decreasing = TRUE, index.return = TRUE)
    if (perm == FALSE) {
      if (pleft > 0)
        newix = candind[rankcoef$ix[1:pleft]] else newix = NULL
    } else {
      randcoef = margcoef(x, y, ix1, family = family, null.model = TRUE, iterind = iterind)
      if (length(which(margcoef >= quantile(randcoef, q))) > 0) {
        if (greedy == FALSE) {
          length1 = length(which(margcoef >= quantile(randcoef, q)))
          above.thresh = candind[rankcoef$ix[1:length1]]
          newix = candind[rankcoef$ix[1:pleft]]
          newix = sort(intersect(newix, above.thresh))
        } else newix = candind[rankcoef$ix[1:greedy.size]]
      } else newix = NULL
    }
  } else {
    margcoef1 = margcoef(x[s1, ], y[s1], ix1, family = family, null.model = FALSE, iterind = iterind)
    margcoef2 = margcoef(x[s2, ], y[s2], ix1, family = family, null.model = FALSE, iterind = iterind)
    rankcoef1 = sort(margcoef1, decreasing = TRUE, index.return = TRUE)
    rankcoef2 = sort(margcoef2, decreasing = TRUE, index.return = TRUE)
    if (perm == FALSE) {
      if (pleft > 0) {
        if (varISIS == "aggr") {
          newix1 = candind[rankcoef1$ix[1:pleft]]
          newix2 = candind[rankcoef2$ix[1:pleft]]
          newix = sort(intersect(newix1, newix2))
        }
        if (varISIS == "cons") {
          iensure = intensure(pleft, l1 = rankcoef1$ix, l2 = rankcoef2$ix, k = pleft)
          newix1 = candind[rankcoef1$ix[1:iensure]]
          newix2 = candind[rankcoef2$ix[1:iensure]]
          newix = sort(intersect(newix1, newix2))
        }
      } else newix = NULL
    } else {
      randcoef1 = margcoef(x[s1, ], y[s1], ix1, family = family, null.model = TRUE, iterind = iterind)
      randcoef2 = margcoef(x[s2, ], y[s2], ix1, family = family, null.model = TRUE, iterind = iterind)
      if (length(which(margcoef1 >= quantile(randcoef1, q))) > 0 && length(which(margcoef2 >= quantile(randcoef2,                                                                          q))) > 0) {
        if (greedy == FALSE) {
          length1 = length(which(margcoef1 >= quantile(randcoef1, q)))
          length2 = length(which(margcoef2 >= quantile(randcoef2, q)))
          above.thresh.1 = candind[rankcoef1$ix[1:length1]]
          above.thresh.2 = candind[rankcoef2$ix[1:length2]]
          newix1 = candind[rankcoef1$ix[1:pleft]]
          newix2 = candind[rankcoef2$ix[1:pleft]]
          newix1 = sort(intersect(newix1, above.thresh.1))
          newix2 = sort(intersect(newix2, above.thresh.2))
          newix = sort(intersect(newix1, newix2))
        } else {
          length1 = length(which(margcoef1 >= quantile(randcoef1, q)))
          length2 = length(which(margcoef2 >= quantile(randcoef2, q)))
          newix1 = candind[rankcoef1$ix[1:length1]]
          newix2 = candind[rankcoef2$ix[1:length2]]
          iensure = intensure(greedy.size, l1 = newix1, l2 = newix2, k = greedy.size)
          if(is.null(iensure)) newix = NULL
          else newix = sort(intersect(newix1[1:iensure], newix2[1:iensure]))
        }
      } else newix = NULL
    }
  }
  return(newix)
}


intensure <- function(i, l1, l2, k) {
  for(j in i:length(l1)){
    if (length(intersect(l1[1:j], l2[1:j])) >= k)
      return(j)
  }
  # if (length(intersect(l1[1:i], l2[1:i])) >= k)
  #    return(i) else return(intensure(i + 1, l1, l2, k))
}

int.size.k <- function(l1, l2, k) {
  iensure = intensure(k, l1 = l1, l2 = l2, k = k)
  ix01 = l1[1:iensure]
  ix02 = l2[1:iensure]
  ix0 = sort(intersect(ix01, ix02))
  return(ix0)
}

calculate.nsis <- function(family, varISIS, n, p) {
  if (varISIS == "aggr")
    nsis = floor(n/log(n)) else {
      if (family == "gaussian") {
        nsis = floor(n/log(n))
      }
      if (family == "binomial") {
        nsis = floor(n/(4 * log(n)))
      }
      if (family == "poisson") {
        nsis = floor(n/(2 * log(n)))
      }
      if (family == "cox") {
        nsis = floor(n/(4 * log(n)))
      }
    }
  if (p < n)
    nsis = p
  return(nsis)
}
loss.ncvsurv <- function(y, eta, total=TRUE) {
  ind <- order(y[,1])
  d <- as.double(y[ind,2])
  if (is.matrix(eta)) {
    eta <- eta[ind, , drop=FALSE]
    r <- apply(eta, 2, function(x) rev(cumsum(rev(exp(x)))))
  } else {
    eta <- eta[ind]
    r <- rev(cumsum(rev(exp(eta))))
  }
  if (total) {
    return(-2*(crossprod(d, eta) - crossprod(d, log(r))))
  } else {
    return(-2*(eta[d==1,] - log(r)[d==1,]))
  }
}


# # ##### new added#######
# x=Xn
# y=y_surv
# family="cox"
# qn=q
# mn=m
# fn=f
# penalty = method
# method
# tune='bic'
# 
# # object <- tune.fit1(Xn, y_surv, family="cox",q,m,f,penalty = method, tune = "bic")
# 
tune.fit1 <- function(x, y, family = c("gaussian", "binomial", "poisson", "cox"),qn,mn,fn, penalty = c("SCAD", "MCP", "lasso"), concavity.parameter = switch(penalty, SCAD = 3.7, 3), tune = c( "aic", "bic", "ebic"), nfolds = 10,
                      type.measure = c("deviance", "class", "auc", "mse", "mae"), gamma.ebic = 1) {
  
  if (is.null(x) || is.null(y))
    stop("The data is missing!")
  
  this.call = match.call()
  family = match.arg(family)
  penalty = match.arg(penalty)
  if (class(concavity.parameter) != "numeric")
    stop("concavity.parameter must be numeric!")
  tune = match.arg(tune)
  if (class(nfolds) != "numeric")
    stop("nfolds must be numeric!")
  type.measure = match.arg(type.measure)
  n = nrow(x)
  # step 1: the splitting stage
  if (penalty == "lasso" ) {
    reg.fit = glmnet(x, y, family = family)
    coef.beta = rbind(reg.fit$a0,as.matrix(reg.fit$beta))  # extract coefficients at all values of lambda,  including the intercept
    dev = deviance(reg.fit)
    reg.df = reg.fit$df
  } else {
    if(family != 'cox'){
      
      reg.fit = ncvreg(x, y, family = family, penalty = penalty, gamma = concavity.parameter)
      coef.beta = reg.fit$beta  # extract coefficients at all values of lambda, including the intercept
      dev = loglik(x, y, coef.beta, family = family)
      reg.df = getdf(coef.beta[-1, , drop = FALSE])
    } else {
      reg.fit = ncvsurv1(x, y, qn,mn,fn, penalty = penalty, gamma = concavity.parameter)
      coef.beta = reg.fit$beta  # extract coefficients at all values of lambda, including the intercept
      dev = 2*reg.fit$loss
      reg.df = getdf(coef.beta)
    }
    
    # formulated the penalized MLE, with different lambda
    
    
    if (tune == "aic") {
      obj = dev + 2 * reg.df
    }
    if (tune == "bic") {
      obj = dev + log(n) * reg.df
    }
    if (tune == "ebic") {
      obj = dev + log(n) * reg.df + 2 * gamma.ebic * log(choose(dim(x)[2], reg.df))
    }
    lambda.ind = which.min(obj)
    coef.beta = coef.beta[, lambda.ind]
    lambda = reg.fit$lambda[lambda.ind]
  }
  
  if(family != 'cox'){
    a0 = coef.beta[1]
    coef.beta = coef.beta[-1]
  } else{
    a0 = NULL
    coef.beta = as.vector(coef.beta)
  }
  ix = which(coef.beta != 0)
  beta = coef.beta[ix]
  return(list(ix = ix, a0 = a0, beta = beta, fit = reg.fit, lambda = lambda, lambda.ind = lambda.ind))
}

css=function(x1,x2,z,obs.time,delta,ngrid){
  p2=dim(x2)[2]
  fn=function(z.grid){
    z_new=ifelse(z>z.grid,1,0)
    int_new=x1*z_new
    data=data.frame(obs.time,delta,x1,x2,int_new)
    coxest=-coxph(Surv(obs.time,delta)~x1+x2+int_new,data=data,na.action=na.exclude)$loglik[2]
  }
  fn1 =function(z.grid){
    z_new=ifelse(z>z.grid,1,0)
    int_new=x1*z_new
    data=data.frame(obs.time,delta,x1,int_new)
    coxest=-coxph(Surv(obs.time,delta)~x1+int_new,data=data,na.action=na.exclude)$loglik[2]
  }
  low<-quantile(z,0.1)
  high<-quantile(z,0.9)
  if(length(p2)==0){
    ma=gridSearch(fn1,lower=low,upper=high,n=ngrid)$minlevels
  }else{
    ma=gridSearch(fn,lower=low,upper=high,n=ngrid)$minlevels
  }
  return(ma)
}

