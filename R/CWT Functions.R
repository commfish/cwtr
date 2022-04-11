#################
# CWT FUNCTIONS #
#################
#cwtEst
#cwtBoot
#cwtBayes (not coded)
#cwt
#print.cwt
#summary.cwt
#print.summary.cwt
#plot.cwt




#' @title CWT Point Estimation
#'
#' @description This function estimates the point estimate for
#' @param N Number of fish harvested in the fishery, optionally including a vector with a length of 2, of
#' 1) Number of fish harvested in fishery and 2) the associated variance. If entered as a single number, the
#' variance will default to NA.
#' @param n Number of fish sampled in the fishery (inspected for marks)
#' @param lambda Vector of length four containing 1) the number of heads with CWT shipped to the taglab, A1,
#' 2) the number of heads received by the taglab, A2, 3) the number of tags detected at the taglab, M1, and
#' 4) the number of tags successfully decoded, M2.
#' @param m Number of marked fish in the sample
#' @param theta the proportion of stock that is marked, optionally can be a vector of length two
#' containing 1) the proportion of stock that is marked and
#' 2) if applicable, the associated variance of 1/theta, else NA

#' @keywords cwt
#' @export
#' @examples
#' Point estimates for the first example from Geiger (1990)
#' cwtEst(N=c(1000,NA), n=200, lambda=c(1,1,1,1), m=10, theta=c(0.5,NA))
#'
#' Alternately coded as:
#' cwtEst(N=1000, n=200, lambda=c(1,1,1,1), m=10, theta=0.5)

cwtEst <- function(N, n, lambda, m, theta) {
  lambda2 <- (lambda[1]/lambda[2]) * (lambda[3]/lambda[4])
  r.est <- (m*N[1])/(n*theta[1]*lambda2)
  t.est <- (m*N[1])/(n*lambda2)
  q.est <- r.est/N[1]
  #variance
  G.p <- (1 - lambda2*(n/N[1])*theta[1]) / m
  G.N <- ifelse(is.na(N[2]), 0, N[2] / N[1]^2)
  G.theta <- ifelse(is.na(theta[2]), 0, theta[2] / (1/theta[1])^2)
  r.var <- r.est^2 * (G.p + G.N + G.theta - G.p*G.N - G.p*G.theta - G.N*G.theta + G.p*G.N*G.theta)
  G.p2 <- (1 - lambda2*(n/N[1])) / m
  t.var <- t.est^2 * (G.p2 + G.N - G.p2*G.N)
  q.var <- ifelse(is.na(N[2]), r.var/N[1]^2, (r.est/N[1])^2 * (N[2]/N[1]^2 + r.var/r.est^2))
  #store results
  PointEstimates <- data.frame(Estimate=c(r.est,t.est,q.est), Variance=c(r.var,t.var,q.var))
  rownames(PointEstimates) <- c("r.est","t.est","q.est")
  out <- list(
    PointEstimates = PointEstimates,
    InputData = list(N=N,n=n,lambda=lambda,m=m,theta=theta)
  )
  #return results
  class(out) <- "cwt"
  return(out)
}





#' @title Bootstrapped CWT Confidence Intervals
#'
#' @description Create a bootstrap confidence interval
#' @param x default cwt input
#' @param nreps the number of fish sampled in the fishery
#' @param method The method to calculate bootstrapped intervals can be either "parametric" or
#' "geiger". Using parametric bootstrapping follows a deterministic approach; using Geiger's
#' method follows the Geiger (1990) approach to bootstrapping. 

#' @keywords cwt
#' @export
#' @examples
#' x <- cwtEst(N=c(125000,NA), n=35000, lambda=c(800,800,600,600), m=6, theta=c(0.019, 102))
#' xBoot <- cwtBoot(x, method="parametric", nreps=10000)

cwtBoot <- function(x, method = "parametric", nreps = 1000) {
  r <- x$PointEstimates[1,1]
  n <- x$InputData$n
  N <- x$InputData$N
  ncwt <- x$InputData$m
  theta  <- x$InputData$theta
  lambda <- x$InputData$lambda
  lambda2 = (lambda[1]/lambda[2]) * (lambda[3]/lambda[4])
  #
  if(method=="geiger") {
    #P known
    if(is.na(theta[2])) {
      t.boot <- rbinom(nreps, round(r,0), theta[1]) #M_b in Geiger
      m.boot <- rep(NA,nreps)
      for(i in 1:nreps) {
        m.boot[i] <- rhyper(1, t.boot[i], N[1], n) #x_b in Geiger
      }
      r.boot <- m.boot*(theta[1]^-1)*(N[1]/n)
    }
    #P not known
    if(!is.na(theta[2])) {
      theta.boot <- rbinom(nreps, theta[2], theta[1])/theta[2]
      t.boot <- rbinom(nreps, round(r,0), theta.boot) #M_b in Geiger
      m.boot <- rep(NA,nreps)
      for(i in 1:nreps) {
        m.boot[i] <- rhyper(1, t.boot[i], N[1], n) #x_b in Geiger
      }
      r.boot <- ifelse(theta.boot==0, 0, m.boot*(theta.boot^-1)*(N[1]/n))
    }
    q.boot = r.boot/N[1]
  } #end Geiger
  # This is my take on how to do this, which is roughly the same as Geiger's method, but
  #  also noting that how I handle lambda in the below simulation is crude
  if(method=="parametric") {
    N.boot <- rnorm(nreps, N[1], ifelse(is.na(N[2]),0,sqrt(N[2])) )
    N.boot <- ifelse(N.boot < n, n, N.boot)
    N.boot <- round(N.boot, 0)
    theta.boot <- rnorm(nreps, 1/theta[1], ifelse(is.na(theta[2]),0,sqrt(theta[2])) )
    t.boot <- rep(NA,nreps)
    aa.boot <- rep(NA,nreps)
    tt.boot <- rep(NA,nreps)
    m.boot <- rep(NA,nreps)
    for(i in 1:nreps) {
      t.boot[i] <- rbinom(1, round(r,0), 1/theta.boot[i]) #M_b in Geiger
      m.boot[i] <- rhyper(1, t.boot[i], N.boot[i], n) #x_b in Geiger
    }
    aa.boot <- rbinom(nreps, lambda[2], lambda[1]/lambda[2])/lambda[2]
    tt.boot <- rbinom(nreps, lambda[4], lambda[3]/lambda[4])/lambda[4]
    lambda.boot = (1/aa.boot) * (1/tt.boot)
    r.boot <- m.boot*theta.boot*(N.boot/n)*lambda.boot
    q.boot <- r.boot / N.boot
  }
  #store results
  x$Bootstrap <- list(r.boot = r.boot, t.boot = t.boot, m.boot = m.boot, q.boot = q.boot)
  #return output
  class(x) <- c("cwt", "cwtBoot")
  return(x)
}




#' @title cwt
#'
#' @description CWT helper function
#' @param x default cwt input
#' @export

cwt <- function(x) UseMethod("cwt")






#' @title Print the CWT Confidence Interval Summary
#'
#' @description Print cWT output
#' @param x Default cwt input
#' @param digits Number of digits to display, default=3
#' @param quiet Silence verbose notes, default = FALSE

#' @keywords cwt
#' @export
#' @examples
#' print.cwt()

print.cwt <- function(x, digits = 3, quiet = FALSE) {
  if(quiet == TRUE){
    cat("Point Estimates:\n")
    print(round(x$PointEstimates,digits))
    
  } else{
    cat("Point Estimates:\n")
    print(round(x$PointEstimates,digits))
    cat("\nNOTE:\nr(hat) is the # of fish caught\nt(hat) is the # of tagged fish caught\nq(hat) is the % of fish caught\n")
    
  }
  
    if(inherits(x, "cwtBoot")) {
    cat("\nBootstrap Estimates:\n")
    y <- x$Bootstrap[1:3]
    outBoot <- cbind(lapply(y, mean),
                     lapply(y, sd))
    colnames(outBoot) <- c("Mean", "SE")
    print(outBoot)
  }
}





#' @title CWT summary
#'
#' @description Display a summary of the
#' @param x default cwt input
#' @param alpha Alpha level for confidence, default=0.1
#' @param quiet Silence verbose notes, default = FALSE
#' @keywords cwt
#' @export
#' @examples
#' summary.cwt()

summary.cwt <- function(x, alpha = 0.1, quiet = FALSE, ...) {
  out <- unclass(x)
  outPoint <- out$PointEstimates
  outPoint$SE <- sqrt(outPoint$Variance)
  outPoint$CV <- outPoint$SE / outPoint$Estimate
  outPoint$CI_lwr <- outPoint$Estimate - qnorm(1 - alpha / 2, 0, 1) * outPoint$SE
  outPoint$CI_upp <- outPoint$Estimate + qnorm(1 - alpha / 2, 0, 1) * outPoint$SE
  out$SummaryPoint <- outPoint
  class(out) <- "summary.cwt"
  if(inherits(x, "cwtBoot")) {
    y <- x$Bootstrap[1:3]
    outBoot <- cbind(lapply(y, mean),
                     lapply(y, var),
                     lapply(y, sd),
                     unlist(lapply(y, sd))/unlist(lapply(y, mean)),
                     unlist(lapply(y, mean)) - x$PointEstimates$Estimate,
                     lapply(y, quantile, alpha/2),
                     lapply(y, quantile, 1-alpha/2))
    colnames(outBoot) <- c("Mean", "Variance", "SE", "CV", "Bias", "CI_lwr", "CI_uppr")
    out$SummaryBoot <- outBoot
    class(out) <- c("summary.cwt", "summary.cwtBoot")
  }
  return(out)
}





#' @title Print CWT summary
#'
#' @description Display a summary of the CWT estimates
#' @param x Default cwt input
#' @param digits Number of digits to display, default=3
#' @param quiet Silence verbose notes, default = FALSE
#' @keywords cwt
#' @export
#' @examples
#' print.summary.cwt()

print.summary.cwt <- function(x, digits = 3, quiet = FALSE, ...) {
  cat("Point Estimates Summary:\n")
  print(round(x$SummaryPoint,digits))
  if(inherits(x, "summary.cwtBoot")) {
    cat("\nBootstrap Estimates Summary:\n")
    print(x$SummaryBoot)
  }
  if(quiet == TRUE){
    cat("")
  } else{
    cat("\nNOTE:\nr is the # of fish caught\nt is the # of tagged fish caught\nq is the % of fish caught\n")
  }
}






#' @title Plot CWT class data
#'
#' @description Four panel plot summary
#' @param x default cwt input
#' @param parameters Vector of length 4 of "r", "t", "q", and "m"
#' @param type Plot type
#' @keywords cwt
#' @export
#' @examples
#' plot.cwt()
plot.cwt <- function(x, parameters = c("r", "t", "q", "m"), type = "histogram") {
  #NOTE, plot method is available for objects that inherit cwtBoot and cwtBayes only
  if(inherits(x, "cwtBoot")) {
    par(mfrow=c(2,2))
    if(type == "histogram") out <- lapply(x$Bootstrap, hist, plot=FALSE)
    else if (type == "density") out <- lapply(x$Bootstrap, density)
    plot(out$r.boot, main="Bootstrap Distribution of r", xlab="r")
    plot(out$t.boot, main="Bootstrap Distribution of t", xlab="t")
    plot(out$q.boot, main="Bootstrap Distribution of q", xlab="q")
    plot(out$m.boot, main="Bootstrap Distribution of m", xlab="m")
  }
}

