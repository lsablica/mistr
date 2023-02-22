#' @title Fitting a GPD-Normal-GPD Model
#' @description \code{GNG_fit} is used to fit three components composite models with components GPD, normal and GPD.
#' @param data vector of values to which the density is optimized.
#' @param start named vector (break1, break2, mean, sd, shape1, shape2) of values that are used to start the optimization,
#'              default: c(break1 = -0.02, break2 = 0.02, mean = 0, sd = 0.0115, shape1 = 0.15,shape2 = 0.15).
#' @param break_fix logical, fix the breakpoints at the values from start?, default: FALSE.
#' @param midd split reals into two subintervals, the first breakpoint is then optimized on the left of \code{midd} and the second on the right, default: mean(data).
#' @param ... further arguments to be passed to the optimizer.
#' @return  A list of class comp_fit.
#' @details The GNG model is the GPD-Normal-GPD model. This means
#'          that a \eqn{-X} transformation of a GPD random variable will be used for the left tail,
#'          normal distribution for the center and again GPD for the right tail.
#'
#'          The code uses the maximum likelihood estimation technique to estimate the six parameters from the start vector
#'          (\code{break1, break2, mean, sd, shape1, shape2}). The other parameters (location and scale parameters of the GPD) are
#'          computed in each step such that the function is continuous. Weights are estimated in every step as a proportion
#'          of points that correspond to each of the truncated region. If the breakpoints are fixed (i.e. \code{break_fix = TRUE}),
#'           the weights are computed before the optimization procedure.
#'
#'          Optimization is handled by the \code{\link[bbmle]{mle2}} function.
#' @examples
#' \dontrun{
#'  GNG_fit(stocks$SAP)
#'
#'  GNG_fit(stocks$MSFT)
#'
#'  autoplot(GNG_fit(stocks$ADS))
#'
#'  GNG_fit(stocks$GSPC, start = c(break1=-0.0075, break2=0.0075, mean=0,
#'          sd=0.0115, shape1=0.15, shape2=0.15), control = list(maxit = 20000))
#'
#'  GNG_fit(stocks$DJI, start = c(break1=-0.0055, break2=0.0055, mean=-0.001,
#'          sd=0.0055,shape1=0.15,shape2=0.15), method = "CG",control = list(maxit = 1000))
#'
#' }
#' @seealso
#'  \code{\link[bbmle]{mle2}}
#' @rdname GNG_fit
#' @export
#' @importFrom bbmle mle2
GNG_fit <- function(data, start = c(break1 = -0.02, break2 = 0.02, mean = 0, sd = 0.0115, shape1 = 0.15, shape2 = 0.15),
                    break_fix = FALSE, midd = mean(data), ...) {
   data <- as.vector(na.omit(data))
   if (length(start) != 6) stop("you need to specify all 6 starting values")
   if (is.null(names(start))) names(start) <- c("break1", "break2", "mean", "sd", "shape1", "shape2")
   if (!all(c("break1", "break2", "mean", "sd", "shape1", "shape2") %in% names(start))) stop("you need to specify names with 'break1','break2','mean','sd','shape1','shape2'")
   if (diff(start[c("break1", "break2")]) <= 0) stop("break1 has to be smaller than break2")
   if (break_fix) {
      GNG_loglike_fix <- function(x, break1, break2, weights, shape1, mean, sd, shape2) {
         sd <- exp(sd) + 1e-10
         scale1 <- (weights[1] * diff(pnorm(breakpoints, mean, sd)))/(weights[2] * dnorm(breakpoints[1], mean, sd))
         scale2 <- (weights[3] * diff(pnorm(breakpoints, mean, sd)))/(weights[2] * dnorm(breakpoints[2], mean, sd))
         if (any(is.na(c(scale1, scale2))) || any(c(scale1, scale2) == 0) || any(is.infinite(c(scale1, scale2)))) return(1e16)
         
         sort <- findInterval(x, c(break1, break2))
         a1 <- sum(dGPD(-x[sort == 0], -break1, scale1, shape1, log = TRUE) + log(weights[1]))
         a2 <- sum(dnorm(x[sort == 1], mean, sd, log = TRUE) + log(weights[2]) - log(diff(pnorm(c(break1,break2), mean, sd))))
         a3 <- sum(dGPD(x[sort == 2], break2, scale2, shape2, log = TRUE) + log(weights[3]))
         h <- -(a1+a2+a3)
         if(is.infinite(h)) return(1e16)
         h
      }
      breakpoints <- start[c("break1", "break2")]
      t <- findInterval(data, breakpoints) + 1
      w <- unlist(lapply(1:3, function(x) sum(t == x)), use.names = FALSE)/length(data)
      dat  <- list(x = data, break1 = start["break1"], break2 = start["break2"], weights = w)
      star <- list(shape1 = start["shape1"], mean = start["mean"],
                    sd = log(start["sd"] - 1e-10), shape2 = start["shape2"])
      fit <- bbmle::mle2(GNG_loglike_fix, start = star, data = dat, ...)
      lik <- -GNG_loglike_fix(data, breakpoints[1], breakpoints[2], w, fit@coef["shape1"], fit@coef["mean"],
                                                                        fit@coef["sd"], fit@coef["shape2"])
      coe <- fit@coef
   } else {
      GNG_loglike <- function(x, break1, break2, shape1, mean, sd, shape2) {
         breakpoints <- c(-exp(break1) - 1e-10 + midd, exp(break2) + 1e-10 + midd)
         sort <- findInterval(x, breakpoints)
         weights <- unlist(lapply(0:2, function(i) sum(sort == i)), use.names = FALSE)/length(x)
         sd <- exp(sd) + 1e-10
         scale1 <- (weights[1] * diff(pnorm(breakpoints, mean, sd)))/(weights[2] * dnorm(breakpoints[1], mean, sd))
         scale2 <- (weights[3] * diff(pnorm(breakpoints, mean, sd)))/(weights[2] * dnorm(breakpoints[2], mean, sd))
         if (any(is.nan(c(scale1, scale2))) || any(c(scale1, scale2) == 0) || any(is.infinite(c(scale1, scale2)))) return(1e16)
         
         a1 <- sum(dGPD(-x[sort == 0], -breakpoints[1], scale1, shape1, log = TRUE) + log(weights[1]))
         a2 <- sum(dnorm(x[sort == 1], mean, sd, log = TRUE) + log(weights[2]) - log(diff(pnorm(breakpoints, mean, sd))))
         a3 <- sum(dGPD(x[sort == 2], breakpoints[2], scale2, shape2, log = TRUE) + log(weights[3]))
         h <- -(a1+a2+a3)
         if(is.infinite(h) || is.nan(h)) return(1e16)
         h
      }
      star <- list(break1 = log(-start["break1"] - 1e-10 + midd), break2 = log(start["break2"] - 1e-10 - midd),
                    shape1 = start["shape1"], mean = start["mean"], sd = log(start["sd"] - 1e-10), shape2 = start["shape2"])
      fit <- bbmle::mle2(GNG_loglike, start = star, data = list(x = data), ...)
      coe <- fit@coef[!(names(start) %in% c("break1", "break2"))]
      breakpoints <- c(-exp(fit@coef["break1"]) - 1e-10 + midd, exp(fit@coef["break2"]) + 1e-10 + midd)
      t <- findInterval(data, breakpoints) + 1
      w <- unlist(lapply(1:3, function(x) sum(t == x)), use.names = FALSE)/length(data)
      lik <- -GNG_loglike(data, fit@coef["break1"], fit@coef["break2"], fit@coef["shape1"], fit@coef["mean"],
                              fit@coef["sd"], fit@coef["shape2"])
   }
   mean <- coe["mean"]
   sd <- exp(coe["sd"]) + 1e-10
   scale1 <- (w[1] * diff(pnorm(breakpoints, mean, sd)))/(w[2] * dnorm(breakpoints[1], mean, sd))
   scale2 <- (w[3] * diff(pnorm(breakpoints, mean, sd)))/(w[2] * dnorm(breakpoints[2], mean, sd))
   coe <- c(-breakpoints[1], scale1, coe["shape1"], mean, sd, breakpoints[2], scale2, coe["shape2"])
   names(coe) <- c("loc1", "scale1", "shape1", "mean", "sd", "loc2", "scale2", "shape2")
   dist <- compdist(-GPDdist(-breakpoints[1], scale1, coe["shape1"]), normdist(mean, sd), 
                    GPDdist(breakpoints[2], scale2, coe["shape2"]), weights = w, breakpoints = breakpoints)
   r <- list(Distribution = dist, convergence = list(convergence = fit@details$convergence, message = fit@details$message),
             fit = fit, params = list(breakpoints = breakpoints, weights = w, coeff = coe), 
             spec = list(name = "GPD-Normal-GPD", lik = lik), data = data)
   class(r) <- c("GNG", "comp_fit")
   r
}


#' @title Fitting a Pareto-Normal-Pareto Model
#' @description \code{GNG_fit} is used to fit three components composite models with components Pareto, normal and Pareto.
#' @param data vector of values to which the density is optimized.
#' @param start named vector (break1, break2, mean, sd) of values that are used to start the optimization,
#'              default: c(break1 = -0.02, break2 = 0.02, mean = 0, sd = 0.012).
#' @param ... further arguments to be passed to optimizer.
#' @return  A list of class comp_fit.
#' @details The PNP model is the Pareto-Normal-Pareto model. This means
#'          that a \eqn{-X} transformation of a Pareto random variable will be used for the left tail,
#'          normal distribution for the center and again Pareto for the right tail.
#'
#'          The code uses the maximum likelihood estimation technique to estimate the four parameters from the start vector
#'          (\code{break1, break2, mean, sd}). The other parameters (shape parameters of Pareto distribution) are
#'          computed in each step such that the function is continuous. Weights are estimated in every step as a proportion
#'          of points that correspond to each of the truncated region.
#'
#'          Optimization is handled by the \code{\link[bbmle]{mle2}} function.
#' @examples
#' \dontrun{
#'  PNP_fit(stocks$SAP)
#'
#'  PNP_fit(stocks$MSFT)
#'
#'  autoplot(PNP_fit(stocks$ADS))
#'
#'  PNP_fit(stocks$GSPC, method = "BFGS")
#'
#'  PNP_fit(stocks$DJI, start = c(-0.01,0.01,0,0.008))
#'
#' }
#' @seealso
#'  \code{\link[bbmle]{mle2}}
#' @rdname PNP_fit
#' @export
#' @importFrom bbmle mle2
PNP_fit <- function(data, start = c(break1 = -0.02, break2 = 0.02, mean = 0, sd = 0.012), ...) {
   data <- as.vector(na.omit(data))
   if (length(start) != 4) stop("you must set all 4 starting values")
   if (is.null(names(start))) names(start) <- c("break1", "break2", "mean", "sd")
   if (!all(c("break1", "break2", "mean", "sd") %in% names(start))) stop("you must specify names with 'break1','break2','mean','sd'")
   
   PNP_loglike <- function(x, break1, break2, mean, sd) {
      sd <- exp(sd) + 1e-10
      breakpoints <- c(-exp(break1) - 1e-10,exp(break2) + 1e-10)
      sort <- findInterval(x, breakpoints)
      weights <- unlist(lapply(0:2, function(i) sum(sort == i)), use.names = FALSE)/length(x)
      shape1 <- -breakpoints[1] * (weights[2] * dnorm(breakpoints[1], mean, sd))/(weights[1] * diff(pnorm(breakpoints, mean, sd)))
      shape2 <- breakpoints[2] * (weights[2] * dnorm(breakpoints[2], mean, sd))/(weights[3] * diff(pnorm(breakpoints, mean, sd)))
      if (any(is.na(c(shape1, shape2))) || any(c(shape1, shape2) == 0) || any(is.infinite(c(shape1, shape2)))) return(1e16)
      Z <- numeric(length(x))
      Z[sort == 0] <- dpareto(-x[sort == 0], -breakpoints[1], shape1, log = TRUE) + log(weights[1])
      Z[sort == 1] <- dnorm(x[sort == 1], mean, sd, log = TRUE) + log(weights[2]) - log(diff(pnorm(breakpoints, mean, sd)))
      Z[sort == 2] <- dpareto(x[sort == 2], breakpoints[2], shape2, log = TRUE) + log(weights[3])
      -sum(Z)
   }
   fit <- bbmle::mle2(PNP_loglike, start = list(break1 = log(-start["break1"] - 1e-10),
                                                break2 = log(start["break2"] - 1e-10),
                                                mean = start["mean"], sd = log(start["sd"] - 1e-10)),
                                   data = list(x = data), ...)
   
   breakpoints <- c(-exp(fit@coef["break1"]) - 1e-10, exp(fit@coef["break2"] + 1e-10))
   t <- findInterval(data, breakpoints) + 1
   w <- unlist(lapply(1:3, function(x) sum(t == x)), use.names = FALSE)/length(data)
   mean = fit@coef["mean"]
   sd = exp(fit@coef["sd"]) + 1e-10
   scale1 <- -breakpoints["break1"]
   scale2 <- breakpoints["break2"]
   shape1 <- -breakpoints[1] * (w[2] * dnorm(breakpoints[1], mean, sd))/(w[1] * diff(pnorm(breakpoints, mean, sd)))
   shape2 <- breakpoints[2] * (w[2] * dnorm(breakpoints[2], mean, sd))/(w[3] * diff(pnorm(breakpoints, mean, sd)))
   coe <- c(scale1, shape1, mean, sd, scale2, shape2)
   names(coe) <- c("scale1", "shape1", "mean", "sd", "scale2", "shape2")
   lik <- -PNP_loglike(data, fit@coef["break1"], fit@coef["break2"], fit@coef["mean"], fit@coef["sd"])
   dist <- compdist(-paretodist(scale1, shape1), normdist(mean, sd), paretodist(scale2, shape2), weights = w, breakpoints = breakpoints)
   r <- list(Distribution = dist, convergence = list(convergence = fit@details$convergence, message = fit@details$message),
             fit = fit, params = list(breakpoints = breakpoints, weights = w, coeff = coe), spec = list(name = "Pareto-Normal-Pareto", lik = lik), data = data)
   class(r) <- c("PNP", "comp_fit")
   r
}

#' @export
print.comp_fit <- function(x, digits = 6, ...) {
  cat("Fitted composite", x$spec$name, "distribution: \n\n")
  cat("Breakpoints:", round(x$params$breakpoints, digits), "\n")
  cat("Weights:", round(x$params$weights, digits), "\n\n")
  cat("Parameters: \n")
  print(round(x$params$coef, digits))
  cat("\nLog-likelihood: ", x$spec$lik, ",  Average log-likelihood: ", round(x$spec$lik/length(x$data), 4), "\n\n", sep = "")
  invisible(x)
}

int_mGPD <- function(alpha, loc, scale, shape, w1, rho1) {
  unname(alpha * (-loc + scale/shape - scale * (alpha * rho1/w1)^(-shape)/(shape * (1 - shape))))
}
int_mP <- function(alpha, scale, shape, w1, rho1) {
  unname(-scale * shape * alpha^(1 - 1/shape)/((rho1/w1)^(1/shape) * (shape - 1)))
}
int_N <- function(alpha, mean, sd, w1, w2, rho2, lam2) {
  r1 <- lam2 + ((alpha - w1)/w2) * (rho2 - lam2)
  unname((w2/(rho2 - lam2)) * (mean * (r1 - lam2) - sd * dnorm(qnorm(r1)) + sd * dnorm(qnorm(lam2))))
}
int_GPD <- function(alpha, loc, scale, shape, w1, w2, w3, lam3) {
  r1 <- lam3 + (alpha - w1 - w2) * (1 - lam3)/w3
  r2 <- numeric(length(alpha)) + lam3
  r1[alpha == 1] <- 1
  r2[alpha == 1] <- 0
  unname((w3/(1 - lam3)) * ((r1 - r2) * (loc - scale/shape) - scale * (1 - r1)^(1 - shape)/(shape * (1 - shape)) + scale *
                              (1 - r2)^(1 - shape)/(shape * (1 - shape))))
}
int_P <- function(alpha, scale, shape, w1, w2, w3, lam3) {
  r1 <- lam3 + (alpha - w1 - w2) * (1 - lam3)/w3
  r2 <- numeric(length(alpha)) + lam3
  r1[alpha == 1] <- 1
  r2[alpha == 1] <- 0
  unname((w3/(1 - lam3)) * scale * shape * ((1 - r2)^(1 - 1/shape) - (1 - r1)^(1 - 1/shape))/(shape - 1))
}
int_PNP <- function(model, alpha) {
  m <- matrix(numeric(3 * length(alpha)), ncol = 3)
  cs <- cumsum(model$params$weights)
  index <- findInterval(alpha, cs, left.open = TRUE) + 1
  m[index == 1, 1] <- alpha[index == 1]
  m[index != 1, 1] <- cs[1]
  m[index == 2, 2] <- alpha[index == 2]
  m[index == 3, 2] <- cs[2]
  m[index == 3, 3] <- alpha[index == 3]
  m[, 1] <- int_mP(alpha = m[, 1], scale = parameters(model)["scale1"], shape = parameters(model)["shape1"], w1 = weights(model)[1],
                   rho1 = distribution(model)$trunc$rho[1])
  if (any(m[, 2] > 0)) {
    m[, 2][m[, 2] != 0] <- int_N(alpha = m[, 2][m[, 2] != 0], mean = parameters(model)["mean"], sd = parameters(model)["sd"],
                                 w1 = weights(model)[1], w2 = weights(model)[2], rho2 = distribution(model)$trunc$rho[2], lam2 = distribution(model)$trunc$lambda[2])
    if (any(m[, 3] > 0)) {
      m[, 3][m[, 3] != 0] <- int_P(alpha = m[, 3][m[, 3] != 0], scale = parameters(model)["scale2"], shape = parameters(model)["shape2"],
                                   w1 = weights(model)[1], w2 = weights(model)[2], w3 = weights(model)[3], lam3 = distribution(model)$trunc$lambda[3])
    }
  }
  rowSums(m)
}
int_GNG <- function(model, alpha) {
  m <- matrix(numeric(3 * length(alpha)), ncol = 3)
  cs <- cumsum(model$params$weights)
  index <- findInterval(alpha, cs, left.open = TRUE) + 1
  m[index == 1, 1] <- alpha[index == 1]
  m[index != 1, 1] <- cs[1]
  m[index == 2, 2] <- alpha[index == 2]
  m[index == 3, 2] <- cs[2]
  m[index == 3, 3] <- alpha[index == 3]
  m[, 1] <- int_mGPD(alpha = m[, 1], loc = parameters(model)["loc1"], scale = parameters(model)["scale1"], shape = parameters(model)["shape1"],
                     w1 = weights(model)[1], rho1 = distribution(model)$trunc$rho[1])
  if (any(m[, 2] > 0)) {
    m[, 2][m[, 2] != 0] <- int_N(alpha = m[, 2][m[, 2] != 0], mean = parameters(model)["mean"], sd = parameters(model)["sd"],
                                 w1 = weights(model)[1], w2 = weights(model)[2], rho2 = distribution(model)$trunc$rho[2], lam2 = distribution(model)$trunc$lambda[2])
    if (any(m[, 3] > 0)) {
      m[, 3][m[, 3] != 0] <- int_GPD(alpha = m[, 3][m[, 3] != 0], loc = parameters(model)["loc2"], scale = parameters(model)["scale2"],
                                     shape = parameters(model)["shape2"], w1 = weights(model)[1], w2 = weights(model)[2], w3 = weights(model)[3],
                                     lam3 = distribution(model)$trunc$lambda[3])
    }
  }
  rowSums(m)
}

expectile_PNP <- function(model, alpha) {
  v <- q(distribution(model), alpha)
  r <- int_PNP(model, alpha) - alpha * v
  E <- int_PNP(model, 1)
  r/(2 * r + v - E)
}

expectile_GNG <- function(model, alpha) {
  v <- q(distribution(model), alpha)
  r <- int_GNG(model, alpha) - alpha * v
  E <- int_GNG(model, 1)
  r/(2 * r + v - E)
}
#' @title Risk Measures of Fitted Objects
#' @description \code{risk} computes the VaR, ES and expectiles at a given level for fitted distribution.
#' @param model output object of \code{GNG_fit()} or \code{PNP_fit()}.
#' @param alpha levels of risk measures.
#' @param expectile logical, if also expectiles should be computed, default: TRUE.
#' @param plot plot the results?, default: FALSE.
#' @param ggplot plot the results with ggplot2?, default: FALSE.
#' @param text_ylim y coordinate for annotation in ggplot2, default: -0.15.
#' @param size size of the text indicating the risk measures in the plot, default: 1.
#' @return  List of class risk_measures.
#' @details VaR are computed using the \code{q()} call of the fitted distribution.
#'
#'          ES is computed directly (i.e. the integrals are precomputed, not numerically)
#'          as an integral of the quantile function.
#'
#'          Expectiles can be obtained as a unit-root solution of the identity between quantiles
#'          and expectiles. These are equivalent for corresponding \eqn{\tau} and \eqn{\alpha}
#'          if \deqn{\tau=(\alpha q(\alpha) -G(\alpha))/(\mu - 2G(\alpha)-(1-2\alpha)q(\alpha))} where \eqn{\mu}
#'          is mean, \eqn{q()} is the quantile function and \eqn{G(\alpha) =\int_{-\infty}^{q(\alpha)} y dF(y)}.
#' @examples
#' \dontrun{
#'  GNG <- GNG_fit(stocks$SAP)
#'  PNP <- PNP_fit(stocks$MSFT)
#'
#'  risk(PNP, alpha = c(0.01,0.05,0.08,0.1))
#'  risk(GNG, alpha = c(0.01,0.05,0.08,0.1), plot = TRUE)
#' }
#' @rdname risk
#' @export
risk <- function(model, alpha, expectile = TRUE, plot = FALSE, ggplot = FALSE, text_ylim = -0.15, size = 1) UseMethod("risk")
#' @rdname risk
#' @export
risk.PNP <- function(model, alpha = 0.05, expectile = TRUE, plot = FALSE, ggplot = FALSE, text_ylim = -0.15, size = 1) {
  if (any(alpha > 100) || any(alpha <= 0))
    stop("alpha out of bounds")
  if (any(alpha > 1))
    alpha <- alpha/100
  if (model$params$coeff["shape1"] <=1 || model$params$coeff["shape2"] <=1) 
    warning("Not all components have finite expectation. ES and Exp can be infinite.")
  
  ES <- int_PNP(model, alpha)/alpha
  VaR <- q(distribution(model), alpha)
  da <- data.frame(level = alpha, VaR = -VaR, ES = -ES)
  
  if (expectile) {
    qq <- sapply(alpha, function(tau) {
      uniroot(function(x) {
        expectile_PNP(model, x) - tau
      }, interval = c(0 + 1e-10, 1 - 1e-10))$root
    })
    da <- cbind(da, Exp = -q(distribution(model), qq))
  }
  
  p <- NULL
  if (plot == TRUE) {
    plot_risk(model, da, size)
  }
  if (ggplot == TRUE) {
    p <- plot_riskgg(model, da, text_ylim, size)
    print(p)
  }
  structure(list(results = da, plot = p), class = "risk_measures")
}
#' @rdname risk
#' @export
risk.GNG <- function(model, alpha = 0.05, expectile = TRUE, plot = FALSE, ggplot = FALSE, text_ylim = -0.15, size = 1) {
  if (any(alpha > 100) || any(alpha <= 0))
    stop("alpha out of bounds")
  if (any(alpha > 1))
    alpha <- alpha/100
  if (model$params$coeff["shape1"] >=1 || model$params$coeff["shape2"] >=1) 
    warning("Not all components have finite expectation. ES and Exp can be infinite.")
  
  ES <- int_GNG(model, alpha)/alpha
  VaR <- q(distribution(model), alpha)
  da <- data.frame(level = alpha, VaR = -VaR, ES = -ES)
  
  if (expectile) {
    qq <- sapply(alpha, function(tau) {
      uniroot(function(x) {
        expectile_GNG(model, x) - tau
      }, interval = c(0 + 1e-10, 1 - 1e-10))$root
    })
    da <- cbind(da, Exp = -q(distribution(model), qq))
  }
  
  p <- NULL
  if (plot == TRUE) {
    plot_risk(model, da, size)
  }
  if (ggplot == TRUE) {
    p <- plot_riskgg(model, da, text_ylim, size)
    print(p)
  }
  structure(list(results = da, plot = p), class = "risk_measures")
}
#' @export
print.risk_measures <- function(x,...){
  print(x$results)
  invisible(x)
}
