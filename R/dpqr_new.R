#' @title The Burr Distribution
#' @description Density, distribution function, quantile function and random generation for the Burr distribution with parameters
#' shape1 and shape2.
#' @param x,q  vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param shape1 shape parameter.
#' @param shape2 shape parameter.
#' @param log,log.p logical; if TRUE, probabilities \eqn{p} are given as \eqn{log(p)}, default: FALSE.
#' @param lower.tail logical; if TRUE, probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}, default: TRUE.
#' @details The Burr distribution function with shape1 parameter c and shape2 parameter k has density given by
#'          \deqn{f(x)=ckx^(c-1)/(1+x^c)^(k+1)} for \eqn{x>0}. The cumulative distribution function is
#'          \deqn{F(x)=1-(1+x^c)^-k} on \eqn{x>0}.
#'
#'          See \url{https://en.wikipedia.org/wiki/Burr_distribution} for more details.
#' @return \code{dburr} gives the density, \code{pburr} gives the distribution function, \code{qburr} gives the quantile function,
#'          and \code{rburr} generates random deviates.
#'
#'          Invalid arguments will result in return value NaN, with a warning.
#' @examples
#' dburr(seq(1, 5), 2, 2)
#' qburr(pburr(seq(1, 5), 2, 2), 2 ,2)
#' rburr(5, 2, 2)
#' @rdname Burr
#' @name Burr
#' @seealso \code{\link{burrdist}}
#' @export
dburr <- function(x, shape1, shape2, log = FALSE) {
  if (shape1 <= 0 || shape2 <= 0) {
    warning("NaNs produced")
    return(rep.int(NaN, length(x)))
  }
  Z <- numeric(length(x))
  nat <- is.na(x)
  Z[nat] <- x[nat]
  x <- x[!nat]
  z <- numeric(length(x))
  z[x > 0] <- shape1 * shape2 * x[x > 0]^(shape1 - 1)/(1 + x[x > 0]^shape1)^(shape2 + 1)
  Z[!nat] <- z
  if (log)
    log(Z) else Z
}

#' @rdname Burr
#' @export
pburr <- function(q, shape1, shape2, lower.tail = TRUE, log.p = FALSE) {
    if (shape1 <= 0 || shape2 <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, length(q)))
    }
    Z <- numeric(length(q))
    nat <- is.na(q)
    Z[nat] <- q[nat]
    q <- q[!nat]
    z <- numeric(length(q))
    z[q > 0] <- 1 - (1 + q[q > 0]^shape1)^(-shape2)
    Z[!nat] <- z
    if (!lower.tail) {
        Z <- 1 - Z
    }
    if (log.p)
        log(Z) else Z
}



#' @rdname Burr
#' @export
qburr <- function(p, shape1, shape2, lower.tail = TRUE, log.p = FALSE) {
    if (shape1 <= 0 || shape2 <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, length(p)))
    }
    Z <- numeric(length(p))
    nat <- is.na(p)
    Z[nat] <- p[nat]
    p <- p[!nat]
    z <- numeric(length(p))
    if (log.p)
        p <- exp(p)
    if (!lower.tail)
        p <- 1 - p
    ok <- p >= 0 & p <= 1
    p <- p[ok]
    z[ok] <- ((1 - p)^(-1/shape2) - 1)^(1/shape1)
    z[!ok] <- NaN
    Z[!nat] <- z
    Z
}

#' @rdname Burr
#' @export
rburr <- function(n, shape1, shape2) {
    if (shape1 <= 0 || shape2 <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, n))
    }
    ((runif(n))^(-1/shape2) - 1)^(1/shape1)
}



#' @title The Gumbel Distribution
#' @description Density, distribution function, quantile function and random generation for the Gumbel distribution with
#' location and scale parameters.
#' @param x,q  vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param loc location parameter.
#' @param scale scale parameter.
#' @param log,log.p logical; if TRUE, probabilities \eqn{p} are given as \eqn{log(p)}, default: FALSE.
#' @param lower.tail logical; if TRUE, probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}, default: TRUE.
#' @details The Gumbel distribution function with location parameter \eqn{\mu} and scale parameter \eqn{\beta} has density given by
#'          \deqn{f(x)=1/\beta e^-(z+e^-z)}, where \eqn{z=(x-\mu)/\beta}. The cumulative distribution function is
#'          \deqn{F(x)=e^(-e^z)} with \eqn{z} as stated above.
#'
#'          See \url{https://en.wikipedia.org/wiki/Gumbel_distribution} for more details.
#' @return \code{dgumbel} gives the density, \code{pgumbel} gives the distribution function, \code{qgumbel} gives the quantile function, and
#'         \code{rgumbel} generates random deviates.
#'
#'         Invalid arguments will result in return value NaN, with a warning.
#' @examples
#' dgumbel(seq(1, 5), 0, 1)
#' qgumbel(pgumbel(seq(1, 5), 0, 1), 0 ,1)
#' rgumbel(5, 0, 1)
#' @rdname Gumbel
#' @name Gumbel
#' @seealso \code{\link{gumbeldist}}
#' @export
dgumbel <- function(x, loc, scale, log = FALSE) {
    if (scale <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, length(x)))
    }
    Z <- numeric(length(x))
    nat <- is.na(x)
    Z[nat] <- x[nat]
    x <- x[!nat]
    z <- exp(-((x - loc)/scale + exp(-(x - loc)/scale)))/scale
    Z[!nat] <- z
    if (log)
        log(Z) else Z
}

#' @rdname Gumbel
#' @export
pgumbel <- function(q, loc, scale, lower.tail = TRUE, log.p = FALSE) {
    if (scale <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, length(q)))
    }
    Z <- numeric(length(q))
    nat <- is.na(q)
    Z[nat] <- q[nat]
    q <- q[!nat]
    z <- exp(-exp(-(q - loc)/scale))
    Z[!nat] <- z
    if (!lower.tail) {
        Z <- 1 - Z
    }
    if (log.p)
        log(Z) else Z
}

#' @rdname Gumbel
#' @export
qgumbel <- function(p, loc, scale, lower.tail = TRUE, log.p = FALSE) {
    if (scale <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, length(p)))
    }
    Z <- numeric(length(p))
    nat <- is.na(p)
    Z[nat] <- p[nat]
    p <- p[!nat]
    z <- numeric(length(p))
    if (log.p)
        p <- exp(p)
    if (!lower.tail)
        p <- 1 - p
    ok <- p >= 0 & p <= 1
    p <- p[ok]
    z[ok] <- loc - scale * log(-log(p))
    z[!ok] <- NaN
    Z[!nat] <- z
    Z
}

#' @rdname Gumbel
#' @export
rgumbel <- function(n, loc, scale) {
    if (scale <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, n))
    }
    loc - scale * log(-log(runif(n)))
}



#' @title The Frechet Distribution
#' @description Density, distribution function, quantile function and random generation for the Frechet distribution with
#' location, scale and shape parameters.
#' @param x,q  vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param loc location parameter.
#' @param scale scale parameter.
#' @param shape shape parameter.
#' @param log,log.p logical; if TRUE, probabilities \eqn{p} are given as \eqn{log(p)}, default: FALSE.
#' @param lower.tail logical; if TRUE, probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}, default: TRUE.
#' @details The Frechet distribution function with location parameter \eqn{m}, scale parameter \eqn{s} and shape parameter
#'          \eqn{\alpha} has density given by
#'          \deqn{f(x)=\alpha/s z^(-\alpha-1) e^-z^-\alpha} for \eqn{x>m}, where \eqn{z=(x-m)/s}. The cumulative distribution function is
#'          \deqn{F(x)=e^-z^-\alpha} for \eqn{x>m}, with \eqn{z} as stated above.
#'
#'          See \url{https://en.wikipedia.org/wiki/Frechet_distribution} for more details.
#' @return \code{dfrechet} gives the density, \code{pfrechet} gives the distribution function, \code{qfrechet} gives the quantile function, and
#'         \code{rfrechet} generates random deviates.
#'
#'         Invalid arguments will result in return value NaN, with a warning.
#' @examples
#' dfrechet(seq(1, 5), 0, 1, 1)
#' qfrechet(pfrechet(seq(1, 5), 0, 1, 1), 0, 1, 1)
#' rfrechet(5, 0, 1, 1)
#' @rdname Frechet
#' @name Frechet
#' @seealso \code{\link{frechetdist}}
#' @export
dfrechet <- function(x, loc = 0, scale = 1, shape = 1, log = FALSE) {
    if (scale <= 0 || shape <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, length(x)))
    }
    Z <- numeric(length(x))
    nat <- is.na(x)
    Z[nat] <- x[nat]
    x <- x[!nat]
    z <- numeric(length(x))
    zz <- (x - loc)/scale
    t <- x > loc
    z[t] <- shape/scale * zz[t]^(-1 - shape) * exp(-zz[t]^(-shape))
    Z[!nat] <- z
    if (log)
        log(Z) else Z
}

#' @rdname Frechet
#' @export
pfrechet <- function(q, loc = 0, scale = 1, shape = 1, lower.tail = TRUE, log.p = FALSE) {
    if (scale <= 0 || shape <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, length(q)))
    }
    Z <- numeric(length(q))
    nat <- is.na(q)
    Z[nat] <- q[nat]
    q <- q[!nat]
    z <- numeric(length(q))
    zz <- (q - loc)/scale
    z[q > loc] <- exp(-zz[q > loc]^(-shape))
    Z[!nat] <- z
    if (!lower.tail) {
        Z <- 1 - Z
    }
    if (log.p)
        log(Z) else Z
}


#' @rdname Frechet
#' @export
qfrechet <- function(p, loc = 0, scale = 1, shape = 1, lower.tail = TRUE, log.p = FALSE) {
    if (scale <= 0 || shape <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, length(p)))
    }
    Z <- numeric(length(p))
    nat <- is.na(p)
    Z[nat] <- p[nat]
    p <- p[!nat]
    z <- numeric(length(p))
    if (log.p)
        p <- exp(p)
    if (!lower.tail)
        p <- 1 - p
    ok <- p >= 0 & p <= 1
    p <- p[ok]
    z[ok] <- loc + scale * (-log(p))^(-1/shape)
    z[!ok] <- NaN
    Z[!nat] <- z
    Z
}

#' @rdname Frechet
#' @export
rfrechet <- function(n, loc = 0, scale = 1, shape = 1) {
    if (scale <= 0 || shape <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, n))
    }
    loc + scale * (-log(runif(n)))^(-1/shape)
}



#' @title The Pareto Distribution
#' @description Density, distribution function, quantile function and random generation for the Pareto distribution with
#' scale and shape parameters.
#' @param x,q  vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param scale scale parameter.
#' @param shape shape parameter.
#' @param log,log.p logical; if TRUE, probabilities \eqn{p} are given as \eqn{log(p)}, default: FALSE.
#' @param lower.tail logical; if TRUE, probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}, default: TRUE.
#' @details The Pareto distribution function with scale parameter \eqn{s} and shape parameter
#'          \eqn{\alpha} has density given by
#'          \deqn{f(x)=\alpha s^\alpha/x^(\alpha+1)} for \eqn{x\ge s}. The cumulative distribution function is
#'          \deqn{F(x)=1-(s/x)^\alpha} for \eqn{x\ge s}. See \url{https://en.wikipedia.org/wiki/Pareto_distribution}
#'          for more details.
#' @return \code{dpareto} gives the density, \code{ppareto} gives the distribution function, \code{qpareto} gives the quantile function, and
#'         \code{rpareto} generates random deviates.
#'
#'         Invalid arguments will result in return value NaN, with a warning.
#' @examples
#' dpareto(seq(1, 5), 1, 1)
#' qpareto(ppareto(seq(1, 5), 1, 1), 1 ,1)
#' rpareto(5, 1, 1)
#' @rdname Pareto
#' @name Pareto
#' @seealso \code{\link{paretodist}}
#' @export
dpareto <- function(x, scale = 1, shape = 1, log = FALSE) {
    if (scale <= 0 || shape <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, length(x)))
    }
    Z <- numeric(length(x))
    nat <- is.na(x)
    Z[nat] <- x[nat]
    x <- x[!nat]
    z <- numeric(length(x))
    t <- x > scale
    z[t] <- shape/x[t] * (scale/x[t])^(shape)
    Z[!nat] <- z
    if (log)
        log(Z) else Z
}

#' @rdname Pareto
#' @export
ppareto <- function(q, scale = 1, shape = 1, lower.tail = TRUE, log.p = FALSE) {
    if (scale <= 0 || shape <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, length(q)))
    }
    Z <- numeric(length(q))
    nat <- is.na(q)
    Z[nat] <- q[nat]
    q <- q[!nat]
    z <- numeric(length(q))
    z[q > scale] <- 1 - (scale/q[q > scale])^(shape)
    Z[!nat] <- z
    if (!lower.tail) {
        Z <- 1 - Z
    }
    if (log.p)
        log(Z) else Z
}

#' @rdname Pareto
#' @export
qpareto <- function(p, scale = 1, shape = 1, lower.tail = TRUE, log.p = FALSE) {
    if (scale <= 0 || shape <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, length(p)))
    }
    Z <- numeric(length(p))
    nat <- is.na(p)
    Z[nat] <- p[nat]
    p <- p[!nat]
    z <- numeric(length(p))
    if (log.p)
        p <- exp(p)
    if (!lower.tail)
        p <- 1 - p
    ok <- p >= 0 & p <= 1
    p <- p[ok]
    z[ok] <- scale * (1 - p)^(-1/shape)
    z[!ok] <- NaN
    Z[!nat] <- z
    Z
}

#' @rdname Pareto
#' @export
rpareto <- function(n, scale = 1, shape = 1) {
    if (scale <= 0 || shape <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, n))
    }
    scale * (runif(n))^(-1/shape)
}



#' @title The Generalized Pareto Distribution
#' @description Density, distribution function, quantile function and random generation for the generalized Pareto distribution with
#' location, scale and shape parameters.
#' @param x,q  vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param loc location parameter.
#' @param scale scale parameter.
#' @param shape shape parameter.
#' @param log,log.p logical; if TRUE, probabilities \eqn{p} are given as \eqn{log(p)}, default: FALSE.
#' @param lower.tail logical; if TRUE, probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}, default: TRUE.
#' @details The generalized Pareto distribution function with location parameter \eqn{\mu}, scale parameter \eqn{\sigma}
#'          and shape parameter \eqn{\xi} has density given by
#'          \deqn{f(x)=1/\sigma  (1 + \xi z)^-(1/\xi + 1)} for \eqn{x\ge \mu} and \eqn{\xi> 0},
#'          or \eqn{\mu-\sigma/\xi \ge x\ge \mu} and \eqn{\xi< 0},
#'          where \eqn{z=(x-\mu)/\sigma}. In the case where \eqn{\xi= 0}, the density is equal to
#'          \eqn{f(x)=1/\sigma  e^-z} for \eqn{x\ge \mu}.
#'          The cumulative distribution function is
#'          \deqn{F(x)=1-(1+\xi z)^(-1/\xi)} for \eqn{x\ge \mu} and \eqn{\xi> 0},
#'          or \eqn{\mu-\sigma/\xi \ge x\ge \mu} and \eqn{\xi< 0},
#'          with \eqn{z} as stated above. If \eqn{\xi= 0} the CDF has form \eqn{F(x)=1-e^-z}.
#'
#'          See \url{https://en.wikipedia.org/wiki/Generalized_Pareto_distribution} for more details.
#' @return \code{dGPD} gives the density, \code{pGPD} gives the distribution function, \code{qGPD} gives the quantile function, and
#'         \code{rGPD} generates random deviates.
#'
#'         Invalid arguments will result in return value NaN, with a warning.
#' @examples
#' dGPD(seq(1, 5), 0, 1, 1)
#' qGPD(pGPD(seq(1, 5), 0, 1, 1), 0, 1 ,1)
#' rGPD(5, 0, 1, 1)
#' @rdname GPD
#' @name GPD
#' @seealso \code{\link{GPDdist}}
#' @export
dGPD <- function(x, loc = 0, scale = 1, shape = 0, log = FALSE) {
    if (scale <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, length(x)))
    }
    Z <- numeric(length(x))
    nat <- is.na(x)
    Z[nat] <- x[nat]
    x <- x[!nat]
    z <- numeric(length(x))
    t <- x >= loc
    zz <- (x - loc)/scale
    if (shape > 0) {
        z[t] <- (1/scale) * (1 + shape * zz[t])^(-(1/shape + 1))
    } else if (shape < 0) {
        tt <- x <= (loc - scale/shape)
        z[t & tt] <- (1/scale) * (1 + shape * zz[t & tt])^(-(1/shape + 1))
    } else {
        z[t] <- (1/scale) * exp(-zz[t])
    }
    Z[!nat] <- z
    if (log)
        log(Z) else Z
}

#' @rdname GPD
#' @export
pGPD <- function(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE) {
    if (scale <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, length(q)))
    }
    Z <- numeric(length(q))
    nat <- is.na(q)
    Z[nat] <- q[nat]
    q <- q[!nat]
    z <- numeric(length(q))
    t <- q >= loc
    zz <- (q - loc)/scale
    if (shape > 0) {
        z[t] <- 1 - (1 + shape * (zz[t]))^(-1/shape)
    } else if (shape < 0) {
        tt <- q <= (loc - scale/shape)
        ttt <- t & tt
        z[!tt] <- 1
        z[ttt] <- 1 - (1 + shape * (zz[ttt]))^(-1/shape)
    } else {
        z[t] <- 1 - exp(-zz[t])
    }
    Z[!nat] <- z
    if (!lower.tail) {
        Z <- 1 - Z
    }
    if (log.p)
        log(Z) else Z
}

#' @rdname GPD
#' @export
qGPD <- function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE) {
    if (scale <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, length(p)))
    }
    Z <- numeric(length(p))
    nat <- is.na(p)
    Z[nat] <- p[nat]
    p <- p[!nat]
    z <- numeric(length(p))
    if (log.p)
        p <- exp(p)
    if (!lower.tail)
        p <- 1 - p
    ok <- p >= 0 & p <= 1
    p <- p[ok]
    if (shape == 0)
        zz <- loc - scale * log(1 - p) else zz <- loc + ((1 - p)^(-shape) - 1) * scale/shape
    if (shape < 0)
        zz[p == 1] <- loc - scale/shape
    z[ok] <- zz
    z[!ok] <- NaN
    Z[!nat] <- z
    Z
}

#' @rdname GPD
#' @export
rGPD <- function(n, loc = 0, scale = 1, shape = 0) {
    if (scale <= 0) {
        warning("NaNs produced")
        return(rep.int(NaN, n))
    }
    if (shape == 0)
        loc - scale * log(runif(n)) else loc + ((runif(n))^(-shape) - 1) * scale/shape
}
