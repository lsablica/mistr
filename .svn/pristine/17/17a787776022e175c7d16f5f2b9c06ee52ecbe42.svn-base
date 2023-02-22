#' @title Creates an Object Representing Uniform Distribution
#' @description The function creates an object which represents the uniform distribution.
#' @param min minimum parameter, default: 0.
#' @param max maximum parameter, default: 1.
#' @return Object of class unifdist.
#' @details See \code{\link[stats]{Uniform}}.
#' @examples
#' U <- unifdist(1, 5)
#' d(U, c(2, 3, 4, NA))
#' r(U, 5)
#' @export
#' @rdname unifdist
#' @seealso \code{\link[stats]{Uniform}}
unifdist <- function(min = 0, max = 1) {
  if (!is.numeric(min) || !is.numeric(max))
    stop("Parameters must be a numeric")
  if (min >= max)
    stop("min must be smaller than max.")
  x <- list(parameters = list(min = min, max = max), type = "Uniform", support = list(from = min, to = max))
  class(x) <- c("unifdist", "contdist", "standist", "univdist", "dist")
  x
}


#' @title Creates an Object Representing Normal Distribution
#' @description The function creates an object which represents the normal distribution.
#' @param mean mean parameter, default: 0.
#' @param sd standard deviation parameter, default: 1.
#' @return Object of class normdist.
#' @details See \code{\link[stats]{Normal}}.
#' @examples
#' N <- normdist(1, 5)
#' d(N, c(2, 3, 4, NA))
#' r(N, 5)
#' @export
#' @rdname normdist
#' @seealso \code{\link[stats]{Normal}}
normdist <- function(mean = 0, sd = 1) {
  if (!is.numeric(mean) || !is.numeric(sd))
    stop("Parameters must be a numeric.")
  if (sd <= 0)
    stop("sd must be greater than 0.")
  x <- list(parameters = list(mean = mean, sd = sd), type = "Normal", support = list(from = -Inf, to = Inf))
  class(x) <- c("normdist", "contdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Student-t Distribution
#' @description The function creates an object which represents the Student-t distribution.
#' @param df degrees of freedom parameter, default: 2.
#' @return Object of class tdist.
#' @details See \code{\link[stats]{TDist}}.
#' @examples
#' t <- tdist(2)
#' d(t, c(2, 3, 4, NA))
#' r(t, 5)
#' @rdname tdist
#' @export
#' @seealso \code{\link[stats]{TDist}}
tdist <- function(df = 2) {
  if (!is.numeric(df))
    stop("Parameters must be a numeric")
  if (df <= 0)
    stop("df must be positive.")
  x <- list(parameters = list(df = df), type = "Student-t", support = list(from = -Inf, to = Inf))
  class(x) <- c("tdist", "contdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing F Distribution
#' @description The function creates an object which represents the F distribution.
#' @param df1 degrees of freedom parameter, default: 2.
#' @param df2 degrees of freedom parameter, default: 2.
#' @return Object of class fdist.
#' @details See \code{\link[stats]{FDist}}.
#' @examples
#' f <- fdist(2, 2)
#' d(f, c(2, 3, 4, NA))
#' r(f, 5)
#' @rdname fdist
#' @export
#' @seealso \code{\link[stats]{FDist}}
fdist <- function(df1 = 2, df2 = 2) {
  if (!is.numeric(df1) || !is.numeric(df2))
    stop("Parameters must be a numeric")
  if (df1 <= 0 || df2 <= 0)
    stop("df must be positive.")
  x <- list(parameters = list(df1 = df1, df2 = df2), type = "F", support = list(from = 0, to = Inf))
  class(x) <- c("fdist", "contdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Beta Distribution
#' @description The function creates an object which represents the beta distribution.
#' @param shape1  shape parameter, default: 2.
#' @param shape2  shape parameter, default: 2.
#' @return Object of class betadist.
#' @details See \code{\link[stats]{Beta}}.
#' @examples
#' B <- betadist(2, 2)
#' d(B, c(2, 3, 4, NA))
#' r(B, 5)
#' @rdname betadist
#' @export
#' @seealso \code{\link[stats]{Beta}}
betadist <- function(shape1 = 2, shape2 = 2) {
  if (!is.numeric(shape1) || !is.numeric(shape2))
    stop("Parameters must be a numeric")
  if (shape1 <= 0 || shape2 <= 0)
    stop("shape must be positive.")
  x <- list(parameters = list(shape1 = shape1, shape2 = shape2), type = "Beta", support = list(from = 0, to = 1))
  class(x) <- c("betadist", "contdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Cauchy Distribution.
#' @description The function creates an object which represents the Cauchy distribution.
#' @param location location parameter, default: 0.
#' @param scale scale parameter, default: 1.
#' @return Object of class cauchydist.
#' @details See \code{\link[stats]{Cauchy}}.
#' @examples
#' C <- cauchydist(0, 1)
#' d(C, c(2, 3, 4, NA))
#' r(C, 5)
#' @rdname cauchydist
#' @export
#' @seealso \code{\link[stats]{Cauchy}}
cauchydist <- function(location = 0, scale = 1) {
  if (!is.numeric(location) || !is.numeric(scale))
    stop("Parameters must be a numeric")
  if (scale <= 0)
    stop("scale must be positive.")
  x <- list(parameters = list(location = location, scale = scale), type = "Cauchy", support = list(from = -Inf, to = Inf))
  class(x) <- c("cauchydist", "contdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Chi-Squared Distribution
#' @description The function creates an object which represents the chi-squared distribution.
#' @param df degrees of freedom parameter, default: 2.
#' @return Object of class chisqdist.
#' @details See \code{\link[stats]{Chisquare}}.
#' @examples
#' C <- chisqdist(2)
#' d(C, c(2, 3, 4, NA))
#' r(C, 5)
#' @rdname chisqdist
#' @export
#' @seealso \code{\link[stats]{Chisquare}}
chisqdist <- function(df = 2) {
  if (!is.numeric(df))
    stop("Parameters must be a numeric")
  if (df <= 0)
    stop("df must be positive.")
  x <- list(parameters = list(df = df), type = "Chi-squared", support = list(from = 0, to = Inf))
  class(x) <- c("chisqdist", "contdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Exponential Distribution
#' @description The function creates an object which represents the exponential distribution.
#' @param rate rate parameter, default: 1.
#' @return Object of class expdist.
#' @details See \code{\link[stats]{Exponential}}.
#' @examples
#' E <- expdist(1)
#' d(E, c(2, 3, 4, NA))
#' r(E, 5)
#' @rdname expdist
#' @export
#' @seealso \code{\link[stats]{Exponential}}
expdist <- function(rate = 1) {
  if (!is.numeric(rate))
    stop("Parameters must be a numeric")
  if (rate <= 0)
    stop("rate must be positive.")
  x <- list(parameters = list(rate = rate), type = "Exponential", support = list(from = 0, to = Inf))
  class(x) <- c("expdist", "contdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Gamma Distribution
#' @description The function creates an object which represents the gamma distribution.
#' @param shape shape parameter, default: 2.
#' @param rate rate parameter, an alternative way to specify the scale.
#' @param scale scale parameter.
#' @return Object of class gammadist.
#' @details See \code{\link[stats]{GammaDist}}.
#' @examples
#' G <- gammadist(shape = 2, scale = 3)
#' d(G, c(2, 3, 4, NA))
#' r(G, 5)
#' @rdname gammadist
#' @export
#' @seealso \code{\link[stats]{GammaDist}}
gammadist <- function(shape = 2, rate, scale) {
  if (missing(scale) & !missing(rate)) {
    if (!is.numeric(shape) || !is.numeric(rate))
      stop("Parameters must be a numeric")
    if (shape <= 0 || rate <= 0)
      stop("shape and rate must be positive.")
    x <- list(parameters = list(shape = shape, rate = rate), type = "Gamma", support = list(from = 0, to = Inf))
  } else {
    if ((!missing(scale) & missing(rate)) || rate == 1/scale) {
      if (!is.numeric(shape) || !is.numeric(scale))
        stop("Parameters must be a numeric")
      if (shape <= 0 || scale <= 0)
        stop("shape and scale must be positive.")
      x <- list(parameters = list(shape = shape, scale = scale), type = "Gamma", support = list(from = 0, to = Inf))
    } else {
      stop("Scale must be the reciprocal of the rate.")
    }
  }
  class(x) <- c("gammadist", "contdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Log Normal Distribution.
#' @description The function creates an object which represents the log normal distribution.
#' @param meanlog mean parameter, default: 0.
#' @param sdlog standard deviation parameter, default: 1.
#' @return Object of class lnormdist.
#' @details See \code{\link[stats]{Lognormal}}.
#' @examples
#' L <- lnormdist(0, 1)
#' d(L, c(2, 3, 4, NA))
#' r(L, 5)
#' @rdname lnormdist
#' @export
#' @seealso \code{\link[stats]{Lognormal}}
lnormdist <- function(meanlog = 0, sdlog = 1) {
  if (!is.numeric(meanlog) || !is.numeric(sdlog))
    stop("Parameters must be a numeric.")
  if (sdlog <= 0)
    stop("sdlog must be greater than 0.")
  x <- list(parameters = list(meanlog = meanlog, sdlog = sdlog), type = "Log-Normal", support = list(from = 0, to = Inf))
  class(x) <- c("lnormdist", "contdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Weibull Distribution
#' @description The function creates an object which represents the Weibull distribution.
#' @param shape shape parameter, default: 1.
#' @param scale scale parameter, default: 1.
#' @return Object of class weibulldist.
#' @details See \code{\link[stats]{Weibull}}.
#' @examples
#' W <- weibulldist(1, 1)
#' d(W, c(2, 3, 4, NA))
#' r(W, 5)
#' @rdname weibulldist
#' @export
#' @seealso  \code{\link[stats]{Weibull}}
weibulldist <- function(shape = 1, scale = 1) {
  if (!is.numeric(shape) || !is.numeric(scale))
    stop("Parameters must be a numeric")
  if (shape <= 0 || scale <= 0)
    stop("shape and scale must be positive.")
  x <- list(parameters = list(shape = shape, scale = scale), type = "Weibull", support = list(from = 0, to = Inf))
  class(x) <- c("weibulldist", "contdist", "standist", "univdist", "dist")
  x
}


#' @title Creates an Object Representing Binomial Distribution.
#' @description The function creates an object which represents the binomial distribution.
#' @param size size parameter, default: 10.
#' @param prob probability parameter, default: 0.5.
#' @return Object of class binomdist.
#' @details See \code{\link[stats]{Binomial}}.
#' @examples
#' B <- binomdist(10, 0.4)
#' d(B, c(2, 3, 4, NA))
#' r(B, 5)
#' @rdname binomdist
#' @export
#' @seealso \code{\link[stats]{Binomial}}
binomdist <- function(size = 10, prob = 0.5) {
  if (!is.numeric(prob) || !is.numeric(size))
    stop("Parameters must be a numeric")
  if (prob < 0 || prob > 1)
    stop("prob must be in [0,1].")
  if (size != floor(size))
    stop("size must be an integer.")
  x <- list(parameters = list(size = size, prob = prob), type = "Binomial", support = list(from = 0, to = size, by = 1))
  class(x) <- c("binomdist", "discrdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Multinomial Distribution
#' @description The function creates an object which represents the multinomial distribution.
#' @param size size parameter, default: 10.
#' @param prob probability parameter vector, default: c(0.5, 0.5).
#' @return  Object of class multinomdist.
#' @details See \code{\link[stats]{Multinomial}}.
#' @examples
#' M <- multinomdist(10, c(0.5, 0.5))
#' d(M, c(7, 3))
#' r(M, 5)
#' @rdname multinomdist
#' @export
#' @seealso \code{\link[stats]{Multinomial}}
multinomdist <- function(size = 10, prob = c(0.5, 0.5)) {
  if (any(!is.numeric(prob)) || !is.numeric(size))
    stop("Parameters must be a numeric")
  if (any(prob < 0) || any(prob > 1))
    stop("prob must be in [0,1].")
  if (size != floor(size))
    stop("size must be an integer.")
  x <- list(parameters = list(size = size, prob = prob), type = "Multinomial", support = list(from = 0, to = size, by = 1))
  class(x) <- c("multinomdist", "discrdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Negative Binomial Distribution
#' @description The function creates an object which represents the negative binomial distribution.
#' @param size size parameter, default: 10.
#' @param prob probability parameter.
#' @param mu alternative parametrization via mean, see \code{\link[stats]{NegBinomial}}.
#' @return Object of class nbinomdist.
#' @details See \code{\link[stats]{NegBinomial}}.
#' @examples
#' N <- nbinomdist(10, 0.5)
#' d(N, c(2, 3, 4, NA))
#' r(N, 5)
#' @rdname nbinomdist
#' @export
#' @seealso \code{\link[stats]{NegBinomial}}
nbinomdist <- function(size = 10, prob, mu) {
  if (missing(mu) & !missing(prob)) {
    if (!is.numeric(prob) || !is.numeric(size))
      stop("Parameters must be a numeric")
    if (prob < 0 || prob > 1)
      stop("prob must be in [0,1].")
    if (size != floor(size))
      stop("size must be an integer.")
    x <- list(parameters = list(size = size, prob = prob), type = "Negative Binomial", support = list(from = 0, to = Inf,
                                                                                                      by = 1))
  } else {
    if (!missing(mu) & missing(prob)) {
      if (!is.numeric(mu) || !is.numeric(size))
        stop("Parameters must be a numeric")
      if (size != floor(size))
        stop("size must be an integer.")
      x <- list(parameters = list(size = size, mu = mu), type = "Negative Binomial", support = list(from = 0, to = Inf,
                                                                                                    by = 1))
    } else {
      stop("Either prob or mu has to be set.")
    }
  }
  class(x) <- c("nbinomdist", "discrdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Geometric Distribution
#' @description The function creates an object which represents the geometric distribution.
#' @param prob probability parameter, default: 0.5.
#' @return Object of class geomdist.
#' @details See \code{\link[stats]{Geometric}}.
#' @examples
#' G <- geomdist(0.5)
#' d(G, c(2, 3, 4, NA))
#' r(G, 5)
#' @rdname geomdist
#' @export
#' @seealso \code{\link[stats]{Geometric}}
geomdist <- function(prob = 0.5) {
  if (!is.numeric(prob))
    stop("Parameters must be a numeric")
  if (prob < 0 || prob > 1)
    stop("prob must be in [0,1].")
  x <- list(parameters = list(prob = prob), type = "Geometric", support = list(from = 0, to = Inf, by = 1))
  class(x) <- c("geomdist", "discrdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Hypergeometric Distribution
#' @description The function creates an object which represents the hypergeometric distribution.
#' @param m the number of white balls in the urn, default: 10.
#' @param n the number of black balls in the urn, default: 10.
#' @param k the number of balls drawn from the urn, default: 5.
#' @return Object of class hyperdist.
#' @details See \code{\link[stats]{Hypergeometric}}.
#' @examples
#' H <- hyperdist(0.5)
#' d(H, c(2, 3, 4, NA))
#' r(H, 5)
#' @rdname hyperdist
#' @export
#' @seealso \code{\link[stats]{Hypergeometric}}
hyperdist <- function(m = 10, n = 10, k = 5) {
  if (!is.numeric(m) || !is.numeric(n) || !is.numeric(k))
    stop("Parameters must be a numeric")
  if (k > m + n)
    stop("k must be smaller or equal to m+n.")
  x <- list(parameters = list(m = m, n = n, k = k), type = "Hypergeometric", support = list(from = max(0, k - n), to = min(k,
                                                                                                                           m), by = 1))
  class(x) <- c("hyperdist", "discrdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Poisson Distribution
#' @description The function creates an object which represents the Poisson distribution.
#' @param lambda mean parameter, default: 1.
#' @return Object of class poisdist.
#' @details See \code{\link[stats]{Poisson}}.
#' @examples
#' P <- poisdist(1)
#' d(P, c(2, 3, 4, NA))
#' r(P, 5)
#' @rdname poisdist
#' @export
#' @seealso \code{\link[stats]{Poisson}}
poisdist <- function(lambda = 1) {
  if (!is.numeric(lambda))
    stop("Parameters must be a numeric")
  if (lambda <= 0)
    stop("Lambda must be positive.")
  x <- list(parameters = list(lambda = lambda), type = "Poisson", support = list(from = 0, to = Inf, by = 1))
  class(x) <- c("poisdist", "discrdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Wilcoxon Distribution
#' @description The function creates an object which represents the Wilcoxon distribution.
#' @param m number of observations in the first sample.
#' @param n number of observations in the second sample.
#' @return Object of class wilcoxdist.
#' @details See \code{\link[stats]{Wilcoxon}}.
#' @examples
#' W <- wilcoxdist(20, 15)
#' d(W, c(2, 3, 4, NA))
#' r(W, 5)
#' @rdname wilcoxdist
#' @export
#' @seealso \code{\link[stats]{Wilcoxon}}
wilcoxdist <- function(m, n) {
  if (!is.numeric(m) || !is.numeric(n))
    stop("Parameters must be a numeric")
  if (m <= 0 || n <= 0)
    stop("m and n must be positive.")
  x <- list(parameters = list(m = m, n = n), type = "Wilcoxon Rank Sum", support = list(from = 0, to = m * n, by = 1))
  class(x) <- c("wilcoxdist", "discrdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Burr Distribution
#' @description The function creates an object which represents the Burr distribution.
#' @param shape1 shape parameter, default: 2.
#' @param shape2 shape parameter, default: 2.
#' @return Object of class burrdist.
#' @details See \code{\link{Burr}}.
#' @examples
#' B <- burrdist(2, 2)
#' d(B, c(2, 3, 4, NA))
#' r(B, 5)
#' @rdname burrdist
#' @export
#' @seealso \code{\link{Burr}}
burrdist <- function(shape1 = 2, shape2 = 2) {
  if (!is.numeric(shape1) || !is.numeric(shape2))
    stop("Parameters must be a numeric")
  if (shape1 <= 0 || shape2 <= 0)
    stop("shape must be positive.")
  x <- list(parameters = list(shape1 = shape1, shape2 = shape2), type = "Burr", support = list(from = 0, to = Inf))
  class(x) <- c("burrdist", "contdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Gumbel Distribution
#' @description The function creates an object which represents the Burr distribution.
#' @param loc location parameter, default: 0.
#' @param scale scale parameter, default: 1.
#' @return Object of class gumbeldist.
#' @details See \code{\link{Gumbel}}.
#' @examples
#' G <- gumbeldist(1, 2)
#' d(G, c(2, 3, 4, NA))
#' r(G, 5)
#' @rdname gumbeldist
#' @export
#' @seealso \code{\link{Gumbel}}
gumbeldist <- function(loc = 0, scale = 1) {
  if (!is.numeric(loc) || !is.numeric(scale))
    stop("Parameters must be a numeric")
  if (scale <= 0)
    stop("scale must be positive.")
  x <- list(parameters = list(loc = loc, scale = scale), type = "Gumbel", support = list(from = -Inf, to = Inf))
  class(x) <- c("gumbeldist", "contdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Frechet Distribution
#' @description The function creates an object which represents the Frechet distribution.
#' @param loc location parameter, default: 0.
#' @param scale scale parameter, default: 1.
#' @param shape shape parameter, default: 1.
#' @return Object of class frechetdist.
#' @details See \code{\link{Frechet}}.
#' @examples
#' Fr <- frechetdist(0, 1, 2)
#' d(Fr, c(2, 3, 4, NA))
#' r(Fr, 5)
#' @rdname frechetdist
#' @export
#' @seealso \code{\link{Frechet}}
frechetdist <- function(loc = 0, scale = 1, shape = 1) {
  if (!is.numeric(loc) || !is.numeric(scale) || !is.numeric(shape))
    stop("Parameters must be a numeric")
  if (scale <= 0 || shape <= 0)
    stop("scale and shape must be positive.")
  x <- list(parameters = list(loc = loc, scale = scale, shape = shape), type = "Frechet", support = list(from = loc, to = Inf))
  class(x) <- c("frechetdist", "contdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Pareto Distribution
#' @description The function creates an object which represents the Pareto distribution.
#' @param scale scale parameter, default: 1.
#' @param shape shape parameter, default: 1.
#' @return Object of class paretodist.
#' @details See \code{\link{Pareto}}.
#' @examples
#' P <- paretodist(1, 1)
#' d(P, c(2, 3, 4, NA))
#' r(P, 5)
#' @rdname paretodist
#' @export
#' @seealso \code{\link{Pareto}}
paretodist <- function(scale = 1, shape = 1) {
  if (!is.numeric(scale) || !is.numeric(shape))
    stop("Parameters must be a numeric")
  if (scale <= 0 || shape <= 0)
    stop("scale and shape must be positive.")
  x <- list(parameters = list(scale = scale, shape = shape), type = "Pareto", support = list(from = scale, to = Inf))
  class(x) <- c("paretodist", "contdist", "standist", "univdist", "dist")
  x
}

#' @title Creates an Object Representing Generalized Pareto Distribution
#' @description The function creates an object which represents the generalized Pareto distribution.
#' @param loc location parameter, default: 0.
#' @param scale scale parameter, default: 1.
#' @param shape shape parameter, default: 0.
#' @return Object of class GPDdist.
#' @details See \code{\link{GPD}}.
#' @examples
#' G <- GPDdist(0, 1, 0)
#' d(G, c(2, 3, 4, NA))
#' r(G, 5)
#' @rdname GPDdist
#' @export
#' @seealso \code{\link{GPD}}
GPDdist <- function(loc = 0, scale = 1, shape = 0) {
  if (!is.numeric(loc) || !is.numeric(scale) || !is.numeric(shape))
    stop("Parameters must be a numeric")
  if (scale <= 0)
    stop("scale must be positive.")
  x <- list(parameters = list(loc = loc, scale = scale, shape = shape), type = "Generalized Pareto", support = list(from = loc,
                                                                                                                    to = ifelse(shape < 0, loc - scale/shape, Inf)))
  class(x) <- c("GPDdist", "contdist", "standist", "univdist", "dist")
  x
}
