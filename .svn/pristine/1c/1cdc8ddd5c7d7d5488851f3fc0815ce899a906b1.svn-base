#' @export
print.standist <- function(x, digits = 4, ...) {
    Parameters <- paste(names(x$parameters), round(unlist(x$parameters), digits), sep = " = ", collapse = ", ")
    Distribution <- x$type
    d <- data.frame(Distribution, Parameters)
    name.width <- unlist(lapply(names(d), nchar), use.names = FALSE)
    print(format(d, trim = TRUE, width = name.width, justify = "centre"), row.names = FALSE)
    invisible(x)
}

#' @title Density Function
#' @description \code{d} is a generic function that evaluates the density function of a distribution object at given values.
#' @param O distribution object.
#' @param x vector of quantiles.
#' @param log logical; if TRUE, probabilities \eqn{p} are given as \eqn{log(p)}, default: FALSE.
#' @return Vector of computed results.
#' @details Methods of \code{d} function evaluates any offered
#'          distribution from the package \code{\link{mistr}}. The function makes use of the d[sufix] functions
#'          as \code{dnorm} or \code{dbeta} and thus, if a new distribution is added,
#'          these functions must be reachable through the search path.
#' @examples
#' N <- normdist(1, 3)
#' d(N, c(NA, 1, 3, 5))
#'
#' C <- cauchydist()
#' M <- mixdist(N, C, weights = c(0.5, 0.5))
#' d(M, c(NA, 1, 3, 5))
#'
#' CC <- compdist(N, C, weights = c(0.5, 0.5), breakpoints = 1)
#' CCC <- 2*C+5
#' d(CCC, c(NA, 1, 3, 5))
#' @rdname d
#' @export
d <- function(O, x, log = FALSE) UseMethod("d")

#' @rdname d
#' @export
d.standist <- function(O, x, log = FALSE) {
    func <- get(paste0("d", strsplit(class(O)[1], "dist")[[1]]), mode = "function")
    suppressWarnings(do.call(func, c(list(x = x), O$parameters, list(log = log))))
}


#' @title Distribution Function
#' @description \code{p} is a generic function that evaluates the distribution function of a distribution object at given values.
#' @param O distribution object.
#' @param q vector of quantiles.
#' @param lower.tail logical; if TRUE, probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}, default: TRUE.
#' @param log.p logical; if TRUE, probabilities \eqn{p} are given as \eqn{log(p)}, default: FALSE.
#' @return Vector of computed results.
#' @details Methods of \code{p} function evaluates any offered
#'          distribution from the package \code{\link{mistr}}. The function makes use of the p[sufix] functions
#'          as \code{pnorm} or \code{pbeta} and thus, if a new distribution is added,
#'          these functions must be reachable through the search path.
#' @examples
#' N <- normdist(1,3)
#' p(N, c(NA,1,3,5))
#'
#' C <- cauchydist()
#' M <- mixdist(N, C, weights = c(0.5, 0.5))
#' p(M, c(NA,1,3,5))
#'
#' CC <- compdist(N, C, weights = c(0.5, 0.5), breakpoints = 1)
#' CCC <- 2*C+5
#' p(CCC, c(NA,1,3,5))
#' @rdname p
#' @export
p <- function(O, q, lower.tail = TRUE, log.p = FALSE) UseMethod("p")

#' @rdname p
#' @export
p.standist <- function(O, q, lower.tail = TRUE, log.p = FALSE) {
    func <- get(paste0("p", strsplit(class(O)[1], "dist")[[1]]), mode = "function")
    do.call(func, c(list(q = q), O$parameters, list(lower.tail = lower.tail, log.p = log.p)))
}

#' @title Quantile Function
#' @description \code{q} is a generic function that evaluates the quantile function of a distribution object at given values.
#' @param O distribution object.
#' @param p vector of probabilities.
#' @param lower.tail logical; if TRUE, probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}, default: TRUE.
#' @param log.p logical; if TRUE, probabilities \eqn{p} are given as \eqn{log(p)}, default: FALSE.
#' @param ... further arguments to be passed.
#' @return Vector of computed results.
#' @details Methods of \code{q} function  evaluates any offered
#'          distribution from package \code{\link{mistr}}. The function makes use of the q[sufix] functions
#'          as \code{qnorm} or \code{qbeta} and thus, if a new distribution is added,
#'          these functions must be reachable through the search path.
#'
#'          The mixture method \code{\link{q.mixdist}} and the default
#'          method \code{\link{q.default}} have its own help page.
#' @examples
#' N <- normdist(1, 3)
#' q(N, c(NA, 1, 3, 5))
#'
#' C <- cauchydist()
#' CC <- compdist(N, C, weights = c(0.5, 0.5), breakpoints = 1)
#' CCC <- 2*C+5
#' q(CCC, c(NA, 1, 3, 5))
#' @rdname q
#' @export
q <- function(O, p, lower.tail = TRUE, log.p = FALSE, ...) UseMethod("q")

#' @rdname q
#' @export
q.standist <- function(O, p, lower.tail = TRUE, log.p = FALSE, ...) {
    func <- get(paste0("q", strsplit(class(O)[1], "dist")[[1]]), mode = "function")
    do.call(func, c(list(p = p), O$parameters, list(lower.tail = lower.tail, log.p = log.p)))
}

#' @title Terminate an R Session
#' @description The default method \code{q.default} terminates the current R session.
#' @param O place holder for generic, by default set to save, default: save.
#' @param p place holder for generic, by default set to status, default: status.
#' @param lower.tail place holder for generic, by default set to runLast, default: runLast.
#' @param log.p place holder for generic, default: FALSE.
#' @param save a character string indicating whether the environment (workspace) should be saved, one of "no", "yes", "ask" or "default", default: 'default'.
#' @param status the (numerical) error status to be returned to the operating system, where relevant. Conventionally 0 indicates successful completion, default: 0.
#' @param runLast should .Last() be executed?, default: TRUE.
#' @param ... further arguments to be passed.
#' @details This method is designed to quit R if the \code{q()} without a distribution is called.
#'          The reason for such an implementation is R-Studio in Linux and Mac systems, where
#'          the software calls \code{q()} (rather than \code{base::q()}) once the R-Studio window
#'          is closed. Such implementation solves the issued with the overwriting of \code{q()}.
#' @rdname q.default
#' @seealso \code{\link[base]{q}}
#' @export
q.default <- function(O = save, p = status, lower.tail = runLast, log.p = FALSE, save = "default", status = 0, runLast = TRUE, ...){
  base::q(save = O, status = p, runLast = lower.tail)
}

# q.default <- function(save = "default", status = 0, runLast = TRUE){
#   base::q(save = save, status = status, runLast = runLast)
# }

#' @title Mistr d/p/q/r Wrappers
#' @description The functions \code{mistr_d}, \code{mistr_p}, \code{mistr_q}, \code{mistr_r} are wrappers
#'               for \code{\link{d}}, \code{\link{p}}, \code{\link{q}} and \code{\link{r}}, respectively.
#' @param O distribution object.
#' @param x,q  vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param log,log.p logical; if TRUE, probabilities \eqn{p} are given as \eqn{log(p)}, default: FALSE.
#' @param lower.tail logical; if TRUE, probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}, default: TRUE.
#' @param ... further arguments to be passed.
#' @return Vector of computed results.
#' @details Wrappers are offered as a consequence of R-Studio in Windows OS
#'          where the \code{q()} calls in the console are caught
#'          and terminate the \code{R} session.
#' @rdname mistr_d_p_q_r
#' @name mistr_d_p_q_r
NULL

#' @rdname mistr_d_p_q_r
#' @export
mistr_d <- function(O, x, log = FALSE) {
  d(O = O, x = x, log = log)
}

#' @rdname mistr_d_p_q_r
#' @export
mistr_p <- function(O, q, lower.tail = TRUE, log.p = FALSE) {
  p(O = O, q = q, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname mistr_d_p_q_r
#' @export
mistr_q <- function(O, p, lower.tail = TRUE, log.p = FALSE, ...) {
  q(O = O, p = p, lower.tail = lower.tail, log.p = log.p, ...)
}

#' @rdname mistr_d_p_q_r
#' @export
mistr_r <- function(O, n) {
  r(O = O, n = n)
}

#' @title Random Generation
#' @description \code{r} is a generic function that generates random deviates of a distribution object.
#' @param O distribution object.
#' @param n number of observations.
#' @return Vector of computed results.
#' @details Methods of \code{r} function generates random deviates of offered
#'          distribution from the package \code{\link{mistr}}. The function makes use of the r[sufix] functions
#'          as \code{rnorm} or \code{rbeta} and thus, if a new distribution is added,
#'          these functions must be reachable through the search path.
#'
#'          For more complicated composite distributions, where one of the components is a mixture distribution,
#'          the function performs a rejection sampling of mixture random numbers to improve the speed.
#' @examples
#' N <- normdist(1, 3)
#' r(N, 5)
#'
#' C <- cauchydist()
#' M <- mixdist(N, C, weights = c(0.5, 0.5))
#' r(M, 5)
#'
#' CC <- compdist(N, C, weights = c(0.5, 0.5), breakpoints = 1)
#' CCC <- 2*C+5
#' r(CCC, 5)
#' @rdname r
#' @export
r <- function(O, n) UseMethod("r")

#' @rdname r
#' @export
r.standist <- function(O, n) {
    func <- get(paste0("r", strsplit(class(O)[1], "dist")[[1]]), mode = "function")
    do.call(func, c(list(n = n), O$parameters))
}

#' @rdname r
#' @export
r.hyperdist <- function(O, n) {
    func <- get(paste0("r", strsplit(class(O)[1], "dist")[[1]]), mode = "function")
    do.call(func, c(list(nn = n), O$parameters))
}

#' @rdname r
#' @export
r.wilcoxdist <- function(O, n) {
    func <- get(paste0("r", strsplit(class(O)[1], "dist")[[1]]), mode = "function")
    do.call(func, c(list(nn = n), O$parameters))
}


#' @title Left-Hand Limit of Distribution Function
#' @description \code{plim} is a generic function that evaluates the left-hand limit of distribution function
#'              for a distribution object at given values.
#' @param O distribution object.
#' @param q vector of quantiles.
#' @param lower.tail logical; if TRUE, probabilities are \eqn{P[X < x]} otherwise, \eqn{P[X \ge x]}, default: TRUE.
#' @param log.p logical; if TRUE, probabilities \eqn{p} are given as \eqn{log(p)}, default: FALSE.
#' @return Vector of computed results.
#' @details Methods of \code{plim} function evaluates the left-hand limit of any offered
#'          distribution from the package \code{\link{mistr}}. The left-hand limit is defined as
#'          \eqn{F(x-)=P(X<x)}.
#'          The function makes use of the p[sufix] and q[sufix] functions
#'          as \code{pnorm} or \code{qbeta} and thus, if a new distribution is added,
#'          these functions must be reachable through the search path.
#' @examples
#' B <- binomdist(10, 0.3)
#' plim(B, c(NA, 1, 3, 5))
#'
#' P <- poisdist()
#' M <- mixdist(B, P, weights = c(0.5, 0.5))
#' plim(M, c(NA, 1, 3, 5))
#'
#' CC <- compdist(B, P, weights = c(0.5, 0.5), breakpoints = 1)
#' CCC <- 2*CC+5
#' plim(CCC, c(NA, 1, 3, 5))
#' @rdname plim
#' @export
plim <- function(O, q, lower.tail = TRUE, log.p = FALSE) UseMethod("plim")

#' @rdname plim
#' @export
plim.discrdist <- function(O, q, lower.tail = TRUE, log.p = FALSE) {
  Z <- numeric(length(q))
  nat <- is.na(q)
  Z[nat] <- q[nat]
  q <- q[!nat]
  z <- numeric(length(q))
  h <- near(q%%O$support$by, O$support$from%%O$support$by) | near(q%%O$support$by, O$support$from%%O$support$by + O$support$by)
  z[h] <- q[h] - O$support$by/2
  z[!h] <- q[!h]
  Z[!nat] <- p.standist(O, z, lower.tail = lower.tail, log.p = log.p)
  Z
}
# plim.discrdist <- function(O, q, lower.tail = TRUE, log.p = FALSE) {
#   Z <- numeric(length(q))
#   nat <- is.na(q)
#   Z[nat] <- q[nat]
#   q <- q[!nat]
#   z <- numeric(length(q))
#   dist <- p.standist(O, z, lower.tail = lower.tail)
#   dens <- d.standist(O, z)
#   res <- dist + if(lower.tail) -dens else +dens
#   res[res<=0] <- 0
#   Z[!nat] <- res
#   Z
# }
# OLD  -  using modulus


#' @rdname plim
#' @export
plim.contdist <- function(O, q, lower.tail = TRUE, log.p = FALSE) {
    p.standist(O, q, lower.tail = lower.tail, log.p = log.p)
}


#' @title Right-Hand Limit of Quantile Function
#' @description \code{qlim} is a generic function that evaluates the right-hand limit of quantile function
#'              for a distribution object at given values.
#' @param O distribution object.
#' @param p vector of probabilities.
#' @param lower.tail logical; if TRUE, probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}, default: TRUE.
#' @param log.p logical; if TRUE, probabilities \eqn{p} are given as \eqn{log(p)}, default: FALSE.
#' @return Vector of computed results.
#' @details Methods of \code{qlim} function evaluates the right-hand limit of any offered
#'          distribution object from the package \code{\link{mistr}}. The right-hand limit
#'          of a quantile function is defined as
#'          \deqn{Q(x+)=inf{x: p<P(X\le x)}.}
#'          The function makes use of the p[sufix] and q[sufix] functions
#'          as \code{pnorm}, \code{pbeta}, \code{qnorm}, \code{qbeta}, and thus, if a new distribution is added,
#'          these functions must be reachable through the search path.
#'
#'          Methods for \code{\link[=qlim.discrmixdist]{mixtures}} have its own help page.
#' @examples
#' B <- binomdist(10, 0.3)
#' qlim(B, plim(B, c(NA, 1, 3, 5)))
#'
#' P <- poisdist()
#' M <- mixdist(B, P, weights = c(0.5, 0.5))
#' qlim(M, plim(M, c(NA, 1, 3, 5)))
#'
#' CC <- compdist(B, P, weights = c(0.5, 0.5), breakpoints = 1)
#' CCC <- 2*CC+5
#' qlim(CCC, plim(CCC, c(NA, 1, 3, 5)))
#' @rdname qlim
#' @export
qlim <- function(O, p, lower.tail = TRUE, log.p = FALSE) UseMethod("qlim")

#' @rdname qlim
#' @export
qlim.discrdist <- function(O, p, lower.tail = TRUE, log.p = FALSE) {
    add <- get_opt("add") 
    Z <- numeric(length(p))
    nat <- is.na(p)
    Z[nat] <- p[nat]
    p <- p[!nat]
    z <- numeric(length(p))
    ok <- p >= 0 & p <= 1
    p <- p[ok]
    zz <- numeric(length(p))
    zz[p == 1] <- unname(sudo_support(O)[2])
    zz[p == 0] <- unname(sudo_support(O)[1])
    supp <- p > 0 & p < 1
    p <- p[supp]
    zzz <- numeric(length(p))
    q <- q(O, p, lower.tail = lower.tail, log.p = log.p)
    p2 <- p(O, q, lower.tail = lower.tail, log.p = log.p)
    jump <- near(p, p2)
    p[jump] <- p[jump] + if (lower.tail) add else -add
    zzz <- q(O, p, lower.tail = lower.tail, log.p = log.p)
    zz[supp] <- zzz
    z[ok] <- zz
    z[!ok] <- NaN
    Z[!nat] <- z
    Z
}

#' @rdname qlim
#' @export
qlim.contdist <- function(O, p, lower.tail = TRUE, log.p = FALSE) {
    q.standist(O, p, lower.tail = lower.tail, log.p = log.p)
}
