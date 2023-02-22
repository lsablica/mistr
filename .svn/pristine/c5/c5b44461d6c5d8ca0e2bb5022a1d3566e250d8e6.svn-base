#' @title Creates an Object Representing Mixture Distribution
#' @description \code{mixdist} creates an object which represents the mixture distribution.
#' @param ... distribution objects.
#' @param weights vector of weights for the components.
#' @param dist vector of distribution names.
#' @param params list of parameters for each component.
#' @return Object of class mixdist.
#' @details A CDF of a mixture distribution function is \deqn{F(A)=\sum w_{i}F_{i}(A)}, where
#'          \eqn{w_{i}} is the weight of the i-th component and \eqn{F_{i}()} is the CDF of the i-th component.
#'
#'          The objects can be specified in two ways, either the user may enter
#'          distribution objects or a vector of names and list of parameters.
#'          See the examples below.
#' @examples
#' # using the objects
#' M <- mixdist(normdist(1, 3), expdist(4), weights = c(0.7, 0.3))
#' M
#'
#' # using the names and parameters
#' M2 <- mixdist(c("norm", "exp"), list(c(mean = 1, sd = 3), c(rate = 4)),
#'               weights = c(0.7, 0.3))
#' M2
#' @export
#' @rdname mixdist
#' @seealso \code{\link{compdist}}
mixdist <- function(..., weights) UseMethod("mixdist")
#' @rdname mixdist
#' @export
mixdist.dist <- function(..., weights) {
    K <- list(...)
    if (length(K) != length(weights))
        stop("you need the same length of probabilities and distributions")
    if (length(K) == 1)
        return(K[[1]])
    if (sum(weights) != 1)
        weights <- weights/sum(weights)
    x <- list(objects = K, weights = weights)
    if (all(unlist(lapply(K, is.contin)))) {
        class(x) <- c("contmixdist", "mixdist", "univdist", "dist")
    } else if (all(unlist(lapply(K, is.discrete)))) {
        class(x) <- c("discrmixdist", "mixdist", "univdist", "dist")
    } else {
        class(x) <- c("contdiscrmixdist", "mixdist", "univdist", "dist")
    }
    x
}
#' @rdname mixdist
#' @export
mixdist.default <- function(dist, params, weights,...) {
    if (length(dist) != length(weights))
        stop("you need the same length of probabilities and distributions")
    if (sum(weights) != 1)
        weights <- weights/sum(weights)
    K <- lapply(1:length(dist), function(i) {
        D <- get(paste0(dist[i], "dist"), mode = "function")
        do.call(D, as.list(params[[i]]))
    })
    if (length(K) == 1)
        return(K[[1]])
    x <- list(objects = K, weights = weights)
    if (all(unlist(lapply(K, function(O) {
        any(class(O) == "contdist")
    })))) {
        class(x) <- c("contmixdist", "mixdist", "univdist", "dist")
    } else if (all(unlist(lapply(K, function(O) {
        any(class(O) == "discrdist")
    })))) {
        class(x) <- c("discrmixdist", "mixdist", "univdist", "dist")
    } else {
        class(x) <- c("contdiscrmixdist", "mixdist", "univdist", "dist")
    }
    x
}

#' @export
print.mixdist <- function(x, digits = 4, ...) {
    Distribution <- numeric(length(x$objects))
    Parameters <- numeric(length(x$objects))
    mix <- unlist(lapply(x$objects, is.mixture), use.names = FALSE)
    comp <- unlist(lapply(x$objects, is.composite), use.names = FALSE)
    dist <- !(mix | comp)
    Distribution[mix] <- "Mixture distribution"
    Distribution[comp] <- "Composite distribution"
    Distribution[dist] <- unlist(lapply(x$objects[dist], function(x) x$type), use.names = FALSE)
    Trafo <- sapply(x$objects, function(x) {
        if (is.transformed(x)) {
            deparse(x$trafo$print)
        } else "none"
    })
    Parameters[dist] <- unlist(lapply(x$objects[dist], function(x) {
        paste(names(x$parameters), round(unlist(x$parameters), digits), sep = " = ", collapse = ", ")
    }), use.names = FALSE)
    Parameters[!dist] <- "   "
    Weight <- round(x$weights, digits)
    if (any(Trafo != "none")) {
        results <- data.frame(Trafo, Distribution, Parameters, Weight)
    } else {
        results <- data.frame(Distribution, Parameters, Weight)
    }
    name.width <- unlist(lapply(names(results), nchar), use.names = FALSE)
    cat("Mixture distribution with: \n \n")
    print(format(results, width = name.width, justify = "centre"))
    invisible(x)
}


#' @rdname d
#' @export
d.mixdist <- function(O, x, log = FALSE) {
    k <- sapply(O$objects, function(H) d(H, x))
    Z <- as.numeric(k %*% O$weights)
    if (log)
        log(Z) else Z
}

#' @rdname p
#' @export
p.mixdist <- function(O, q, lower.tail = TRUE, log.p = FALSE) {
    k <- sapply(O$objects, function(H) p(H, q))
    Z <- as.numeric(k %*% O$weights)
    if (!lower.tail)
        Z <- 1 - Z
    if (log.p)
        log(Z) else Z
}

#' @rdname plim
#' @export
plim.mixdist <- function(O, q, lower.tail = TRUE, log.p = FALSE) {
    k <- sapply(O$objects, function(H) plim(H, q))
    Z <- as.numeric(k %*% O$weights)
    if (!lower.tail)
        Z <- 1 - Z
    if (log.p)
        log(Z) else Z
}

#' @title Quantile Function of a Mixture Model
#' @description \code{q.mixdist} is a method that evaluates the quantile function of a mixture distribution object at given values.
#' @param O mixture distribution object.
#' @param p vector of probabilities.
#' @param lower.tail logical; if TRUE, probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}, default: TRUE.
#' @param log.p logical; if TRUE, probabilities \eqn{p} are given as \eqn{log(p)}, default: FALSE.
#' @param ... further arguments to be passed.
#' @return Vector of computed results.
#' @details Methods of \code{q} function evaluates any offered
#'          distribution from the package \code{\link{mistr}}. The function makes use of the p[sufix] and q[sufix] functions
#'          as \code{pnorm} or \code{qbeta} and thus, if a new distribution is added,
#'          these functions must be reachable through the search path.
#'
#'          The values are numerically found using the \code{\link[stats]{uniroot}} function, while the starting intervals are found
#'          automatically. The option parameter \code{tol} specifies the tolerance for the \code{\link[stats]{uniroot}}. 
#'          Options parameter \code{sub} is used to test whether the CDF at computed values minus \code{sub} is not the same and thus the given value is not an
#'          infimum. In such case, the root is found one more time for the value \code{p - sub}.
#'
#'          Other methods \code{\link{q}} and the default
#'          method \code{\link{q.default}} have its own help page.
#' @examples
#' DM <- mixdist(3*binomdist(12, 0.4), -2*poisdist(2)+12, weights=c(0.5, 0.5))
#' y <- c(0.4, p(DM, c(5, 10, 15, 18)), 0.95)
#' x <- q(DM, y)
#' plot(DM, which = "cdf", only_mix=TRUE, xlim1 = c(0, 37))
#' points(x, y)
#' @seealso \code{\link{set_opt}}
#' @rdname q_mixdist
#' @export
q.mixdist <- function(O, p, lower.tail = TRUE, log.p = FALSE,...) {
    tol <- get_opt("tol")
    sub <- get_opt("sub")
    
    Z <- numeric(length(p))
    nat <- is.na(p)
    Z[nat] <- p[nat]
    p <- p[!nat]
    if (log.p)
        p <- exp(p)
    if (!lower.tail)
        p <- 1 - p
    zz <- numeric(length(p))
    ok <- p >= 0 & p <= 1
    p <- p[ok]
    z <- numeric(length(p))
    z[p == 0] <- unname(sudo_support(O)[1])
    z[p == 1] <- unname(sudo_support(O)[2])
    inside <- p < 1 & p > 0
    if (sum(inside) > 0) {
        p <- p[inside]
        qf <- function(p) {
            interval <- range(unlist(lapply(O$objects, function(X) q(X, p)), use.names = FALSE))
            if (interval[1] == interval[2]) {return(interval[1])}
            uniroot(function(x) p(O, x) - p, interval = interval, tol = tol, extendInt = "upX")$root
        }
        qq <- unlist(lapply(p, qf))
        TF <- p(O, qq - sub) == p
        qq[TF] <- unlist(lapply(p[TF] - sub, qf))
        bad_unitroot <- p(O, qq) < p
        qq[bad_unitroot] <- qq[bad_unitroot] + tol
        z[inside] <- qq
    }
    zz[ok] <- z
    zz[!ok] <- NaN
    Z[!nat] <- zz
    Z
}
#' @title Right-Hand Limit of Mixture Quantile Function
#' @description  \code{qlim.mixdist} is a method that evaluates the right-hand limit of quantile function for a mixture distribution object at given values.
#' @param O mixture distribution object.
#' @param p vector of probabilities.
#' @param lower.tail logical; if TRUE, probabilities are \eqn{P[X \le x]} otherwise, \eqn{P[X > x]}, default: TRUE.
#' @param log.p logical; if TRUE, probabilities \eqn{p} are given as \eqn{log(p)}, default: FALSE.
#' @param ... further arguments to be passed.
#' @return Vector of computed results.
#' @details Methods of \code{qlim} function evaluates the right-hand limit of a quantile function for any offered
#'          distribution object from the package \code{\link{mistr}}. The right-hand limit
#'          of a quantile function is defined as
#'          \deqn{Q(x+)=inf{x: p<P(X\le x)}.}
#'          The function makes use of the p[sufix] and q[sufix] functions
#'          as \code{pnorm}, \code{pbeta}, \code{qnorm}, \code{qbeta}, and thus, if a new distribution will be added,
#'          these functions must be reachable through the search path.
#'
#'          The values are numerically found using the \code{\link[stats]{uniroot}} function, while the starting intervals are found
#'          automatically. The option parameter \code{tol} specifies the tolerance for the \code{\link[stats]{uniroot}}. 
#'          Options parameter \code{sub} is used to test whether the CDF at computed value minus \code{sub} is not the same and thus the given value is not an
#'          infimum. In such case, the root is found one more time for the value \code{p - sub}.
#'
#'          Other methods \code{\link{qlim}} have its own help page.
#' @examples
#' # q() of a negative transformed random variable uses qlim()
#' DM <- mixdist(3*binomdist(12,0.4), -2*poisdist(2)+12, weights=c(0.5, 0.5))
#' y <- c(0.05, 0.4, p(-DM, c(-5, -10, -15, -18)), 0.95)
#' x <- q(-DM, y)
#' plot(-DM, which = "cdf", only_mix=TRUE, xlim1 = c(-37, 0))
#' points(x, y)
#' @seealso \code{\link{set_opt}}
#' @rdname qlim_mixdist
#' @export
qlim.discrmixdist <- function(O, p, lower.tail = TRUE, log.p = FALSE) {
    add <- get_opt("add")
    tol <- get_opt("tol")
    sub <- get_opt("sub")
    
    Z <- numeric(length(p))
    nat <- is.na(p)
    Z[nat] <- p[nat]
    p <- p[!nat]
    if (log.p)
        p <- exp(p)
    if (!lower.tail)
        p <- 1 - p
    zz <- numeric(length(p))
    ok <- p >= 0 & p <= 1
    p <- p[ok]
    z <- numeric(length(p))
    z[p == 0] <- unname(sudo_support(O)[1])
    z[p == 1] <- unname(sudo_support(O)[2])
    inside <- p < 1 & p > 0
    if (sum(inside) > 0) {
        p <- p[inside]
        qf <- function(p) {
            interval <- range(unlist(lapply(O$objects, function(X) q(X, p)), use.names = FALSE))
            if (interval[1] == interval[2]){
              return(interval[1])
            }
            ju <- jumps(O, interval)
            if (!is.null(ju) && any(near(p, p(O, ju)))) p <- p + add
            uniroot(function(x) p(O, x) - p, interval = interval, tol = tol, extendInt = "upX")$root
        }
        qq <- unlist(lapply(p, qf))
        TF <- p(O, qq - sub) == p
        qq[TF] <- unlist(lapply(p[TF] - sub, qf))
        bad_unitroot <- p(O, qq) < p
        qq[bad_unitroot] <- qq[bad_unitroot] + tol
        z[inside] <- qq
    }
    zz[ok] <- z
    zz[!ok] <- NaN
    Z[!nat] <- zz
    Z
}
#' @rdname qlim_mixdist
#' @export
qlim.contdiscrmixdist <- function(O, p, lower.tail = TRUE, log.p = FALSE) {
    qlim.discrmixdist(O, p, lower.tail = lower.tail, log.p = log.p)
}
#' @rdname qlim_mixdist
#' @export
qlim.contmixdist <- function(O, p, lower.tail = TRUE, log.p = FALSE) {
  q.mixdist(O, p, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname r
#' @export
r.mixdist <- function(O, n) {
    Z <- numeric(n)
    l <- length(O$objects)
    t <- sample(l, n, replace = TRUE, prob = O$weights)
    v <- unlist(lapply(1:l, function(x) sum(t == x)), use.names = FALSE)
    lapply(1:l, function(i) {
        Z[t == i] <<- r(O$objects[[i]], v[i])
    })
    Z
}

#' @export
print.trans_mixdist <- function(x, digits = 4, ...) {
    Distribution <- numeric(length(x$objects))
    Parameters <- numeric(length(x$objects))
    mix <- unlist(lapply(x$objects, is.mixture), use.names = FALSE)
    comp <- unlist(lapply(x$objects, is.composite), use.names = FALSE)
    dist <- !(mix | comp)
    Distribution[mix] <- "Mixture distribution"
    Distribution[comp] <- "Composite distribution"
    Distribution[dist] <- unlist(lapply(x$objects[dist], function(x) x$type), use.names = FALSE)
    Trafo <- sapply(x$objects, function(x) {
        if (is.transformed(x)) {
            deparse(x$trafo$print)
        } else "none"
    })
    Parameters[dist] <- unlist(lapply(x$objects[dist], function(x) {
        paste(names(x$parameters), round(unlist(x$parameters), digits), sep = " = ", collapse = ", ")
    }), use.names = FALSE)
    Parameters[!dist] <- "   "
    Weight <- round(x$weights, digits)
    if (any(Trafo != "none")) {
        results <- data.frame(Trafo, Distribution, Parameters, Weight)
    } else {
        results <- data.frame(Distribution, Parameters, Weight)
    }
    name.width <- unlist(lapply(names(results), nchar), use.names = FALSE)
    cat("Monotonically transformed mixture distribution with: \n \n")
    cat("Trafo: ", deparse(x$trafo$print), "\n \n")
    print(format(results, width = name.width, justify = "centre"))
    invisible(x)
}


#' @rdname d
#' @export
d.trans_mixdist <- function(O, x, log = FALSE) {
    V <- lapply(O$objects, function(x) eval(O$trafo$print, list(X = x)))
    k <- sapply(V, function(H) d(H, x))
    Z <- as.numeric(k %*% O$weights)
    if (log)
        log(Z) else Z
}

