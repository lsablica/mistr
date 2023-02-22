#' @export
print.trans_standist <- function(x, digits = 4, ...) {
    Trafo <- deparse(x$trafo$print)
    Parameters <- paste(names(x$parameters), round(unlist(x$parameters), digits), sep = " = ", collapse = ", ")
    Distribution <- x$type
    d <- data.frame(Trafo, Distribution, Parameters)
    name.width <- unlist(lapply(names(d), nchar), use.names = FALSE)
    print(format(d, trim = TRUE, width = name.width, justify = "centre"), row.names = FALSE)
    invisible(x)
}

#' @rdname p
#' @export
p.trans_univdist <- function(O, q, lower.tail = TRUE, log.p = FALSE) {
    Z <- numeric(length(q))
    nat <- is.na(q)
    Z[nat] <- q[nat]
    q <- q[!nat]
    z <- numeric(length(q))
    r <- findInterval(q, sudo_support(O))
    q <- eval(O$trafo$invtrans, list(X = q[r == 1]))
    if (lower.tail)
        z[r == 2] <- 1 else z[r == 0] <- 1
    if (log.p)
        z[r != 1] <- log(z[r != 1])
    if (length(q) > 0) {
        if (monot(O) == 1) {
            z[r == 1] <- p(untrafo.trans_standist(O), q, lower.tail = lower.tail, log.p = log.p)
        } else {
            z[r == 1] <- plim(untrafo.trans_standist(O), q, lower.tail = !lower.tail, log.p = log.p)
        }
    }
    Z[!nat] <- z
    Z
}

#' @rdname d
#' @export
d.trans_contdist <- function(O, x, log = FALSE) {
    Z <- numeric(length(x))
    nat <- is.na(x)
    Z[nat] <- x[nat]
    x <- x[!nat]
    z <- numeric(length(x))
    on_support <- x <= sudo_support(O)["To"] & x >= sudo_support(O)["From"]
    x <- x[on_support]
    x1 <- eval(O$trafo$invtrans, list(X = x))
    if (log) {
        z[on_support] <- d(untrafo.trans_standist(O), x1, log = TRUE) + log(abs(eval(O$trafo$deriv, list(X = x))))
        z[!on_support] <- -Inf
    } else {
        z[on_support] <- d(untrafo.trans_standist(O), x1) * abs(eval(O$trafo$deriv, list(X = x)))
    }
    Z[!nat] <- z
    Z
}

#' @rdname d
#' @export
d.trans_discrdist <- function(O, x, log = FALSE) {
    Z <- numeric(length(x))
    nat <- is.na(x)
    Z[nat] <- x[nat]
    x <- x[!nat]
    z <- numeric(length(x))
    on_support <- x <= sudo_support(O)["To"] & x >= sudo_support(O)["From"]
    x <- x[on_support]
    x1 <- eval(O$trafo$invtrans, list(X = x))
    z[on_support] <- d(untrafo.trans_standist(O), x1, log = log)
    z[!on_support] <- if (log)
        -Inf else 0
    Z[!nat] <- z
    Z
}

#' @rdname q
#' @export
q.trans_univdist <- function(O, p, lower.tail = TRUE, log.p = FALSE, ...) {
    if (monot(O) == 1) {
        Z <- q(untrafo.trans_standist(O), p, lower.tail = lower.tail, log.p = log.p)
    } else {
        Z <- qlim(untrafo.trans_standist(O), p, lower.tail = !lower.tail, log.p = log.p)
    }
    eval(O$trafo$trans, list(X = Z))
}

#' @rdname r
#' @export
r.trans_univdist <- function(O, n) {
    D <- r(untrafo(O), n)
    eval(O$trafo$trans, list(X = D))
}

#' @rdname plim
#' @export
plim.trans_univdist <- function(O, q, lower.tail = TRUE, log.p = FALSE) {
    Z <- numeric(length(q))
    nat <- is.na(q)
    Z[nat] <- q[nat]
    q <- q[!nat]
    z <- numeric(length(q))
    r <- findInterval(q, sudo_support(O), left.open = TRUE)
    q <- eval(O$trafo$invtrans, list(X = q[r == 1]))
    if (lower.tail)
        z[r == 2] <- 1 else z[r == 0] <- 1
    if (log.p)
        z[r != 1] <- log(z[r != 1])
    if (length(q) > 0) {
        if (monot(O) == 1) {
            z[r == 1] <- plim(untrafo(O), q, lower.tail = lower.tail, log.p = log.p)
        } else {
            z[r == 1] <- p(untrafo(O), q, lower.tail = !lower.tail, log.p = log.p)
        }
    }
    Z[!nat] <- z
    Z
}

#' @rdname qlim
#' @export
qlim.trans_univdist <- function(O, p, lower.tail = TRUE, log.p = FALSE) {
    if (monot(O) == 1) {
        Z <- qlim(untrafo(O), p, lower.tail = lower.tail, log.p = log.p)
    } else {
        Z <- q(untrafo(O), p, lower.tail = !lower.tail, log.p = log.p)
    }
    eval(O$trafo$trans, list(X = Z))
}


