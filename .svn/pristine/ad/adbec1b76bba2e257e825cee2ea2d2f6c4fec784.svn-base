#' @title mistr: A Computational Framework for Univariate Mixture and Composite Distributions
#' @description A system offering object oriented handling of univariate distributions with focus on composite models.
#' @author Lukas Sablica, \email{lsablica@@wu.ac.at}
#'
#'         Kurt Hornik, \email{Kurt.Hornik@@wu.ac.at}
#'
#' \strong{Maintainer}: Lukas Sablica, \email{lsablica@@wu.ac.at}
#' @rdname mistr
#' @docType package
#' @name mistr-package
NULL
#' @rdname mistr
#' @docType package
#' @name mistr
#' @import stats
#' @import graphics
NULL

opt <- new.env(parent = emptyenv())
opt$sub <- 1e-10
opt$add <- 1e-8
opt$tol <- .Machine$double.eps^0.5

#' @title Get Parameters 
#' @description Function can be used to extract the parameters used in \code{\link{mistr}}.
#' @param ... characteristic strings of desired parameters. Possible values "sub", "add", "tol".
#' @return named vector with values.
#' @examples 
#' get_opt("sub", "tol")
#' @rdname get_opt
#' @seealso \code{\link{set_opt}}
#' @export 
get_opt <- function(...){
   dots <- list(...)
   filter.dots <- dots[dots %in% c("sub","add", "tol")]
   val <- unlist(lapply(filter.dots, function(x) opt[[x]]))
   names(val) <- filter.dots
   val
}


#' @title Creates New Distribution Object
#' @description The function creates distribution objects that satisfy the naming convention used in package mistr.
#' @param name string containing the name of the distribution.
#' @param from numeric representing where the support of distribution starts.
#' @param to numeric representing where the support of distribution ends.
#' @param by numeric representing the deterministic step between support values.
#'           If NULL: continuous distribution is assumed. If the value is specified: 
#'           discrete distribution with specified step is assumed, default: NULL.
#' @param parameters named list of parameters of the distribution, default: mget(names(eval(quote(match.call()),parent)[-1]),parent).
#' @param class class of the distribution, this should be set in [name]dist convention (e.g. normdist, tdist), 
#'              default: deparse(sys.calls()[[sys.nframe() - 1]][[1]]).
#' @param parent parent environment, default: parent.frame().
#' @return distribution object.
#' @details The function can be used in two ways. Either it can be called from the creator functions as for example 
#'          \code{\link{normdist}} or \code{\link{unifdist}}, or directly from any function or enviroment. In the former,
#'          only arguments "name", "from" and "to" must be set. Other arguments will be filled according to the parent calls.
#'           If this function is called directly, the arguments "parameters" and "class" have to be specified also. 
#' @examples 
#' \dontrun{
#' # using creator function
#' unifdist <- function(min = 0, max = 1) { 
#'    if (!is.numeric(min) || !is.numeric(max))   stop("Parameters must be a numeric")
#'    if (min >= max)   stop("min must be smaller than max.")
#'    new_dist(name = "Uniform", from = min, to = max)
#' }
#' 
#' #directly
#' U <- new_dist(name = "Uniform", from = 1, to = 6, 
#'               parameters =  list(min = 1, max = 6), class = "unifdist")
#' }
#' @rdname new_dist
#' @export 
new_dist <- function(name, from, to , by = NULL, parameters = mget(names(eval(quote(match.call()),parent)[-1]),parent),
                     class = deparse(sys.calls()[[sys.nframe() - 1]][[1]]), parent = parent.frame()){
   class2 <- if(is.null(by)) "contdist" else "discrdist"
   parameters <- lapply(parameters, function(a) unname(a))
   x <- list(parameters = parameters, type = name,
             support = list(from = from, to = to))
   x$support$by <- by
   class(x) <- c(class, class2, "standist", "univdist", "dist")
   return(x)
}

# dunifdist <- function(min = 1, max = 6){
#    if (!is.numeric(min) || !is.numeric(max)){
#       stop("parameters must be a numeric")
#    } 
#    if (min >= max){
#       stop("min must be smaller than max")
#    } 
#    if (min%%1 != 0 || max%%1 != 0){
#       stop("min and max must be integers")
#    } 
#    
#    new_dist(name = "Discrete uniform", from = min, to = max, by = 1)
# }

#' @title Set Parameters
#' @description Function can be used to set the parameters used in \code{\link{mistr}}.
#' @param ... arguments in tag = value form, or a list of tagged values. 
#' @return When parameters are set, their previous values are returned in an invisible named list. 
#' @details The function can set the values for:
#' 
#'          \strong{sub} parameter: small value that is used in mixture quantile function 
#'                                  to test if the computed value is infimum, default: 1e-10. 
#'                                  
#'          \strong{add} parameter: small value that is added to values that are in the image of CDF in \code{\link{qlim}}
#'                                  function, default: 1e-08.
#'                                  
#'          \strong{tol} parameter: tolerance for uniroot used in mixture quantile function, default: .Machine$double.eps^0.5.                        
#'          
#' @examples 
#' a <- set_opt(sub = 1e-5, tol = 1e-10)
#' get_opt("sub", "tol")
#' set_opt(a) 
#' @rdname set_opt
#' @export 
set_opt <- function(...){
   dots <- list(...)
   if (length(dots) == 1 && is.list(dots[[1L]])) dots <- dots[[1L]]
   filter.dots <- dots[names(dots) %in% c("sub","add", "tol")]
   old <- lapply(names(filter.dots), function(x) opt[[x]])
   names(old) <- names(filter.dots)
   lapply(names(filter.dots), function(x) opt[[x]] <- filter.dots[[x]])
   invisible(old)
}

#' @title Reports whether O is a Distribution Object
#' @description Reports whether O is a distribution object.
#' @param O an object to test.
#' @rdname is.dist
#' @export
is.dist <- function(O) {
    inherits(O, "dist")
}

#' @title Reports whether O is a Mixture Distribution Object
#' @description Reports whether O is a mixture distribution object.
#' @param O an object to test.
#' @rdname is.mixture
#' @export
is.mixture <- function(O) {
    inherits(O, "mixdist") || inherits(O, "trans_mixdist")
}

#' @title  Reports whether O is a Composite Distribution Object
#' @description Reports whether O is a composite distribution object.
#' @param O an object to test.
#' @rdname is.composite
#' @export
is.composite <- function(O) {
    inherits(O, "compdist") || inherits(O, "trans_compdist")
}

#' @title Reports whether O is a Standard Distribution Object
#' @description Reports whether O is a standard distribution object.
#' @param O an object to test.
#' @rdname is.standard
#' @export
is.standard <- function(O) {
    inherits(O, "standist") || inherits(O, "trans_standist")
}

#' @title Reports whether O is a Continuous Distribution Object
#' @description Reports whether O is a continuous distribution object.
#' @param O an object to test.
#' @rdname is.contin
#' @export
is.contin <- function(O) {
    inherits(O, "contdist") || inherits(O, "contmixdist") || inherits(O, "trans_contdist") || inherits(O, "trans_contmixdist") ||
        inherits(O, "contcompdist") || inherits(O, "trans_contcompdist")
}

#' @title Reports whether O is a Discrete Distribution Object
#' @description Reports whether O is a discrete distribution object.
#' @param O an object to test.
#' @rdname is.discrete
#' @export
is.discrete <- function(O) {
    inherits(O, "discrdist") || inherits(O, "discrmixdist") || inherits(O, "trans_discrdist") || inherits(O, "trans_discrmixdist") ||
        inherits(O, "discrcompdist") || inherits(O, "trans_discrcompdist")
}

#' @title Reports whether O is a Transformed Distribution Object
#' @description Reports whether O is a transformed distribution object.
#' @param O an object to test.
#' @rdname is.transformed
#' @export
is.transformed <- function(O) {
    inherits(O, "trans_univdist")
}


near <- function(x, y, tol = get_opt("tol")) {
    abs(x - y) < tol
}

#' @title Extract Model Parameters
#' @description \code{parameters} is a generic function which extracts parameters from \code{\link{mistr}} distribution objects.
#' @param O an object for which the extraction of model parameters is meaningful.
#' @return Vector (for standard distributions) or list (in the case of mixture/composite distribution)
#'         of parameters extracted from the object.
#'
#'         For a fitted object of class comp_fit returns vector of fitted parameters.
#' @examples
#' N <- normdist(1, 3)
#' parameters(N)
#'
#' C <- cauchydist()
#' M <- mixdist(N, C, weights = c(0.5, 0.5))
#' parameters(M)
#' @seealso \code{\link[stats]{weights}}, \code{\link{breakpoints}}
#' @rdname parameters
#' @export
parameters <- function(O) UseMethod("parameters")
#' @rdname parameters
#' @export
parameters.standist <- function(O) unlist(O$parameters)
#' @rdname parameters
#' @export
parameters.trans_standist <- function(O) unlist(O$parameters)
#' @rdname parameters
#' @export
parameters.mixdist <- function(O) lapply(O$objects, parameters)
#' @rdname parameters
#' @export
parameters.trans_mixdist <- function(O) lapply(O$objects, parameters)
#' @rdname parameters
#' @export
parameters.compdist <- function(O) lapply(O$objects, parameters)
#' @rdname parameters
#' @export
parameters.trans_compdist <- function(O) lapply(O$objects, parameters)
#' @rdname parameters
#' @export
parameters.comp_fit <- function(O) O$params$coeff


#' @export
weights.mixdist <- function(object, ...) object$weights
#' @export
weights.trans_mixdist <- function(object, ...) object$weights
#' @export
weights.compdist <- function(object, ...) object$weights
#' @export
weights.trans_compdist <- function(object, ...) object$weights
#' @export
weights.comp_fit <- function(object, ...) object$params$weights

#' @export
`[.mixdist` <- function(O, i, ...) O$objects[[i]]
#' @export
`[.trans_mixdist` <- function(O, i, ...) O$objects[[i]]
#' @export
`[.compdist` <- function(O, i, ...) O$objects[[i]]
#' @export
`[.trans_compdist` <- function(O, i, ...) O$objects[[i]]

#' @title Extract Model Breakpoints
#' @description \code{breakpoints} is a generic function which extracts breakpoints from \code{\link{mistr}} composite distribution objects.
#' @param O an object for which the extraction of model breakpoints is meaningful.
#' @return Vector of extracted breakpoints form object.
#' @seealso \code{\link{parameters}}, \code{\link[stats]{weights}}
#' @examples
#' N <- normdist(1, 3)
#' C <- cauchydist()
#'
#' CC <- compdist(N, C, weights = c(0.5, 0.5), breakpoints = 1)
#' breakpoints(CC)
#' @rdname breakpoints
#' @export
breakpoints <- function(O) UseMethod("breakpoints")
#' @rdname breakpoints
#' @export
breakpoints.compdist <- function(O) O$breakpoints
#' @rdname breakpoints
#' @export
breakpoints.trans_compdist <- function(O) sort(eval(O$trafo$trans, list(X = O$breakpoints)))
#' @rdname breakpoints
#' @export
breakpoints.comp_fit <- function(O) O$params$breakpoints

#' @title Extract Distribution of Fitted Model
#' @description \code{distribution} is a generic function which extracts the distribution with fitted parameters from fitted objects.
#' @param O an object for which the extraction of distribution is meaningful.
#' @return Object representing the distribution.
#' @rdname distribution
#' @export
distribution <- function(O) UseMethod("distribution")
#' @rdname distribution
#' @export
distribution.comp_fit <- function(O) O$Distribution

#' @title Displays a Useful Description of a Distribution Object
#' @description Displays a useful description of a distribution object from \code{\link{mistr}}.
#' @param object distribution object to summarize.
#' @param level adds 3*(level-1) spaces before the print, default: 1.
#' @param space number of blank lines between outputs, default: 2.
#' @param additional_list,truncation,... additional information that may be passed to summary.
#' @details \code{summary} prints useful description of a distribution object. This feature might
#'          be useful when working with a more complicated distribution that contains
#'          mixture and composite distributions as components and the print function does not
#'          offer enough information.
#'
#'          Arguments \code{level}, \code{additional_list} and truncation
#'          are present for recursive usage that is done for more complicated models
#'          automatically by the function.
#' @rdname summary
#' @name Distribution_summary
NULL
#' @rdname summary
#' @method summary standist
#' @export
summary.standist <- function(object, level = 1, space = 2, additional_list, truncation, ...) {
    if (missing(additional_list)) {
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), object$type, " Distribution:"), "\n")
        cat(sep = "", paste(rep("   ", level - 1), collapse = ""), paste(rep("-", nchar(object$type) + 14), collapse = ""), " \n")
        Parameters <- paste(names(object$parameters), unlist(object$parameters), sep = " = ", collapse = ", ")
        Support <- paste(c("From", "To"), c(object$support$from, object$support$to), sep = ": ", collapse = ", ")
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Parameters: ", Parameters), rep("\n", space))
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Support: ", Support, sep = ""), "\n")
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "_____________________________", sep = ""), rep("\n", space))
    } else {
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "[", additional_list$n, "] ", object$type, " Distribution:"),
            "\n")
        cat(paste(rep("   ", level), collapse = ""), paste(rep("-", nchar(object$type) + 14), collapse = ""), " \n")
        Parameters <- paste(names(object$parameters), unlist(object$parameters), sep = " = ", collapse = ", ")
        Support <- paste(c("From", "To"), c(object$support$from, object$support$to), sep = ": ", collapse = ", ")
        cat(paste0(paste(rep("   ", level), collapse = ""), " Parameters: ", Parameters), rep("\n", space))
        cat(paste0(paste(rep("   ", level), collapse = ""), " Support: ", Support, sep = ""), rep("\n", space))
        if (!missing(truncation))
            cat(paste0(paste(rep("   ", level), collapse = ""), " Truncated to: ", truncation), rep("\n", space))
        cat(paste0(paste(rep("   ", level), collapse = ""), " Weight in ", additional_list$model, " model = ", round(additional_list$prob,
            4), ", Overall weight in model = ", round(additional_list$cumprob, 4), sep = ""), "\n")
        cat(paste0(paste(rep("   ", level), collapse = ""), " __________________________________________________________________",
            sep = ""), rep("\n", space))
    }
    invisible(object)
}

#' @rdname summary
#' @method summary trans_standist
#' @export
summary.trans_standist <- function(object, level = 1, space = 2, additional_list, truncation, ...) {
    if (missing(additional_list)) {
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), object$type, " Distribution:"), "\n")
        cat(sep = "", paste(rep("   ", level - 1), collapse = ""), paste(rep("-", nchar(object$type) + 14), collapse = ""), " \n")
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Transformation: ", deparse(object$trafo$print)), rep("\n", space))
        Parameters <- paste(names(object$parameters), unlist(object$parameters), sep = " = ", collapse = ", ")
        Support <- paste(c("From", "To"), sudo_support(object), sep = ": ", collapse = ", ")
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Parameters: ", Parameters), rep("\n", space))
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Support: ", Support, sep = ""), "\n")
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "_____________________________", sep = ""), rep("\n", space))

    } else {
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "[", additional_list$n, "] ", object$type, " Distribution:"),
            "\n")
        cat(paste(rep("   ", level), collapse = ""), paste(rep("-", nchar(object$type) + 14), collapse = ""), " \n")
        cat(paste0(paste(rep("   ", level), collapse = ""), " Transformation: ", deparse(object$trafo$print)), rep("\n", space))
        Parameters <- paste(names(object$parameters), unlist(object$parameters), sep = " = ", collapse = ", ")
        Support <- paste(c("From", "To"), sudo_support(object), sep = ": ", collapse = ", ")
        cat(paste0(paste(rep("   ", level), collapse = ""), " Parameters: ", Parameters), rep("\n", space))
        cat(paste0(paste(rep("   ", level), collapse = ""), " Support: ", Support, sep = ""), rep("\n", space))
        if (!missing(truncation))
            cat(paste0(paste(rep("   ", level), collapse = ""), " Truncated to: ", truncation), rep("\n", space))
        cat(paste0(paste(rep("   ", level), collapse = ""), " Weight in ", additional_list$model, " model = ", round(additional_list$prob,
            4), ", Overall weight in model = ", round(additional_list$cumprob, 4), sep = ""), " \n")
        cat(paste0(paste(rep("   ", level), collapse = ""), " __________________________________________________________________",
            sep = ""), rep("\n", space))
    }
    invisible(object)
}


#' @rdname summary
#' @method summary mixdist
#' @export
summary.mixdist <- function(object, level = 1, space = 2, additional_list, truncation, ...) {
    g <- object$weights
    if (missing(additional_list)) {
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Mixture Distribution:"), "\n")
        cat(sep = "", paste(rep("   ", level - 1), collapse = ""), paste(rep("-", 21), collapse = ""), " \n")
        Support <- paste(c("From", "To"), sudo_support(object), sep = ": ", collapse = ", ")
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Support: ", Support), rep("\n", space))
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Components:"), rep("\n", space))
        n <- lapply(seq_along(object$objects), function(i) summary(object$objects[[i]], level = level + 1, space = space, additional_list = list(n = i,
            prob = g[i], cumprob = g[i], model = "Mixture")))
    } else {
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "[", additional_list$n, "] ", "Mixture Distribution:"), "\n")
        cat(paste(rep("   ", level), collapse = ""), paste(rep("-", 21), collapse = ""), " \n")
        cat(paste0(paste(rep("   ", level), collapse = ""), " Support: ", paste(c("From", "To"), sudo_support(object), sep = ": ",
            collapse = ", ")), rep("\n", space))
        if (!missing(truncation))
            cat(paste0(paste(rep("   ", level), collapse = ""), " Truncated to: ", truncation), rep("\n", space))
        cat(paste0(paste(rep("   ", level), collapse = ""), " Components:"), rep("\n", space))
        n <- lapply(seq_along(object$objects), function(i) summary(object$objects[[i]], level = level + 1, space = space, additional_list = list(n = i,
            prob = g[i], cumprob = g[i] * additional_list$cumprob, model = "Mixture")))
        cat(paste0(paste(rep("   ", level), collapse = ""), " Weight in ", additional_list$model, " model = ", round(additional_list$prob,
            4), ", Overall weight in model = ", round(additional_list$cumprob, 4), sep = ""), " \n")
        cat(paste0(paste(rep("   ", level), collapse = ""), " _______________________________________________________________________",
            sep = ""), rep("\n", space))
    }
    invisible(object)
}

#' @rdname summary
#' @method summary trans_mixdist
#' @export
summary.trans_mixdist <- function(object, level = 1, space = 2, additional_list, truncation, ...) {
    g <- object$weights
    if (missing(additional_list)) {
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Mixture Distribution:"), "\n")
        cat(sep = "", paste(rep("   ", level - 1), collapse = ""), paste(rep("-", 21), collapse = ""), " \n")
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Transformation: ", deparse(object$trafo$print)), rep("\n", space))
        Support <- paste(c("From", "To"), sudo_support(object), sep = ": ", collapse = ", ")
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Support: ", Support), rep("\n", space))
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Components:"), rep("\n", space))
        n <- lapply(seq_along(object$objects), function(i) summary(object$objects[[i]], level = level + 1, space = space, additional_list = list(n = i,
            prob = g[i], cumprob = g[i], model = "Mixture")))
    } else {
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "[", additional_list$n, "] ", "Mixture Distribution:"), "\n")
        cat(paste(rep("   ", level), collapse = ""), paste(rep("-", 21), collapse = ""), " \n")
        cat(paste0(paste(rep("   ", level), collapse = ""), " Transformation: ", deparse(object$trafo$print)), rep("\n", space))
        cat(paste0(paste(rep("   ", level), collapse = ""), " Support: ", paste(c("From", "To"), sudo_support(object), sep = ": ",
            collapse = ", ")), rep("\n", space))
        if (!missing(truncation))
            cat(paste0(paste(rep("   ", level), collapse = ""), " Truncated to: ", truncation), rep("\n", space))
        cat(paste0(paste(rep("   ", level), collapse = ""), " Components:"), rep("\n", space))
        n <- lapply(seq_along(object$objects), function(i) summary(object$objects[[i]], level = level + 1, space = space, additional_list = list(n = i,
            prob = g[i], cumprob = g[i] * additional_list$cumprob, model = "Mixture")))
        cat(paste0(paste(rep("   ", level), collapse = ""), " Weight in ", additional_list$model, " model = ", round(additional_list$prob,
            4), ", Overall weight in model = ", round(additional_list$cumprob, 4), sep = ""), " \n")
        cat(paste0(paste(rep("   ", level), collapse = ""), " _______________________________________________________________________",
            sep = ""), rep("\n", space))
    }
    invisible(object)
}


#' @rdname summary
#' @method summary compdist
#' @export
summary.compdist <- function(object, level = 1, space = 2, additional_list, truncation, ...) {
    g <- object$weights
    if (missing(additional_list)) {
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Composite Distribution:"), "\n")
        cat(sep = "", paste(rep("   ", level - 1), collapse = ""), paste(rep("-", 23), collapse = ""), " \n")
        Support <- paste(c("From", "To"), sudo_support(object), sep = ": ", collapse = ", ")
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Support: ", Support), rep("\n", space))
        vv <- unlist(lapply(object$interval, function(x) if (x == "R")
            c(")", "[") else c("]", "(")))
        b <- paste(c("(", vv[seq.int(1L, length(vv), 2L) + 1]), c("-Inf", object$breakpoints), ",", c(object$breakpoints, "Inf"), c(vv[seq.int(1L,
            length(vv), 2L)], ")"), sep = "")
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Components:"), rep("\n", space))
        n <- lapply(seq_along(object$objects), function(i) summary(object$objects[[i]], level = level + 1, space = space, additional_list = list(n = i,
            prob = g[i], cumprob = g[i], model = "Composite"), truncation = b[i]))
    } else {
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "[", additional_list$n, "] ", "Composite Distribution:"),
            "\n")
        cat(paste(rep("   ", level), collapse = ""), paste(rep("-", 23), collapse = ""), " \n")
        cat(paste0(paste(rep("   ", level), collapse = ""), " Support: ", paste(c("From", "To"), sudo_support(object), sep = ": ",
            collapse = ", ")), rep("\n", space))
        if (!missing(truncation))
            cat(paste0(paste(rep("   ", level), collapse = ""), " Truncated to: ", truncation), rep("\n", space))
        cat(paste0(paste(rep("   ", level), collapse = ""), " Components:"), rep("\n", space))
        vv <- unlist(lapply(object$interval, function(x) if (x == "R")
            c(")", "[") else c("]", "(")))
        b <- paste(c("(", vv[seq.int(1L, length(vv), 2L) + 1]), c("-Inf", object$breakpoints), ",", c(object$breakpoints, "Inf"), c(vv[seq.int(1L,
            length(vv), 2L)], ")"), sep = "")
        n <- lapply(seq_along(object$objects), function(i) summary(object$objects[[i]], level = level + 1, space = space, additional_list = list(n = i,
            prob = g[i], cumprob = g[i] * additional_list$cumprob, model = "Composite"), truncation = b[i]))
        cat(paste0(paste(rep("   ", level), collapse = ""), " Weight in ", additional_list$model, " model = ", round(additional_list$prob,
            4), ", Overall weight in model = ", round(additional_list$cumprob, 4), sep = ""), " \n")
        cat(paste0(paste(rep("   ", level), collapse = ""), " _______________________________________________________________________",
            sep = ""), rep("\n", space))
    }
    invisible(object) 
}

#' @rdname summary
#' @method summary trans_compdist
#' @export
summary.trans_compdist <- function(object, level = 1, space = 2, additional_list, truncation, ...) {
    g <- object$weights
    if (missing(additional_list)) {
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Composite Distribution:"), "\n")
        cat(sep = "", paste(rep("   ", level - 1), collapse = ""), paste(rep("-", 23), collapse = ""), " \n")
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Transformation: ", deparse(object$trafo$print)), rep("\n", space))
        Support <- paste(c("From", "To"), sudo_support(object), sep = ": ", collapse = ", ")
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Support: ", Support), rep("\n", space))
        vv <- unlist(lapply(object$interval, function(x) if (x == "R")
            c(")", "[") else c("]", "(")))
        b <- paste(c("(", vv[seq.int(1L, length(vv), 2L) + 1]), c("-Inf", object$breakpoints), ",", c(object$breakpoints, "Inf"), c(vv[seq.int(1L,
            length(vv), 2L)], ")"), sep = "")
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "Components:"), rep("\n", space))
        n <- lapply(seq_along(object$objects), function(i) summary(object$objects[[i]], level = level + 1, space = space, additional_list = list(n = i,
            prob = g[i], cumprob = g[i], model = "Composite"), truncation = b[i]))
    } else {
        cat(paste0(paste(rep("   ", level - 1), collapse = ""), "[", additional_list$n, "] ", "Composite Distribution:"),
            "\n")
        cat(paste(rep("   ", level), collapse = ""), paste(rep("-", 23), collapse = ""), " \n")
        cat(paste0(paste(rep("   ", level), collapse = ""), " Transformation: ", deparse(object$trafo$print)), rep("\n", space))
        cat(paste0(paste(rep("   ", level), collapse = ""), " Support: ", paste(c("From", "To"), sudo_support(object), sep = ": ",
            collapse = ", ")), rep("\n", space))
        if (!missing(truncation))
            cat(paste0(paste(rep("   ", level), collapse = ""), " Truncated to: ", truncation), rep("\n", space))
        cat(paste0(paste(rep("   ", level), collapse = ""), " Components:"), rep("\n", space))
        vv <- unlist(lapply(object$interval, function(x) if (x == "R")
            c(")", "[") else c("]", "(")))
        b <- paste(c("(", vv[seq.int(1L, length(vv), 2L) + 1]), c("-Inf", object$breakpoints), ",", c(object$breakpoints, "Inf"), c(vv[seq.int(1L,
            length(vv), 2L)], ")"), sep = "")
        n <- lapply(seq_along(object$objects), function(i) summary(object$objects[[i]], level = level + 1, space = space, additional_list = list(n = i,
            prob = g[i], cumprob = g[i] * additional_list$cumprob, model = "Composite"), truncation = b[i]))
        cat(paste0(paste(rep("   ", level), collapse = ""), " Weight in ", additional_list$model, " model = ", round(additional_list$prob,
            4), ", Overall weight in model = ", round(additional_list$cumprob, 4), sep = ""), " \n")
        cat(paste0(paste(rep("   ", level), collapse = ""), " _______________________________________________________________________",
            sep = ""), rep("\n", space))
    }
    invisible(object)
}

#' @title Displays a Useful Description of a Fitted Object
#' @description Displays a useful description of a fitted object.
#' @param object distribution object to summarize.
#' @param ...	 additional arguments.
#' @return Function returns summary of the fit, offered by bbmle package for class \code{\link[bbmle]{mle2-class}}.
#' @rdname summary_comp_fit
#' @seealso \code{\link[bbmle]{mle2-class}}
#' @export
#' @method summary comp_fit
#' @importFrom bbmle summary
summary.comp_fit <- function(object, ...) bbmle::summary(object$fit)




