binomialNoWarn<-function (link = "logit") 
{
    linktemp <- substitute(link)
    if (!is.character(linktemp)) {
        linktemp <- deparse(linktemp)
        if (linktemp == "link") 
            linktemp <- eval(link)
    }
    if (any(linktemp == c("logit", "probit", "cloglog", "cauchit", 
        "log"))) 
        stats <- make.link(linktemp)
    else stop(gettextf("link \"%s\" not available for binomial family, available links are \"logit\", \"\"probit\", \"cloglog\", \"cauchit\" and \"log\"", 
        linktemp), domain = NA)
    variance <- function(mu) mu * (1 - mu)
    validmu <- function(mu) all(mu > 0) && all(mu < 1)
    dev.resids <- function(y, mu, wt) .Call("binomial_dev_resids", 
        y, mu, wt, PACKAGE = "stats")
    aic <- function(y, n, mu, wt, dev) {
        m <- if (any(n > 1)) 
            n
        else wt
        -2 * sum(ifelse(m > 0, (wt/m), 0) * dbinom(round(m * 
            y), round(m), mu, log = TRUE))
    }
    initialize <- expression({
        if (NCOL(y) == 1) {
            if (is.factor(y)) y <- y != levels(y)[1]
            n <- rep.int(1, nobs)
            if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
            mustart <- (weights * y + 0.5)/(weights + 1)
            m <- weights * y
        } else if (NCOL(y) == 2) {
            n <- y[, 1] + y[, 2]
            y <- ifelse(n == 0, 0, y[, 1]/n)
            weights <- weights * n
            mustart <- (n * y + 0.5)/(n + 1)
        } else stop("for the binomial family, y must be a vector of 0 and 1's\n", 
            "or a 2 column matrix where col 1 is no. successes and col 2 is no. failures")
    })
    structure(list(family = "binomial", link = linktemp, linkfun = stats$linkfun, 
        linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, 
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
        validmu = validmu, valideta = stats$valideta), class = "family")
}
