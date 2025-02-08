# empPvals function
empPvals <- function(stat, stat0) {
    m <- length(stat)
    m0 <- length(stat0)
    if (is.matrix(stat0)) 
      stat0 <- as.vector(stat0)
    v <- c(rep(TRUE, m), rep(FALSE, m0))
    v <- v[order(c(stat, stat0), decreasing = TRUE)]
    u <- 1:length(v)
    w <- 1:m
    p <- (u[v == TRUE] - w)/m0
    p <- p[rank(-stat)]
    p <- pmax(p, 1/m0)
    if (anyNA(stat)) 
      p[is.na(stat)] <- NA
    return(p)
}

# FSTAT function
FSTAT <- function (dat, LV, ALV = NULL, covariate = NULL, parametric = FALSE) 
{
  if (missing(dat)) 
    stop("`dat` is required!")
  if (!is.matrix(dat)) 
    stop("`dat` must be a matrix!")
  if (missing(LV)) 
    stop("`LV` is required!")
  if (!is.matrix(LV)) 
    stop("`LV` must be a matrix!")
  m <- nrow(dat)
  n <- ncol(dat)
  if (nrow(LV) != n) 
    stop("`LV` number of rows (", nrow(LV), ") does not equal the number of columns of `dat` (", 
         n, ")!")
  if (is.null(ALV)) {
    if (is.null(covariate)) {
      model.alt <- stats::model.matrix(seq(n) ~ LV)
      model.null <- stats::model.matrix(seq(n) ~ 1)
    }
    if (!is.null(covariate)) {
      model.alt <- stats::model.matrix(seq(n) ~ LV + covariate)
      model.null <- stats::model.matrix(seq(n) ~ 1 + covariate)
    }
  }
  else if (is.matrix(ALV) || is.vector(ALV)) {
    if (is.null(covariate)) {
      model.alt <- stats::model.matrix(seq(n) ~ LV + ALV)
      model.null <- stats::model.matrix(seq(n) ~ 1 + ALV)
    }
    if (!is.null(covariate)) {
      model.alt <- stats::model.matrix(seq(n) ~ LV + ALV + 
                                         covariate)
      model.null <- stats::model.matrix(seq(n) ~ 1 + ALV + 
                                          covariate)
    }
  }
  else stop("Invalid arguments into a function 'FSTAT'. Adjustment latent variable must be either a matrix (n rows) or a vector (size n)")
  RSS.alt <- RSS(dat, model.alt)
  RSS.null <- RSS(dat, model.null)
  fstat <- (RSS.null - RSS.alt)/(ncol(model.alt) - ncol(model.null))/(RSS.alt/(n - 
                                                                                 ncol(model.alt)))
  if (!parametric) {
    return(list(fstat = fstat))
  }
  else {
    fstat.pval <- 1 - stats::pf(fstat, ncol(model.alt) - 
                                  ncol(model.null), n - ncol(model.alt))
    return(list(fstat = fstat, p.value = fstat.pval))
  }
}

# RSS function
RSS <- function (dat, mod) 
{
  if (missing(dat)) 
    stop("`dat` is required!")
  if (missing(mod)) 
    stop("`mod` is required!")
  if (!is.matrix(mod)) 
    stop("`mod` must be a matrix!")
  if (is.matrix(dat)) {
    n <- ncol(dat)
  }
  else if (is.vector(dat)) {
    n <- length(dat)
  }
  else stop("`dat` must be vector or matrix!")
  Id <- diag(n)
  res <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
  rss <- rowSums(res^2)
  return(rss)
}