#' Fit Bayesian Univariate Box-Cox Transformed AR(1)
#' 
#' This is a function that fits a univariate Box-Cox model for energy usage
#' 
#' @author Daeyoung Lim, \email{daeyoung.lim@uconn.edu}
#' @param yobs a vector of measured energy usage to be used as the response variable
#' @param miss a vector of binary values (0 or 1), each element of which indicates whether the response variable is missing
#' @param a the shift parameter for the Box-Cox. By default, set to the minimum value minus one.
#' @param xobs the covariate matrix for the main regression equation
#' @param zobs the covariate matrix for the Box-Cox transformation parameter lambda
#' @param uobs the covariate matrix for the time-dependent variances sigma2
#' @param mcmc list of MCMC-related parameters: number of burn-ins (ndiscard), number of thinning(nskip), and posterior sample size (nkeep)
#' @param prior list of hyperparameters; when not given, algorithm will run in default setting
#' @param verbose logical variable for printing progress bar. Default to FALSE.
#' @return a dataframe with input arguments, posterior samples, Metropolis algorithm acceptance rates, etc
#' @examples
#' \dontrun{
#' data(df)
#' library(magrittr)
#' library(dplyr)
#' library(readr)
#' library(lubridate)
#' x <- df %>% filter(TimeStamp >= ymd_hms("2018-07-01 00:00:00"),
#'                    TimeStamp < ymd_hms("2018-08-01 00:00:00"))
#' y <- x$Reading %>% scale(center=FALSE) %>% drop
#' n <- length(y)
#' x1 <- x$TempRead %>% scale %>% drop
#' x2 <- x$HumidityRead %>% scale %>% drop
#' x3 <- 1*(hour(x$TimeStamp) %in% c(0:6, 19:23))
#' 
#' xobs <- cbind(1, x1, x2, x3)
#' zobs <- cbind(1, x1, x2)
#' uobs <- cbind(1, x1, x2, x3)
#' 
#' fit <- boxcoxar(y, integer(length(y)), mcmc=list(ndiscard=10000, nkeep = 20000),
#'          xobs = xobs, zobs = zobs, uobs = uobs, verbose=TRUE)
#' }
#' @export
boxcoxar <- function(yobs, miss, a, xobs, zobs, uobs, mcmc=list(), prior=list(), verbose = FALSE) {
  if (missing(a)) a <- min(yobs, na.rm=TRUE) - 1

  xcols <- ncol(xobs)
  zcols <- ncol(zobs)
  ucols <- ncol(uobs)

  privals <- list(mu_beta0 = numeric(xcols), sig_beta0 = diag(100, xcols),
                  mu_alpha0 = numeric(zcols), sig_alpha0 = diag(100, zcols),
                  mu_eta0 = numeric(ucols), sig_eta0 = diag(100, ucols))
  privals[names(prior)] <- prior
  mu_beta0 <- privals$mu_beta0
  sig_beta0 <- privals$sig_beta0
  mu_alpha0 <- privals$mu_alpha0
  sig_alpha0 <- privals$sig_alpha0
  mu_eta0 <- privals$mu_eta0
  sig_eta0 <- privals$sig_eta0

  mcvals <- list(ndiscard = 5000L, nskip = 1L, nkeep = 20000L)
  mcvals[names(mcmc)] = mcmc
  ndiscard <- mcvals$ndiscard
  nskip <- mcvals$nskip
  nkeep <- mcvals$nkeep


  mcmctime <- system.time({
    fout <- .Call(`boxcoxar_amh`,
                  PACKAGE="TSEnergy",
                  as.double(yobs),
                  as.integer(miss),
                  as.double(a),
                  as.matrix(xobs),
                  as.matrix(zobs),
                  as.matrix(uobs),
                  as.double(mu_beta0),
                  as.matrix(sig_beta0),
                  as.double(mu_alpha0),
                  as.matrix(sig_alpha0),
                  as.double(mu_eta0),
                  as.matrix(sig_eta0),
                  as.integer(ndiscard),
                  as.integer(nskip),
                  as.integer(nkeep),
                  as.logical(verbose))
  })
  out <- list(yobs = yobs,
              miss = miss,
              a = a,
              xobs = xobs,
              zobs = zobs,
              uobs = uobs,
              prior = privals,
              mcmctime = mcmctime,
              mcmc = mcvals,
              mcmc.draws = fout)
  class(out) <- "boxcoxar"
  out
}

