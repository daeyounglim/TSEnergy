#' Fit Bayesian Multivariate Box-Cox Transformed AR(1)
#' 
#' This is a function that fits a univariate Box-Cox model for energy usage
#' 
#' @author Daeyoung Lim, \email{daeyoung.lim@uconn.edu}
#' @param yobs a matrix of measured energy usage to be used as the response variable; each row is a different time series of length T
#' @param miss a integer matrix of binary values (0 or 1), each element of which indicates whether the response variable is missing
#' @param a the shift parameter for the Box-Cox. By default, set to the minimum value minus one.
#' @param xobs the covariate matrix for the systemic equations for both the regression mean and Box-Cox lambda
#' @param mcmc list of MCMC-related parameters: number of burn-ins (ndiscard), number of thinning(nskip), and posterior sample size (nkeep)
#' @param prior list of hyperparameters; when not given, algorithm will run in default setting
#' @param verbose logical variable for printing progress bar. Default to FALSE.
#' @return a dataframe with input arguments, posterior samples, Metropolis algorithm acceptance rates, etc
#' @examples
#' \dontrun{
#'
#' y<-matrix(rnorm(100*3),3,100)
#' x1<-rnorm(100)
#' x2<-rnorm(100)
#' 
#' fit <- mvboxcoxar(y, matrix(0L, 3, 100), xobs = cbind(x1, x2), mcmc=list(ndiscard=10000, nkeep = 20000),
#'          xobs = xobs, zobs = zobs, uobs = uobs, verbose=TRUE)
#' }
#' @export
mvboxcoxar <- function(yobs, miss, a, xobs, mcmc=list(), prior=list(), verbose = FALSE) {
  if (missing(a)) a <- min(yobs, na.rm=TRUE) - 1

  xcols <- ncol(xobs)

  privals <- list(sigb2 = 10, sigc2 = 10)
  privals[names(prior)] <- prior
  sigb2 <- privals$sigb2
  sigc2 <- privals$sigc2

  mcvals <- list(ndiscard = 5000L, nskip = 1L, nkeep = 20000L)
  mcvals[names(mcmc)] = mcmc
  ndiscard <- mcvals$ndiscard
  nskip <- mcvals$nskip
  nkeep <- mcvals$nkeep


  mcmctime <- system.time({
    fout <- .Call('_TSEnergy_boxcoxvar_amh',
                  as.matrix(yobs),
                  as.matrix(miss),
                  as.double(a),
                  as.matrix(xobs),
                  as.double(sigb2),
                  as.double(sigc2),
                  as.integer(ndiscard),
                  as.integer(nskip),
                  as.integer(nkeep),
                  as.logical(verbose))
  })
  out <- list(yobs = yobs,
              miss = miss,
              a = a,
              xobs = xobs,
              prior = privals,
              mcmctime = mcmctime,
              mcmc = mcvals,
              mcmc.draws = fout)
  class(out) <- "mvboxcoxar"
  out
}

