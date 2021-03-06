% 2013-09-06 A. Papritz
% R CMD Rdconv -t html -o bla.html georob.control.Rd ; open bla.html; R CMD Rd2pdf --force georob.control.Rd; 
\encoding{macintosh}
\name{georob.control}
\alias{georob.control}
\alias{param.transf}
\alias{fwd.transf}
\alias{dfwd.transf}
\alias{bwd.transf}
\alias{rq.control}
% \alias{bbsolve.control}
\alias{nleqslv.control}
\alias{optim.control}

\title{Tuning Parameters for georob}

\description{
  This page documents the tuning parameters for \code{\link{georob}}.  It
  describes the arguments of the functions \code{control.georob},
  \code{param.transf}, \code{fwd.transf}, \code{dfwd.transf},
  \code{bwd.transf}, \code{rq.control}, \code{nleqslv.control} and
  \code{optim.control}, which all serve to control the behaviour of
  \code{\link{georob}}.}

\usage{
georob.control(initial.method = c("lmrob", "rq", "lm"), bhat = NULL, 
    param.tf = param.transf(),fwd.tf = fwd.transf(), 
    deriv.fwd.tf = dfwd.transf(), bwd.tf = bwd.transf(), 
    safe.param = 1.e12, psi.func = c("logistic", "t.dist", "huber"),
    tuning.psi.nr = 1000, min.rweight = 0.25,
    irwls.initial = TRUE, irwls.maxiter = 50, 
    irwls.reltol = sqrt(.Machine[["double.eps"]]),
    force.gradient = FALSE, zero.dist = sqrt(.Machine[["double.eps"]]),
    cov.bhat = FALSE, full.cov.bhat = FALSE, cov.betahat = TRUE, 
    cov.bhat.betahat = FALSE, 
    cov.delta.bhat = TRUE, full.cov.delta.bhat = TRUE,
    cov.delta.bhat.betahat = TRUE,
    cov.ehat = TRUE, full.cov.ehat = FALSE,
    cov.ehat.p.bhat = FALSE, full.cov.ehat.p.bhat = FALSE,
    aux.cov.pred.target = FALSE, min.condnum = 1.e-12,
    rq = rq.control(), lmrob = lmrob.control(),
    nleqslv = nleqslv.control(), 
%     bbsolve = bbsolve.control(), 
    optim = optim.control(), full.output = TRUE)
    
param.transf( variance = "log", snugget = "log", nugget = "log", scale = "log", 
    alpha = "identity", beta = "log", delta = "identity", 
    gamma = "identity", kappa = "identity", lambda = "identity", 
    mu = "log", nu = "log",
    f1 = "log", f2  ="log", omega = "rad", phi = "rad", zeta = "rad"
  )
  
fwd.transf(...)

dfwd.transf(...)

bwd.transf(...)

rq.control(tau = 0.5, rq.method = "br", rq.alpha = 0.1, ci = FALSE, iid = TRUE, 
    interp = TRUE, tcrit = TRUE, rq.beta = 0.99995, eps = 1e-06, 
    Mm.factor = 0.8, max.bad.fixup = 3)
               
nleqslv.control(nleqslv.method = c("Broyden", "Newton"), 
    global = c( "dbldog", "pwldog", "qline", "gline", "none" ),
    xscalm = c( "fixed", "auto" ), nleqslv.control = NULL)
    
% bbsolve.control(bbsolve.method = c("2", "3", "1"), bbsolve.control = NULL)
    
optim.control(optim.method = c("BFGS", "Nelder-Mead", "CG", 
        "L-BFGS-B", "SANN", "Brent"), lower = -Inf, upper = Inf, 
    optim.control = NULL, hessian = TRUE)
 

}

\arguments{

  \item{initial.method}{character keyword defining whether the function
  \code{\link[robustbase]{lmrob}} or \code{\link[quantreg]{rq}} is used to
  compute robust initial estimates of the regression parameters
  \eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta} (default \code{"lmrob"}).
  If the fixed effects model matrix has not full columns rank, then
  \code{\link[stats]{lm}} is used to compute initial values of the
  regression coefficients.  Note that this requires robust estimation.}
  
  \item{bhat}{initial values for the spatial random effects
  \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB}, with
  \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}=\mbox{\boldmath$0$\unboldmath}}{hatB=0}
  if \code{bhat} is equal to \code{NULL} (default).}
  
  \item{param.tf}{a function such as \code{param.transf}, which returns a
  named vector of character strings that define the transformations to be
  applied to the variogram parameters for model fitting, see
  \emph{Details}.}
  
  \item{fwd.tf}{a function such as \code{fwd.transf}, which returns a named
  list of invertible functions to be used to transform variogram
  parameters, see \emph{Details}.}

  \item{deriv.fwd.tf}{a function such as \code{dfwd.transf}, which
  returns a named list of functions corresponding to the first derivatives
  of \code{fwd.tf}, see \emph{Details}.}

  \item{bwd.tf}{a function such as \code{bwd.transf}, which returns the
  named list of inverse functions corresponding to \code{fwd.tf}, see
  \emph{Details}.}

  \item{safe.param}{maximum acceptable value for any variogram parameter.
  If trial parameter values generated by \code{\link{optim}} or
  \code{\link{nleqslv}} exceed \code{safe.param} then an error is signalled
  to force \code{optim} or \code{nleqslv} to update the trial values
  (default \code{1.e12}).}

  \item{psi.func}{character keyword defining what \eqn{\psi_c}-function should be
  used for robust model fitting.  Possible values are \code{"logistic"} (a
  scaled and shifted logistic cdf, default), \cr \code{"t.dist"} (re-descending
  \eqn{\psi_c}-function associated with Student \eqn{t}-distribution with
  \eqn{c} degrees of freedom) and \code{"huber"} (Huber's
  \eqn{\psi_c}-function).}

  \item{tuning.psi.nr}{positive numeric.  If \code{tuning.psi} is less than
  \code{tuning.psi.nr} then the model is fitted robustly by solving the
  robustified estimating equations, and for \code{tuning.psi} equal to or
  larger than \code{tuning.psi.nr} the Gaussian restricted loglikelihood is
  maximized (default \code{1000}).}
  
  \item{min.rweight}{positive numeric.  \dQuote{Robustness weight} of
  the initial \code{\link[robustbase]{lmrob}} fit that observations must
  exceed to be used for computing robust initial estimates of variogram
  parameters by setting \code{initial.param = TRUE} (see
  \code{\link{georob}}; default \code{0.25}).}
  
  \item{irwls.initial}{logical.  If \code{TRUE} (default) the estimating
  equations of \eqn{\mbox{\boldmath$B$\unboldmath}}{B} and
  \eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta} are always solved by
  IRWLS from the initial estimates of
  \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB} and
  \eqn{\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\beta}.  If
  \code{FALSE} then IRWLS starts from respective estimates computed for the
  variogram parameter estimates of the previous iteration of \code{nleqslv}
  or \code{optim}.}
  
  \item{irwls.maxiter}{positive integer equal to the maximum number of
  IRWLS iterations to solve the estimating equations of
  \eqn{\mbox{\boldmath$B$\unboldmath}}{B} and
  \eqn{\mbox{\boldmath$\beta$\unboldmath}}{\beta} (default \code{50}).}
   
  \item{irwls.reltol}{numeric convergence criterion for IRWLS.
  Convergence is assumed if\cr
  \code{max(abs(oldres-newres)) < sqrt(reltol*nugget)}, 
  where \code{oldres} and \cr \code{newres} are the
  estimated residuals \eqn{\widehat{\varepsilon}_i}{\epsilon_i} of the
  previous and current iterations, respectively.}
  
  \item{force.gradient}{logical controlling whether the estimating
  equations or the gradient of the Gaussian restricted loglikelihood are
  evaluated even if all variogram parameters are fixed (default
  \code{FALSE}).}

  \item{zero.dist}{positive numeric equal to the maximum distance, separating two
  sampling locations that are still considered as being coincident.}
  
  \item{cov.bhat}{logical controlling whether the covariances of
  \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB} are returned by
  \code{georob} (default \code{FALSE}).}
  
  \item{full.cov.bhat}{logical controlling whether the full covariance
  matrix (\code{TRUE}) or only the variance vector of
  \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB} is returned (default
  \code{FALSE}).}
  
  \item{cov.betahat}{logical controlling whether the covariance matrix of
  \eqn{\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\beta} is returned
  (default \code{TRUE}).}
  
  \item{cov.bhat.betahat}{logical controlling whether the covariance matrix
  of \eqn{\widehat{\mbox{\boldmath$B$\unboldmath}}}{hatB} and
  \eqn{\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\beta} is returned
  (default \code{FALSE}).}
  
  \item{cov.delta.bhat}{logical controlling whether the covariances of
  \eqn{\mbox{\boldmath$B$\unboldmath}-
  \widehat{\mbox{\boldmath$B$\unboldmath}}}{B-hatB} are returned (default \code{TRUE}).}
  
  \item{full.cov.delta.bhat}{logical controlling whether the full covariance
  matrix (\code{TRUE}) or only the variance vector of
  \eqn{\mbox{\boldmath$B$\unboldmath}-
  \widehat{\mbox{\boldmath$B$\unboldmath}}}{B-hatB} is returned (default \code{TRUE}).}
  
  \item{cov.delta.bhat.betahat}{logical controlling whether the covariance
  matrix of \eqn{\mbox{\boldmath$B$\unboldmath}-
  \widehat{\mbox{\boldmath$B$\unboldmath}}}{B-hatB} and
  \eqn{\widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\beta} is returned
  (default \code{TRUE}).}
  
  \item{cov.ehat}{logical controlling whether the covariances of
  \eqn{\widehat{\mbox{\boldmath$\varepsilon$\unboldmath}}=\mbox{\boldmath$Y$\unboldmath}
  - \mbox{\boldmath$X$\unboldmath}
  \widehat{\mbox{\boldmath$\beta$\unboldmath}} -
  \widehat{\mbox{\boldmath$B$\unboldmath}}}{hat\epsilon=Y-X hat\beta -
  hatB} are returned (default \code{TRUE}).}
  
  \item{full.cov.ehat}{logical controlling whether the full covariance
  matrix (\code{TRUE}) or only the variance vector of
  \eqn{\widehat{\mbox{\boldmath$\varepsilon$\unboldmath}}=\mbox{\boldmath$Y$\unboldmath}
  - \mbox{\boldmath$X$\unboldmath}
  \widehat{\mbox{\boldmath$\beta$\unboldmath}} -
  \widehat{\mbox{\boldmath$B$\unboldmath}}}{hat\epsilon=Y-X hat\beta -
  hatB} is returned (default \code{FALSE}).}
  
  \item{cov.ehat.p.bhat}{logical controlling whether the covariances of
  \eqn{\widehat{\mbox{\boldmath$\varepsilon$\unboldmath}} +
  \widehat{\mbox{\boldmath$B$\unboldmath}} =\mbox{\boldmath$Y$\unboldmath}
  - \mbox{\boldmath$X$\unboldmath}
  \widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\epsilon+ hatB=Y-X
  hat\beta} are returned (default \code{FALSE}).}
  
  \item{full.cov.ehat.p.bhat}{logical controlling whether the full
  covariance matrix (\code{TRUE}) or only the variance vector
  of \eqn{\widehat{\mbox{\boldmath$\varepsilon$\unboldmath}}
  + \widehat{\mbox{\boldmath$B$\unboldmath}}
  =\mbox{\boldmath$Y$\unboldmath} - \mbox{\boldmath$X$\unboldmath}
  \widehat{\mbox{\boldmath$\beta$\unboldmath}}}{hat\epsilon+ hatB=Y-X
  hat\beta} is returned (default \code{FALSE}).}
  
  \item{aux.cov.pred.target}{logical controlling whether a covariance term
  required for the back-transformation of kriging predictions of
  log-transformed data is returned (default \code{FALSE}).}
      
  \item{min.condnum}{positive numeric.  Minimum acceptable ratio of smallest to
  largest singular value of the model matrix
  \eqn{\mbox{\boldmath$X$\unboldmath}}{X} (default \code{1.e-12}).}
  
  \item{rq}{a list of arguments passed to \code{rq} or a function such as
  \code{rq.control} that generates such a list (see
  \code{\link[quantreg]{rq}} for allowed arguments).}
  
  \item{lmrob}{a list of arguments passed to the \code{control} argument of
  \code{lmrob} or a function such as
  \code{lmrob.control} that generates such a list (see
  \code{\link[robustbase]{lmrob.control}} for allowed arguments).}
  
  \item{nleqslv}{a list of arguments passed to
  \code{{nleqslv}} or a function such as
  \code{nleqslv.control} that generates such a list (see
  \code{\link[nleqslv]{nleqslv}} for allowed arguments).}
  
  \item{optim}{a list of arguments passed to \code{optim} or a function
  such as \code{optim.control} that generates such a list (see
  \code{\link[stats]{optim}} for allowed arguments).}
  
  \item{full.output}{logical controlling how much output is returned in
  \code{georob} object.  Proper functioning requires default (\code{TRUE})}.
  
  \item{...}{named vectors of functions, extending the definition of
  transformations for variogram parameters (see \emph{Details}).}
  
  \item{variance, snugget, nugget, scale, alpha, beta, delta, gamma, kappa, lambda, mu, nu}{character
  strings with names of transformation functions of the variogram
  parameters.}
  
  \item{f1, f2, omega, phi, zeta}{character
  strings with names of transformation functions of the variogram
  parameters.}

  \item{tau, rq.method, rq.alpha, ci, iid, interp, tcrit}{arguments passed
  as \code{\ldots} to \code{\link[quantreg]{rq}}.  Note that \code{rq.} is
  stripped on passing from the argument names on passing.}
  
  \item{rq.beta, eps, Mm.factor, max.bad.fixup}{arguments passed as
  \code{\ldots} to \code{\link[quantreg]{rq}}.  Note that \code{rq.} is
  stripped on passing from the argument names on passing.}
  
  \item{nleqslv.method, global, xscalm, nleqslv.control}{arguments passed
  to related arguments of \code{\link[nleqslv]{nleqslv}}.  Note that
  \code{nleqslv.} is stripped from the argument names on passing.}
  
  \item{optim.method, lower, upper, hessian, optim.control}{arguments
  passed to related arguments of \code{\link[stats]{optim}}. Note that
  \code{optim.} is stripped from the argument names on passing.}
  
}
%   \item{}

\details{

  The arguments \code{param.tf}, \code{fwd.tf}, \code{deriv.fwd.tf},
  \code{bwd.tf} define the transformations of the variogram parameters for
  robust REML estimation.  Implemented are currently \code{"log"},
  \code{"rad"} (conversion from degree to radian) and \code{"identity"} (= no)
  transformations.  These are the possible values that the many arguments
  of the function \code{param.transf} accept (as quoted character strings)
  and these are the names of the list components returned by
  \code{fwd.transf}, \code{dfwd.transf} and \code{bwd.transf}.  Additional
  transformations can be implemented by:
  
  \enumerate{
  
    \item Extending the function definitions by arguments like
    
    \code{fwd.tf = fwd.transf(c(my.fun = function(x) your transformation))},\cr
    \code{deriv.fwd.tf = dfwd.transf(c(my.fun = function(x) your derivative))},\cr
    \code{bwd.tf = bwd.transf(c(my.fun = function(x) your back-transformation))},
    
    \item Assigning to a given argument of \code{param.transf} the name of
    the new function, e.g.\cr \code{variance = "my.fun"}.  }
  
  Note the values given for the arguments of \code{param.transf} must match
  the names of the functions returned by \code{fwd.transf},
  \code{dfwd.transf} and \code{bwd.transf}.  }

\examples{
\dontrun{
data( meuse )

r.logzn.rob <- georob(log(zinc) ~ sqrt(dist), data = meuse, locations = ~ x + y,
    variogram.model = "exponential",
    param = c( variance = 0.15, nugget = 0.05, scale = 200 ),
    tuning.psi = 1, initial.method = "rq", 
    control = georob.control(cov.bhat = TRUE, cov.ehat.p.bhat = TRUE), verbose = 2)
  
qqnorm(rstandard(r.logzn.rob, level = 0)); abline(0, 1)
qqnorm(ranef(r.logzn.rob, standard = TRUE)); abline(0, 1)
}
}

\author{
   Andreas Papritz \email{andreas.papritz@env.ethz.ch}
}

\seealso{
  \code{\link{georobIntro}} for a description of the model and a brief summary of the algorithms; 
  \code{\link{georob}} for (robust) fitting of spatial linear models;
  \code{\link{georobObject}} for a description of the class \code{georob};
  \code{\link{plot.georob}} for display of REML variogram estimates;
  \code{\link{predict.georob}} for computing robust kriging predictions; and finally
  \code{\link{georobMethods}} for further methods for the class \code{georob}.
}

% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.
% 
% }

\keyword{models}
\keyword{spatial}
\keyword{robust}
