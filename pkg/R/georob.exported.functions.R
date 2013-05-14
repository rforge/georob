georob <- 
  function( 
    formula, data, subset, weights, na.action,
    model = TRUE, x = FALSE, y = FALSE, 
    contrasts = NULL, offset, 
    locations,
    variogram.model = c( "exponential", "bessel", "cauchy", "cauchytbm",
      "circular", "cubic", "dagum", "dampedcosine", "DeWijsian", "fractalB",
      "gauss", "genB", "gencauchy", "gengneiting", "gneiting", "lgd1",
      "matern", "penta", "power", "qexponential", "spherical", "stable",
      "wave", "whittle"
    ), 
    param,
    fit.param = c( 
      variance = TRUE, snugget = FALSE, nugget = TRUE, scale = TRUE, 
      a = FALSE, alpha = FALSE, beta = FALSE, delta = FALSE, 
      gamma = FALSE, lambda = FALSE, n = FALSE, nu = FALSE
    )[ names(param) ],
    aniso = c( f1 = 1., f2 = 1., omega = 90., phi = 90., zeta = 0. ),
    fit.aniso = c( f1 = FALSE, f2 = FALSE, omega = FALSE, phi = FALSE, zeta = FALSE ),
    tuning.psi = 2, initial.param  = TRUE,
    control = georob.control( ... ), verbose = 0,
    ...
  )
{
  
  ## wrapper function to georob.fit with "standard" interface for
  ## statistical modelling
  
  ## ToDos:
  
  ## - ...$xy durch ...[["xy"]] ersetzen
  ## - Variante fuer Klasse SpatialPointsDataFrame{sp}
  ## - Ueberpruefung der Konsistenz des Inputs hier statt in georob.fit
  ## - ausgewaehlte Elemente von control in Resultate speichern
  
  
  ## History:
  
  ## georob:
  
  ## - Implementierung von Schaetzung von Startwerten der Kovarianzparameter mittels REML
  ##   von um Ausreisser bereinigten Datensatz
  ## - Implementierung von Schaetzung von Nugget und Varianz von mikroskaliger Variation (snugget )
  ## - Implementierung der Schaetzungen fuer replizierte Messungen an gleicher Messtelle    
  ## - neue Varianten um Startwerte von betahat zu rechnen (Ersatz von use.lm durch 
  ##   Argument initial.method = c( "rq", "lmrob", "rlm" ), vgl. georob.control)
  ## - Kontrolle, ob Designmatrix vollen Spaltenrang hat
  ## - Modifikation fuer Fall, dass alle Variogrammparameter fixiert sind
  ## - IRWLS Berechung von betahat und bhat entweder von Werten im initial.object 
  ##   oder von Schaetzwerten aus vorangehender Iteration
  ## - Berechung der Kovarianzen zwischen betahat und (z - bhat)
  ## - Steuerung der Berechnung der verschiedenen Kovarianzen via georob.control
  ## - vollstaendige Implementierung von Standard Interfaces fuer Input und
  ##   Output fuer statistische Modellierung (analog lm, lmrob)
  ## - korrekte Behandlung von NAs und Implementierung von subset
  ## - Verzicht auf Berechnung von initialem Variogramm
  ## - Startwerte fuer bhat alle gleich Null gesetzt
  ## - neue Struktur von initial.objects
  ## - teilweise geaenderte Namen der Argumente (model -> variogram.model)
  
  ## f.glsrob803:
  
  ## - Umbenennung einiger Argument
  
  ## f.glsrob801:
  
  ## - teilweise Implementation von Standard Interfaces fuer Input und
  ##   Output fuer statistische Modellierung 
  ## - Berechnung des Objekts mit den Startwerten
  
  
  ## 2012-04-21 AP
  ## 2012-05-03 AP bounds for safe parameter values
  ## 2012-05-07 AP correction of error for constant trend
  ## 2012-05-28 AP handle missing names of coefficients after calling update
  ## 2013-04-23 AP new names for robustness weights
  
  ## check whether input is complete
  
  if( any( missing( formula ) || missing( locations ) || missing( param ) ) )
    stop( "some mandatory arguments are missing" )
  
  ## get model frame, response vector, weights, offset and design
  ## matrix (cf.  lm, lmrob)
  
  ret.x <- x
  ret.y <- y
  
  ## vector with row number of included observations
  
  in.subset <- 1:NROW( data )
  if( !missing( subset ) ) in.subset <- in.subset[subset]
  
  ## start setting-up model frames
  
  cl <- match.call()
  mf <- match.call( expand.dots = FALSE )
  m <- match( 
    c( "formula", "data", "subset", "weights", "na.action", "offset"),
    names(mf), 0L 
  )
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name( "model.frame" )
  
  ## spot missing observations for fixed effects
  
  mf.try <- eval( mf, parent.frame() )
  if( !is.null( attr( mf.try, "na.action" ) ) ){
    ex.mf <- unname( unclass( attr( mf.try, "na.action" ) ) )
  } else ex.mf <- integer(0)
  ##     if( verbose > 3 ) cat( "ex.mf    :", ex.mf, "\n" )
  
  ## spot missing observations for locations
  
  locations <-
    as.formula( paste( deparse( locations ), "-1" ), env = parent.frame() )
  mf.loc <- mf
  mf.loc$formula <- locations
  mf.loc.try <- eval( mf.loc, parent.frame() )
  if( !is.null( attr( mf.loc.try, "na.action" ) ) ){
    ex.mf.loc <- unname( unclass( attr( mf.loc.try, "na.action" ) ) )
  } else ex.mf.loc <- integer( 0 )
  ##     if( verbose > 3 ) cat( "ex.mf.loc:", ex.mf.loc, "\n" )
  
  ## missing observations either for fixed effects or locations
  
  ex <- sort( unique( c( ex.mf.loc, ex.mf ) ) )
  ##     if( verbose > 3 ) cat( "ex       :", ex, "\n" )
  
  ## new subset argument taking account of all missing observations
  
  if( length( ex ) ) in.subset <- in.subset[-ex]
  ##     if( verbose > 3 ) cat( "in.subset:", in.subset, "\n" )
  mf$subset <- mf.loc$subset <- in.subset
  
  ## final model frames for fixed effects and locations
  
  mf     <- eval( mf,     parent.frame() )
  mf.loc <- eval( mf.loc, parent.frame() )
  
  mt     <- attr( mf,     "terms" )
  mt.loc <- attr( mf.loc, "terms" )
  
  ## check whether 'empty' models have been entered
  
  if( is.empty.model( mt ) )
    stop( "an 'empty' fixed effects model has been specified" )
  if( is.empty.model( mt.loc ) )
    stop( "an 'empty' locations model has been specified" )
  
  ## check whether fixed effects model includes an intercept if an
  ## intrinsic variogram model is used
  
  if( identical( attr( mt, "intercept" ), 0L ) && 
    variogram.model %in% georob.control()$irf.model )
    stop(
      "the fixed effects model must include an intercept ",
      "if an unbounded variogram model is used"
    )
  
  ## extract fixed effects response variable, weights, offset and design matrix
  
  y <- model.response( mf, "numeric" )
  
  w <- as.vector( model.weights( mf ) )
  if( !is.null(w) )
    stop( "weights are not yet implemented for this estimator" )
  
  offset <- as.vector( model.offset(mf) )
  if( !is.null(offset) ) {
    if( length( offset ) != NROW(y) )
      stop( gettextf(
        "number of offsets is %d, should equal %d (number of observations)", 
        length(offset), NROW(y) ), domain = NA )
  }
  
  ## check if optionally provided bhat has correct length
  
  if( !is.null( control$bhat ) && length( y ) != length( control$bhat ) ) stop(
    "lengths of response vector and 'bhat' do not match"    
  )
  
  x <- model.matrix( mt, mf, contrasts )
  
  ## check whether design matrix has full column rank
  
  min.max.sv <- range( svd( crossprod( x ) )$d )
  condnum <- min.max.sv[1] / min.max.sv[2] 
  
  
  if( condnum <= control$min.condnum ) stop( 
    "design matrix has not full column rank (condition number of X^T X: ", 
    signif( condnum, 2 ), ")"
  )
  
  ## subtract offset
  
  yy <- y
  if( !is.null( offset ) ) yy <- yy - offset
  
  ## adjust choice for initial.method to compute regression coefficients
  
  if( initial.param ) control$initial.method <- "lmrob"
  
  ## compute initial guess of fixed effects parameters (betahat)
  
  r.initial.fit <- switch(
    control$initial.method,
    rq = {
      
      Rho <- function( u, tau) u * (tau - (u < 0))
      tau <- control$rq$tau
      process <- (tau < 0 || tau > 1)
      
      f.rq.fit <- rq.fit
      formals( f.rq.fit ) <- c( alist( x=, y= ), control$rq, alist( ...= ) )
      
      fit <- f.rq.fit( x = x, y = yy, ... ) 
      
      if( process ) {
        rho <- list(x = fit$sol[1, ], y = fit$sol[3, ])
      } else {
        dimnames(fit$residuals) <- list( dimnames( x )[[1]], NULL )
        rho <- sum( Rho( fit$residuals, tau ) )
      }
      if( control$rq$method == "lasso" ){
        class(fit) <- c("lassorq", "rq")
      } else if( control$rq$method == "scad"){
        class(fit) <- c("scadrq", "rq")
      } else {
        class(fit) <- ifelse(process, "rq.process", "rq")
      }
      fit$na.action <- attr( mf, "na.action" )
      fit$formula <- formula
      fit$terms <- mt
      fit$xlevels <- .getXlevels(mt, mf)
      fit$call <- call
      fit$tau <- tau
      fit$weights <- w
      fit$residuals <- drop( fit$residuals )
      fit$rho <- rho
      fit$method <- control$rq$method
      fit$fitted.values <- drop( fit$fitted.values )
      attr(fit, "na.message") <- attr( m, "na.message" )
      if( model ) fit$model <- mf
      fit
      
    },
    lmrob = {
      
      fit <- lmrob.fit( x, yy, control = control$lmrob )
      fit$na.action <- attr(mf, "na.action")
      fit$offset <- offset
      fit$contrasts <- attr(x, "contrasts")
      fit$xlevels <- .getXlevels(mt, mf)
      fit$call <- cl
      fit$terms <- mt
      if( control$lmrob$compute.rd && !is.null(x) )
      fit$MD <- robustbase:::robMD( x, attr( mt, "intercept" ) )
      if( !is.null( offset ) ) fit$fitted.values + offset
      fit
      
    }
  )
  
  ## match variogram model
  
  variogram.model <- match.arg( variogram.model )  
  
  ## compute coordinates of locations and distance object
  
  locations.coords <- model.matrix( mt.loc, mf.loc )
  
  if( 
    !( missing( aniso ) || missing( fit.aniso ) ) && 
    ( NCOL( locations.coords ) < 2 || NCOL( locations.coords ) > 3 )
  ) stop( 
    "anisotropic variogram models are implemented only for 2 or 3 dimensions" 
  )

  names( yy ) <- rownames( mf )
  
  ##  create environment to store items required to compute likelihood and
  ##  estimating equations that are provided by
  ##  prepare.likelihood.calculations
  
  envir <- new.env()
  lik.item <- list()
  assign( "lik.item", lik.item, pos = as.environment( envir ) )
  
  ## check whether anisotropy parameters were passed to georob
  
  aniso.missing <- missing( aniso ) && missing( fit.aniso )
  
  ## prune data set for computing initial values of variogram and
  ## anisotropy parameters by reml
  
  if( initial.param && tuning.psi < control$tuning.psi.nr ){
    
    t.sel <- switch(
      control$initial.method,
      lmrob = r.initial.fit$rweights > control$min.rweight,
      stop(
        "computing initial values of covariance parameters requires 'lmrob' as ",
        "method for computing initial regression coefficients"
      )
    )
    
    if( verbose > 0 ) cat( 
      "\ndiscarding", sum( !t.sel ), "of", length( t.sel ), 
      "observations for computing initial estimates of variogram\nand anisotropy parameter by gaussian reml\n"
    )
    
    ## collect.items for initial object
    
    initial.objects <- list(
      x = as.matrix( x[t.sel, ] ),
      y = yy[t.sel],
      betahat = coef( r.initial.fit ),
      bhat = if( is.null( control$bhat ) ){
        rep( 0., length( yy ) )[t.sel]
      } else {
        control$bhat[t.sel]
      },
      initial.fit = r.initial.fit,
      locations.objects = list(
        locations = locations,
        coordinates = locations.coords[t.sel, ]
      ),
      isotropic = aniso.missing
    )
    
    ## estimate model parameters with pruned data set
    
    t.georob <- georob.fit(
      envir = envir,
      initial.objects = initial.objects,
      variogram.model = variogram.model, param = param, fit.param = fit.param,
      aniso = aniso, fit.aniso = fit.aniso,
      param.tf = control$param.tf,
      fwd.tf = control$fwd.tf, 
      deriv.fwd.tf = control$deriv.fwd.tf,
      bwd.tf = control$bwd.tf, 
      safe.param = control$safe.param,
      tuning.psi = control$tuning.psi.nr,
      cov.bhat = control$cov.bhat, full.cov.bhat = control$full.cov.bhat,
      cov.betahat = control$cov.betahat, 
      cov.bhat.betahat = control$cov.bhat.betahat,
      cov.delta.bhat = control$cov.delta.bhat,
      full.cov.delta.bhat = control$full.cov.delta.bhat,
      cov.delta.bhat.betahat = control$cov.delta.bhat.betahat,
      cov.ehat = control$cov.ehat, full.cov.ehat = control$full.cov.ehat,
      cov.ehat.p.bhat = control$cov.ehat.p.bhat, full.cov.ehat.p.bhat = control$full.cov.ehat.p.bhat,
      aux.cov.pred.target = control$aux.cov.pred.target,
      psi.func = control$psi.func,
      tuning.psi.nr = control$tuning.psi.nr,
      irwls.initial = control$irwls.initial,
      irwls.maxiter = control$irwls.maxiter,
      irwls.reltol = control$irwls.reltol,
      force.gradient = control$force.gradient,
      zero.dist = control$zero.dist,
      nleqslv.method =  control$nleqslv$method,
      nleqslv.control = control$nleqslv$control,
      optim.method =  control$optim$method,
      optim.lower = control$optim$lower,
      optim.upper = control$optim$upper,
      hessian =       control$optim$hessian,
      optim.control = control$optim$control,
      full.output = control$full.output,
      verbose = verbose
    )
    
    param = t.georob$param[names(fit.param)]
    aniso = t.georob$aniso$aniso[names(fit.aniso)] * c( 1, 1, rep( 180/pi, 3 ) )
    
  }
  
  ## collect.items for initial object
  
  initial.objects <- list(
    x = as.matrix( x ),
    y = yy,
    betahat = coef( r.initial.fit ),
    bhat = if( is.null( control$bhat ) ){
      rep( 0., length( yy ) )
    } else {
      control$bhat
    },
    initial.fit = r.initial.fit,
    locations.objects = list(
      locations = locations,
      coordinates = locations.coords
    ),
    isotropic = aniso.missing
  )
  
  ## estimate model parameters
  
  r.georob <- georob.fit(
    envir = envir,
    initial.objects = initial.objects,
    variogram.model = variogram.model, param = param, fit.param = fit.param,
    aniso = aniso, fit.aniso = fit.aniso,
    param.tf = control$param.tf,
    fwd.tf = control$fwd.tf, 
    deriv.fwd.tf = control$deriv.fwd.tf,
    bwd.tf = control$bwd.tf, 
    safe.param = control$safe.param,
    tuning.psi = tuning.psi,
    cov.bhat = control$cov.bhat, full.cov.bhat = control$full.cov.bhat,
    cov.betahat = control$cov.betahat, 
    cov.bhat.betahat = control$cov.bhat.betahat,
    cov.delta.bhat = control$cov.delta.bhat,
    full.cov.delta.bhat = control$full.cov.delta.bhat,
    cov.delta.bhat.betahat = control$cov.delta.bhat.betahat,
    cov.ehat = control$cov.ehat, full.cov.ehat = control$full.cov.ehat,
    cov.ehat.p.bhat = control$cov.ehat.p.bhat, full.cov.ehat.p.bhat = control$full.cov.ehat.p.bhat,
    aux.cov.pred.target = control$aux.cov.pred.target,
    psi.func = control$psi.func,
    tuning.psi.nr = control$tuning.psi.nr,
    irwls.initial = control$irwls.initial,
    irwls.maxiter = control$irwls.maxiter,
    irwls.reltol = control$irwls.reltol,
    force.gradient = control$force.gradient,
    zero.dist = control$zero.dist,
    nleqslv.method =  control$nleqslv$method,
    nleqslv.control = control$nleqslv$control,
    optim.method =  control$optim$method,
    optim.lower = control$optim$lower,
    optim.upper = control$optim$upper,
    hessian =       control$optim$hessian,
    optim.control = control$optim$control,
    full.output = control$full.output,
    verbose = verbose
  )
  
  ## add offset to fitted values
  
  if( !is.null( offset ) )
    r.georob$fitted.values <- r.georob$fitted.values + offset
  
  ## add remaining items to output
  
  if( control$full.output ){
    
    r.georob$initial.objects$initial.param <- initial.param
    
    if( control$lmrob$compute.rd && !is.null( x ) )
       r.georob$MD <- robustbase:::robMD( x, attr(mt, "intercept") )
    if( model ) r.georob$model <- mf
    if( ret.x ) r.georob$x <- x
    if( ret.y ) r.georob$y <- y
    r.georob$df.residual <- -diff( dim( initial.objects$x ) )
    
    r.georob$na.action <- attr(mf, "na.action")
    r.georob$offset <- offset
    r.georob$contrasts <- attr(x, "contrasts")
    r.georob$xlevels <- .getXlevels(mt, mf)
    r.georob$rank <- ncol( initial.objects$x )
    r.georob$call <- cl
    r.georob$terms <- mt
    
  }
  
  ## set missing names of coefficients (bug of update)
  
  if( length( r.georob$coefficients ) == 1 && is.null( names( r.georob$coefficients ) ) ){
    names( r.georob$coefficients ) <- "(Intercept)"
  }

  
  class( r.georob ) <- c( "georob" )
  
  invisible( r.georob )
  
}

#  ##############################################################################

georob.control <- 
  function(
    initial.method = c("lmrob", "rq"),
    bhat = NULL,
    param.tf = param.transf(),
    fwd.tf = fwd.transf(), 
    deriv.fwd.tf = dfwd.transf(), 
    bwd.tf = bwd.transf(),
    safe.param = 1.e12,
    psi.func = c( "logistic", "t.dist", "huber" ),
    tuning.psi.nr = 1000,
    min.rweight = 0.25,
    irwls.initial = TRUE,
    irwls.maxiter = 50, irwls.reltol = sqrt( .Machine$double.eps ),
    force.gradient = FALSE,
    zero.dist = sqrt( .Machine$double.eps ),
    cov.bhat = FALSE, full.cov.bhat = FALSE,
    cov.betahat = TRUE, 
    cov.bhat.betahat = FALSE,
    cov.delta.bhat = TRUE, full.cov.delta.bhat = TRUE,
    cov.delta.bhat.betahat = TRUE,
    cov.ehat = TRUE, full.cov.ehat = FALSE,
    cov.ehat.p.bhat = FALSE, full.cov.ehat.p.bhat = FALSE,
    aux.cov.pred.target = FALSE,
    min.condnum = 1.e-12,
    rq = rq.control(),
    lmrob = lmrob.control(),
    nleqslv = nleqslv.control(),
    optim = optim.control(),
    full.output = TRUE
  )
{
  
  ## auxiliary function to set meaningful default values for the
  ## arguments of the function georob.fit
  
  ## Arguments:
  
  ## initial.method    character scalar, controlling how the intitial estimate of the fixed-effects
  ##                   parameters are computed, possible values are 
  ##                   "rq"    to use rq{quantreg},
  ##                   "lmrob" to use lmrob{robustbase},
  ## param.tf       list, used to pass arguents to param.tf{georob}
  ##                   parameters, implemented values are "log" or "identity" (no transformation)
  ## fwd.tf
  ## rho.function      character, defining the rho/psi functions family
  ## tuning.psi.nr numeric, if tuning.psi exceeds tuning.psi.nr for
  ##                   logistic or huber rho.function then only one IRWLS iteration is executed
  ##                   to estimate beta and z
  ## min.rweight  minimum robustness weights of lmrob fit required for 
  ##                   including an observations into the pruned data set from which initial values
  ##                   of variogram and anisotropy parameters are computed by Gaussian REML
  ## irwls.initial  logical, flag controlling whether IRWLS starts from the lmrob 
  ##                   estimates of beta and from z=0 (TRUE) or from the previous IRWLS results
  ## irwls.maxiter     integer, maximum number of IRWLS steps
  ## irwls.reltol      numeric, relative convergence tolerance for IRWLS, see optim{stats}
  ## force.gradient    logical, flag controlling whether the gradient (REML) or the
  ##                   estimation equations should be evaluated if all variogram parameter are 
  ##                   fixed
  ## zero.dist         observations from sampling locations less than zero.dist apart will be 
  ##				   considered as multiple observations from same location
  ## cov.bhat        logical, flag controlling whether the covariances of bhat should be computed
  ## full.cov.bhat   logical, flag controlling whether the full covariance matrix of bhat 
  ##                   is computed (TRUE) or only the diagonal elements (FALSE)
  ## cov.betahat     logical, flag controlling whether the covariance matrix of betahat 
  ##                   should be computed
  ## cov.bhat.betahat  logical, flag controlling whether the covariance matrix of 
  ##                   bhat and betahat should be computed
  ## cov.delta.bhat      logical, flag controlling whether the covariances of z-bhat should be computed
  ## full.cov.delta.bhat logical, flag controlling whether the full covariance matrix of z-bhat 
  ##                   is computed (TRUE) or only the diagonal elements (FALSE)
  ## cov.delta.bhat.betahat    logical, flag controlling whether the covariance matrix of z-bhat
  ##                    and betahat should be computed
  ## cov.ehat        logical, flag controlling whether the covariances of the resdiuals should be computed
  ## full.cov.ehat   logical, flag controlling whether the full covariance matrix of the residuals 
  ##                   is computed (TRUE) or only the diagonal elements (FALSE)
  ## cov.ehat.p.bhat       logical, flag controlling whether the covariances of the resdiuals+bhat should be computed
  ## full.cov.ehat.p.bhat  logical, flag controlling whether the full covariance matrix of the resdiuals+bhat 
  ##                   is computed (TRUE) or only the diagonal elements (FALSE)
  ## aux.cov.pred.target  logical, flag controlling whether the auxiliary matrix for computing the covariances
  ##                   of the predicted and true y should be computed
  ##                   is computed (TRUE) or only the diagonal elements (FALSE)
  ## min.condnum       minimum condition number for a matrix to be numerically non-singular
  ## rq                list, see rq{quantreg}
  ## lmrob             list, see lmrob.control{robustbase}
  ## nleqslv           list, used to pass arguent to nleqslv, see nleqslv.control{georob}
  ## optim             list, used to pass arguments to optim, see optim.control{georob}
  ## full.output       logical, flag used to control the amount of output returned by georob, warning:
  ##                   is TRUE then the output will not contain all required items required by some methods
  
  
  ## 2012-04-21 A. Papritz   
  ## 2012-05-03 AP bounds for safe parameter values
  ## 2012-05-04 AP modifications for lognormal block kriging
  ## 2013-04-23 AP new names for robustness weights
  
  if( 
    !( all( param.tf %in% names( fwd.tf ) ) &&
       all( param.tf %in% names( deriv.fwd.tf ) ) &&
       all( param.tf %in% names( bwd.tf ) )
    )
  ) stop( 
    "undefined transformation of variogram parameters; extend respective function definitions" 
  )
  
  list(
    initial.method = match.arg( initial.method ),
    bhat = bhat,
    param.tf = param.tf, fwd.tf = fwd.tf, deriv.fwd.tf = deriv.fwd.tf, bwd.tf = bwd.tf,
    safe.param = safe.param,
    psi.func = match.arg( psi.func ), 
    tuning.psi.nr = tuning.psi.nr,
    min.rweight = min.rweight,
    irwls.initial = irwls.initial,
    irwls.maxiter = irwls.maxiter, irwls.reltol = irwls.reltol,
    force.gradient = force.gradient,
    zero.dist = zero.dist,
    cov.bhat = cov.bhat, full.cov.bhat = full.cov.bhat,
    cov.betahat = cov.betahat, 
    cov.bhat.betahat = cov.bhat.betahat,
    cov.delta.bhat = cov.delta.bhat, full.cov.delta.bhat = full.cov.delta.bhat,
    cov.delta.bhat.betahat = cov.delta.bhat.betahat,
    cov.ehat = cov.ehat, full.cov.ehat = full.cov.ehat,
    cov.ehat.p.bhat = cov.ehat.p.bhat, full.cov.ehat.p.bhat = full.cov.ehat.p.bhat,
    aux.cov.pred.target = aux.cov.pred.target,
    min.condnum = min.condnum,
    irf.models = c( "DeWijsian", "fractalB", "genB" ),
    rq = rq, lmrob = lmrob, nleqslv = nleqslv, optim = optim, 
    full.output = full.output,
    stored.items = list( 
      Valpha.objects = c( "Valpha.inverse", "Valpha.ilcf", "Valpha.ucf", "gcr.constant" )
    )
  )
  
}

## ======================================================================
param.transf <-
  function( 
    variance = "log", snugget = "log", nugget = "log", scale = "log", 
    a = "identity", alpha = "identity", beta = "identity", delta = "identity", 
    gamma = "identity", lambda = "identity", n = "identity", nu = "identity",
    f1 = "log", f2  ="log", omega = "identity", phi = "identity", zeta = "identity"
  )
{
  
  ## function sets meaningful defaults for transformation of variogram
  ## parameters
  
  ## 2012-11-27 A. Papritz
  
  c( 
    variance = variance, snugget = snugget, nugget = nugget, scale = scale,
    a = a, alpha = alpha, beta = beta, delta = delta, gamma = gamma, 
    lambda = lambda, n = n, nu = nu, 
    f1 = f1, f2 = f2, omega = omega, phi = phi, zeta = zeta
  )
  
}

## ======================================================================
fwd.transf <- 
  function( 
    ...
  )
{
  
  ## definition of forward transformation of variogram parameters
  
  ## 2012-11-27 A. Papritz
  
  list( log = function(x) log(x),  identity = function(x) x, ... )
}

## ======================================================================
dfwd.transf<-
  function( 
    ...
  )
{
  
  ## definition of first derivative of forward transformation of variogram
  ## parameters
  
  ## 2012-11-27 A. Papritz
  
  list( log = function(x) 1/x, identity = function(x) rep(1, length(x)), ... )  
  
}

## ======================================================================
bwd.transf <-
  function( 
    ...
  )
{
  
  ## definition of backward transformation of variogram parameters
  
  ## 2012-11-27 A. Papritz
  
  list( log = function(x) exp(x), identity = function(x) x, ... )
  
}


## ======================================================================
rq.control <-
  function(
    ## arguments for rq
    tau = 0.5, rq.method = "br",
    ## arguments for rq.fit.br
    rq.alpha = 0.1, ci = FALSE, iid = TRUE, interp = TRUE, tcrit = TRUE,
    ## arguments for rq.fit.fnb
    rq.beta = 0.99995, eps = 1e-06,
    ## arguments for rq.fit.pfn
    Mm.factor = 0.8, max.bad.fixup = 3
  )
{
  
  ## function sets meaningful defaults for selected arguments of function
  ## rq{quantreg}
  
  ## 2012-12-14 A. Papritz
  
  list( 
    tau = tau,  method = rq.method,
    alpha = rq.alpha, ci = ci, iid = iid, interp = interp, tcrit = tcrit,
    beta = rq.beta, eps = eps,
    Mm.factor = Mm.factor, max.bad.fixup = max.bad.fixup
  )
}


## ======================================================================
nleqslv.control <-
  function( 
    nleqslv.method = c( "Broyden", "Newton"), 
    global = c( "dbldog", "pwldog", "qline", "gline", "none" ),
    xscalm = c( "fixed", "auto" ),
    nleqslv.control = NULL
  )
{
  
  ## function sets meaningful defaults for selected arguments of function
  ## nleqslv{nleqslv} 
  
  ## 2012-12-14 A. Papritz
  
  aux <- function( trace = 0, ... ) list( trace = trace, ... )
  
  list( 
    method = match.arg( nleqslv.method ),
    global = match.arg( global ),
    xscalm = match.arg( xscalm ),
    control = nleqslv.control
  )
}

## ======================================================================
optim.control <-
  function( 
    optim.method = c( "BFGS", "Nelder-Mead", "CG", "L-BFGS-B", "SANN", "Brent" ),
    lower = -Inf, upper = Inf,
    optim.control = NULL,
    hessian = TRUE
  )
{
  
  ## function sets meaningful defaults for selected arguments of function optim
  ## 2012-12-14 A. Papritz
  
  aux <- function( trace = 0, maxit = 500, ... ) list( trace = trace, maxit = maxit, ... )
  
  list( 
    method = match.arg( optim.method ),
    lower = lower, upper = upper,
    hessian = hessian,
    control = optim.control
  )
}


## ======================================================================

compress <- 
  function( m )
{
  
  ## function stores a list of or a single lower, upper triangular or
  ## symmetric matrix compactly
  
  aux <- function( x ){
    struc <- attr( x, "struc" )
    if( !is.null( struc ) ){
      switch( 
        struc,
        sym = {
          aux <- list( diag = diag( x ), tri = x[lower.tri(x)] )
          attr( aux, "struc" ) <- "sym"
          aux
        },
        lt = {
          aux <- list( diag = diag( x ), tri = x[lower.tri(x)] )
          attr( aux, "struc" ) <- "lt"
          aux
        }
        ,
        ut = {
          aux <- list( diag = diag( x ), tri = x[upper.tri(x)] )
          attr( aux, "struc" ) <- "ut"
          aux
        }
      )
    } else {
      x
    }
  }
  
  if( is.list( m ) ){
    lapply( m, aux )
  } else {
    aux ( m )
  }
  
  
  
}

## ======================================================================

expand <- 
  function( object )
{
  
  ## function expands a list of or a compactly stored lower, upper
  ## triangular or symmetric matrices
  
  aux <- function( x ){
    struc <- attr( x, "struc" )
    if( !is.null( struc ) ){
      switch( 
        struc,
        sym = {
          n <- length( x$diag )
          dn <- names( x$diag )
          aux <- matrix( 0., n, n )
          aux[lower.tri( aux )] <- x$tri
          aux <- aux + t( aux )
          diag( aux ) <- x$diag
          dimnames( aux ) <- list( dn, dn )
          attr( aux, "struc" ) <- "sym"
          aux
        },
        lt = {
          n <- length( x$diag )
          dn <- names( x$diag )
          aux <- matrix( 0., n, n )
          aux[lower.tri( aux )] <- x$tri
          diag( aux ) <- x$diag
          dimnames( aux ) <- list( dn, dn )
          attr( aux, "struc" ) <- "lt"
          aux
        }
        ,
        ut = {
          n <- length( x$diag )
          dn <- names( x$diag )
          aux <- matrix( 0., n, n )
          aux[upper.tri( aux )] <- x$tri
          diag( aux ) <- x$diag
          dimnames( aux ) <- list( dn, dn )
          attr( aux, "struc" ) <- "ut"
          aux
        }
      )
    } else {
      x
    }
  }
  
  ln <- names( object )
  if( is.list( object ) ){
    if( length( ln ) == 2 && all( ln == c( "diag", "tri" ) ) ){
      aux( object )
    } else {
      lapply( object, aux )
    }
  } else {
    object
  }  
  
}

## ======================================================================

param.names <- 
  function( model )
{
  
  ## function returns names of extra parameters of implemented variogram
  ## models (cf. Variogram{RandomFields})
  
  ## 2012-01-24 A. Papritz
  
  switch(
    model,
    "bessel"        = "nu",
    "cauchy"        = "beta",
    "cauchytbm"     = c( "alpha", "beta", "gamma" ),
    "circular"      = NULL,
    "cubic"         = NULL,
    "dagum"         = c( "beta", "gamma" ),
    "dampedcosine"  = "lambda",
    "DeWijsian"     = "alpha",
    "exponential"   = NULL,
    "fractalB"      = "alpha",
    "gauss"         = NULL,
    "genB"          = c( "alpha", "delta" ),
    "gencauchy"     = c( "alpha", "beta" ),
    "gengneiting"   = c( "n", "alpha" ),
    "gneiting"      = NULL,
    "lgd1"          =  c( "alpha", "beta" ),
    "matern"        = "nu",
    "penta"         = NULL,
    "power"         = "a",
    "qexponential"  = "alpha",
    "spherical"     = NULL,
    "stable"        = "alpha",
    "wave"          = NULL,
    "whittle"       = "nu",
    stop( model, " variogram not implemented" )
  )
}

##  ##############################################################################

param.bounds <- 
  function( model, d, param )
{
  
  ## function returns range of parameters for which variogram models are
  ## valid (cf.  Variogram{RandomFields})
  
  ## 2012-03-30 A. Papritz
  
  switch(
    model,
    "bessel"        = list( nu = c( 0.5 * (d - 2.), Inf ) ),
    "cauchy"        = list( beta = c( 1.e-18, Inf ) ),
    "cauchytbm"     = list( alpha = c( 1.e-18, 2.), beta = c( 1.e-18, Inf), gamma = c( d, Inf) ),
    "circular"      = NULL,
    "cubic"         = NULL,
    "dagum"         = list( beta = c( 1.e-18, 1.), gamma = c( 1.e-18, 1.-1.e-18) ),
    "dampedcosine"  = list( lambda = c( if( d > 2 ) sqrt(3.) else 1., Inf ) ),
    "DeWijsian"     = list( alpha = c( 1.e-18, 2 ) ),
    "exponential"   = NULL,
    "fractalB"      = list( alpha = c( 1.e-18, 2.) ),
    "gauss"         = NULL,
    "genB"          = list( alpha = c(1.e-18, 2.), delta = c(1.e-18, 1.-1.e-18) ),
    "gencauchy"     = list( alpha = c(1.e-18, 2.), beta = c(1.e-18, Inf) ),
    "gengneiting"   = list( n = c(1, 3), alpha = c( 0.5 * (d + 2 * unname(param["n"]) + 1), Inf ) ),
    "gneiting"      = NULL,
    "lgd1"          = list( 
                       alpha = c( 
                         1.e-18, 
                         if( d <= 2 ) 0.5 * (3-d) else stop("dimension > 3 not allowed for lgd1 model" ) 
                       ), 
                       beta = c(1.e-18, Inf)
                     ),
    "matern"        = list( nu = c(1.e-18, Inf) ),
    "penta"         = NULL,
    "power"         = list( a = c( 0.5 * (d + 1), Inf ) ),
    "qexponential"  = list( alpha = c(0., 1.) ),
    "spherical"     = NULL,
    "stable"        = list( alpha = c(1.e-18, 2.) ),
    "wave"          = NULL,
    "whittle"       = list( nu = c(1.e-18, Inf) ),
    stop( model, " variogram not implemented" )
  )
}
