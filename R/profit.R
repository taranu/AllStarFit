.repl <- function(val, torepl)
{
	return(rep(val, length(torepl)))
}

timeinmins <- function() {
  return(proc.time()["elapsed"]/60)
}

exceededmaxtime <- function(tinit, tfinal, iter=1e10) {
	t = timeinmins()
	if(t > tfinal) return(TRUE)
  return((t - tinit)*(iter+1)/iter > (tfinal-tinit))
}

profitGetExtraDefaults <- function(profile)
{
	extraparams = c()
	if(profile == "brokenexp")
	{
		extraparams = c(a=0,rb=0)
	}
	else if(profile == "coresersic")
	{
		extraparams = c(a=0,b=0,rb=0)
	}
	else if(profile == "ferrer" || profile == "ferrers")
	{
		extraparams = c(b=0)
	}
	else if((profile == "moffat") || (profile == "pointsource") || (profile == "sersic"))
	{
	}
	else
	{
		stop(paste0("Unknown profile=['",profile,"']; aborting."))
	}
	return(extraparams)
}

profitGetExtraIntervalDefaults <- function(profile)
{
	extraparams = c()
	if(profile == "brokenexp")
	{
		extraparams = c(a=list(lim=c(0,1e10)),rb=list(lim=c(0,1e10)))
	}
	else if(profile == "coresersic")
	{
		extraparams = c(a=0,b=0,rb=0)
	}
	else if(profile == "ferrer" || profile == "ferrers")
	{
		extraparams = c(b=0)
	}
	else if((profile == "moffat") || (profile == "pointsource") || (profile == "sersic"))
	{
	}
	else
	{
		stop(paste0("Unknown profile=['",profile,"']; aborting."))
	}
	return(extraparams)
}

profitRepExtra <- function(profile, values=profitGetExtraDefaults(profile), nrep=1)
{
	extraparams = list()
	if(profile == "brokenexp")
	{
		extraparams = c("a","rb")
	}
	else if(profile == "coresersic")
	{
		extraparams = c("a","b","rb")
	}
	else if(profile == "ferrer" || profile == "ferrers")
	{
		extraparams = "b"
	}
	else if((profile == "moffat") || (profile == "pointsource") || (profile == "sersic"))
	{
	}
	else
	{
		stop(paste0("Unknown profile=['",profile,"']; aborting."))
	}
	rv = list()
	if(length(extraparams) > 0)
	{
		for(param in extraparams) rv[[param]] = rep(values[[param]],nrep)
	}
	return(rv)
}

profitGetProfileParamNames <- function(profile)
{
	extraparams = list()
	if(profile == "brokenexp")
	{
		sizename = "h1"
		slopename = "h2"
		extraparams = c("a","rb")
	}
	else if(profile == "coresersic")
	{
		sizename = "re"
		slopename = "nser"
		extraparams = c("a","b","rb")
	}
	else if(profile == "ferrer" || profile == "ferrers")
	{
		sizename = "rout"
		slopename = "a"
		extraparams = "b"
	}
	else if(profile == "moffat")
	{
		sizename = "fwhm"
		slopename = "con"
	}
	else if(profile == "pointsource")
	{
		sizename = NULL
		slopename = NULL
	}
	else if(profile == "sersic")
	{
		sizename = "re"
		slopename = "nser"
	}
	else
	{
		stop(paste0("Unknown profile=['",profile,"']; aborting."))
	}
	rv = list(size=sizename,slope=slopename,profile=profile,
		extraparams=extraparams)
	return(rv)
}

profitGetProfileSlopeLimits <- function(profile)
{
	lims = c()
	errdesc = "Unknown"
	if((profile == "brokenexp") || (profile == "coresersic") || (profile == "sersic"))
	{
		lims = c(0,Inf)
	}
	else if(profile == "moffat")
	{
		lims = c(1,Inf)
	}
	else if(profile == "pointsource")
	{
		lims = c(-Inf,Inf)
	}
	else if(profile == "ferrer")
	{
		errdesc="Unsupported"
	}
	if(length(lims) != 2) stop(paste0(errdesc," profile=['",profile,"']; aborting."))
	return(lims)
}

profitMakeList <- function(xcen=0, ycen=0, mag=0, size=0, slope=0,
	ang=0, axrat=1, box=0, profile="sersic", ispsf=FALSE,
	extraparams=profitGetExtraDefaults(profile))
{
	if(profile == "pointsource") return(list(
		xcen=xcen, ycen=ycen,mag=mag
	))

	varnames = profitGetProfileParamNames(profile)
	sizename = varnames$size
	slopename = varnames$slope
	stopifnot(identical(length(varnames$extraparams), length(extraparams)))

	plist = list(mag=mag)
	if(!ispsf) plist = c(list(xcen=xcen,ycen=ycen),plist)
	plist = c(plist, list(
    size = size,
    slope = slope,
    ang = ang,
    axrat = axrat,
    box = box)
  )
	nml = names(plist)
	names(plist)[match(c("size","slope"),nml)] = c(sizename,slopename)
	for(name in names(extraparams)) plist[[name]] = extraparams[[name]]
	return(plist)
}

profitMakeLists <- function(xcen=0, ycen=0, mag=0, size=0, slope=0,
	ang=0, axrat=1, box=0, skybg=0, profile="sersic", extra=profitRepExtra(profile,nrep=length(xcen)),
	ispsf=FALSE, addsky=FALSE,
	fitx=.repl(TRUE,xcen), fity=fitx, fitmag=.repl(TRUE,mag),
	fitsize=.repl(TRUE,size), fitslope=.repl(TRUE,slope),
	fitang=.repl(TRUE,ang), fitaxrat=.repl(TRUE,axrat), fitbox=.repl(FALSE,box), fitsky=FALSE,
	fitextra=profitRepExtra(profile,values=FALSE & profitGetExtraDefaults("brokenexp"),nrep=length(xcen)),
	logsize=.repl(TRUE,size), logslope=.repl(TRUE,slope),
	logang=.repl(FALSE,ang), logaxrat=.repl(TRUE,axrat),
	logextra=profitRepExtra(profile,values=FALSE & profitGetExtraDefaults("brokenexp"),nrep=length(xcen)),
	xmin=-Inf, xmax=Inf, ymin=xmin, ymax=xmax, magmin=-Inf, magmax=Inf,
	sizemin=0, sizemax=Inf, slopemin=profitGetProfileSlopeLimits(profile)[1],
	slopemax=profitGetProfileSlopeLimits(profile)[2], angmin=-180, angmax=360,
	axratmin=0, axratmax=1, skybgmin=-Inf, skybgmax=Inf, boxmin=-1, boxmax=1,
	extramin=profitRepExtra(profile,nrep=length(xcen)), extramax=profitRepExtra(profile,nrep=length(xcen)),
	xmean=xcen, ymean=ycen, magmean=mag, sizemean=size,
	slopemean=slope, angmean=ang, axratmean=axrat, boxmean=box, skybgmean=skybg,
	extramean=profitRepExtra(profile,nrep=length(xcen)),
	xsd=Inf, ysd=Inf, magsd=Inf, sizesd=Inf, slopesd=Inf,
	angsd=Inf, axratsd=Inf, boxsd=Inf, skybgsd=Inf,
	extrasd=profitRepExtra(profile,Inf+profitGetExtraDefaults(profile),nrep=length(xcen)),
	priortofitonly=FALSE, unlog=FALSE
)
{
	isps = profile=="pointsource"
	profilenames = profitGetProfileParamNames(profile)
	extranames = profilenames$extraparams
	nextranames = length(extranames)
	tocheck = list(extra=extra, fitextra=fitextra, logextra=logextra, extramin=extramin,
		extramax = extramax, extramean=extramean, extrasd=extrasd)
	nprofiles = length(xcen)
	for(check in tocheck)
	{
		stopifnot(identical(length(check),nextranames))
		for(name in extranames)
		{
			stopifnot(name %in% names(check))
			stopifnot(identical(length(check[[name]]),nprofiles))
		}
	}
	ncomp = length(xcen)
	tocheck = list(mag=mag)
	if(!ispsf) tocheck = c(list(xcen=xcen,ycen=ycen),tocheck)
	if(!isps) tocheck = c(tocheck,list(mag=mag,size=size,slope=slope,ang=ang,axrat=axrat,box=box))
	stopifnot(all(ncomp==unlist(lapply(tocheck,length))))
	tocheck = list(mag)
	if(!ispsf) tocheck = c(list(fitx,fity,fitmag),tocheck)
	if(!isps) tocheck = c(tocheck,list(fitsize,fitslope,fitang,fitaxrat,fitbox))
	stopifnot(all(ncomp==unlist(lapply(tocheck,length))))
	modellist = list(profitMakeList(xcen=xcen,ycen=ycen,mag=mag,size=size,
		slope=slope,ang=ang,axrat=axrat,box=box,profile=profile,ispsf=ispsf,extraparams=extra))
	names(modellist) = profile
	if(addsky) modellist$sky = list(bg=skybg)
	tofit = list(profitMakeList(xcen=fitx,ycen=fity,mag=fitmag,size=fitsize,
		slope=fitslope,ang=fitang,axrat=fitaxrat,box=fitbox,profile=profile,ispsf=ispsf,extraparams=fitextra))
	names(tofit) = profile
	if(addsky) tofit$sky = list(bg=fitsky)
	tolog = list(profitMakeList(xcen=logical(ncomp),ycen=logical(ncomp),mag=logical(ncomp),
		size=logsize,	slope=logslope,ang=logang,axrat=logaxrat,box=logical(ncomp),profile=profile,ispsf=ispsf,
		extraparams=logextra))
	names(tolog) = profile

	if(unlog)
	{
		parm = unlist(modellist)
		tounlog = which(unlist(tolog))
		parm[tounlog] = 10^parm[tounlog]
		modellist = relist(parm, skeleton = modellist)
	}

	if(addsky) tolog$sky = list(bg=FALSE)
	alllims = list(magmin,magmax)
	if(!ispsf) alllims = c(list(xmin,xmax,ymin,ymax),alllims)
	if(!isps) alllims = c(alllims, list(sizemin,sizemax,slopemin,slopemax,angmin,angmax,axratmin,axratmax,boxmin,boxmax))
	for(par in extranames) alllims = c(alllims, extramin[[par]],extramax[[par]])
	sigmas = list(mag=magsd)
	if(!ispsf) sigmas = c(list(xcen=xsd,ycen=ysd),sigmas)
	if(!isps) sigmas = c(sigmas,list(size=sizesd,slope=slopesd,ang=angsd,axrat=axratsd,box=boxsd))
	sigmas = c(sigmas, structure(extrasd,names=extranames))
	means = list(mag=magmean)
	if(!ispsf) means = c(list(xcen=xmean,ycen=ymean),means)
	if(!isps) means = c(means,list(size=sizemean,slope=slopemean,ang=angmean,axrat=axratmean,box=boxmean))
	means = c(means, structure(extramean,names=extranames))

	indlims = all(ncomp==unlist(lapply(alllims,length)))
	samelims = all(1==unlist(lapply(alllims,length)))
	stopifnot(indlims || samelims)

	indsigmas = all(ncomp==unlist(lapply(sigmas,length)))
	samesigmas = all(1==unlist(lapply(sigmas,length)))
	stopifnot(indsigmas || samesigmas)

	indmeans = all(ncomp==unlist(lapply(means,length)))
	samemeans = all(1==unlist(lapply(means,length)))
	stopifnot(indmeans || samemeans)

	# This is hideous :(((
	if(samelims)
	{
		alllims = list(xcen=c(xmin,xmax),ycen=c(ymin,ymax),mag=c(magmin,magmax),size=c(sizemin,sizemax),
			slope=c(slopemin,slopemax), ang=c(angmin,angmax), axrat=c(axratmin,axratmax), box=c(boxmin,boxmax))
		alllims = lapply(alllims, function(x) { return(rep(list(lim=x),ncomp))})
		extralims = list()
		for(name in names(extramin)) extralims[[name]] = c(extramin[[name]],extramax[[name]])
		intervals = list(profitMakeList(xcen=alllims$xcen, ycen=alllims$ycen, mag=alllims$mag, size=alllims$size,
			slope=alllims$slope,ang=alllims$ang, axrat=alllims$axrat, box=alllims$box,profile=profile,ispsf=ispsf,
			extraparams = lapply(extralims, function(x) { return(rep(list(lim=x),ncomp))})))
		names(intervals) = profile
	} else {
		stop("Independent parameter limits not supported yet")
	}
	if(addsky) intervals$sky = list(bg=list(lim=c(skybgmin,skybgmax)))
	#print(sigmas)
	if(samesigmas)
	{
		sigmas = lapply(sigmas, function(x) { return(rep(x,ncomp))})
	}
	sigmas = profitMakeList(xcen=sigmas$xcen, ycen=sigmas$ycen, mag=sigmas$mag, size=sigmas$size,
		slope=sigmas$slope,ang=sigmas$ang, axrat=sigmas$axrat, box=sigmas$box, ispsf = ispsf, profile = profile,
		extraparams = extrasd)
	if(samemeans)
	{
		means = lapply(means, function(x) { return(rep(x,ncomp))})
	}
	if(addsky) sigmas$sky = sky = list(bg=skybgsd)
	means = profitMakeList(xcen=means$xcen, ycen=means$ycen, mag=means$mag, size=means$size,
		slope=means$slope,ang=means$ang, axrat=means$axrat, box=means$box, ispsf = ispsf, profile = profile,
		extraparams = extramean)
	if(addsky) means$sky = sky = list(bg=skybgmean)
	meanvec = unlist(means)
	tologvec = unlist(tolog)
	meanvec[tologvec] = 10^meanvec[tologvec]
	means = relist(flesh = meanvec, skeleton = means)
	tofitprior = tofit
	if(priortofitonly) tofitprior=NULL
	rv = list(modellist=modellist,
		tofit=tofit, tolog=tolog,
		intervals=intervals)
	if(ispsf) for(n in names(rv)) rv[[n]] = list(psf=rv[[n]])
	rv$priors = profitMakePriors(modellist=modellist,sigmas = sigmas,tolog = tolog, tofit=tofitprior, means = means, allowflat = TRUE)
	return(rv)
}

profitCombineLists <- function(lists)
{
	combined = list()
	lnames = c("modellist","tofit","tolog","intervals")
	priors = list(sigmas=list(),means=list(),tolog=list(),tofit=list())
	noadd = c("sky","psf")
	for(lname in lnames)
	{
		combined[[lname]] = list()
		for(mlist in lists)
		{
			stopifnot(all(lname %in% names(mlist)))
			compnames = names(mlist[[lname]])
			matched = compnames %in% names(combined[[lname]])
			for(comp in 1:length(compnames))
			{
				compname = compnames[comp]
				if(matched[comp])
				{
					if(!(compname %in% noadd))
					{
						for(subcomp in names(mlist[[lname]][[compname]]))
						{
							combined[[lname]][[compname]][[subcomp]] = c(combined[[lname]][[compname]][[subcomp]],
								mlist[[lname]][[compname]][[subcomp]])
						}
					}
				} else combined[[lname]][[compname]] = c(combined[[lname]][[compname]], mlist[[lname]][[compname]])
			}
		}
	}
	for(pname in names(priors))
	{
		for(mlist in lists)
		{
			newformals = relist(formals(mlist$priors)[[pname]], mlist$modellist)
			compnames = names(newformals)
			matched = compnames %in% names(priors[[pname]])
			for(comp in 1:length(compnames))
			{
				compname = compnames[comp]
				if(matched[comp])
				{
					if(!(compname %in% noadd))
					{
						for(subcomp in names(newformals[[compname]]))
						{
							priors[[pname]][[compname]][[subcomp]] = c(priors[[pname]][[compname]][[subcomp]],
								newformals[[compname]][[subcomp]])
						}
					}
				} else priors[[pname]][[compname]] = c(priors[[pname]][[compname]],newformals[[compname]])
			}
		}
		#priors[[pname]] = unlist(priors[[pname]])
	}
	combined$priors = profitMakePriors(modellist=combined$modellist,
			sigmas = priors$sigmas, tolog = priors$tolog, means = priors$means, tofit = priors$tofit,
			allowflat = TRUE)
	return(combined)
}

mlFit <- function(Data, init=Data$init, maxeval=1000, maxiter=100, algo.func="NEWUOA",
	cmasigma = formals(Data$priors)$sigmas[formals(Data$priors)$tofit], cmathreads=1,
	tolerance=1e-3, ftol_rel=1e-3, maxwalltime = Inf)
{
	stopifnot(is.numeric(maxwalltime))
	if(maxwalltime <= 0) maxwalltime = 1e-10
	tfinal = timeinmins() + maxwalltime
	algos = c("NEWUOA","BOBYQA","HAR","LM","NM","SPG")
	stopifnot(algo.func %in% algos || algo.func == "CMA")

  # Setup limits for CMA, which for now needs to know them since it doesn't take back the modified params from the fit function
  # We could modify CMA to do that...
  lims = unlist(Data$intervals)
  whichfit = which(unlist(Data$tofit))
  whichlog = which(unlist(Data$tolog))
  lims = list(
  	lower = lims[endsWith(names(lims),"lim1")],
  	upper = lims[endsWith(names(lims),"lim2")]
  )
  for(lim in names(lims))
  {
  	lims[[lim]][whichlog] = log10(lims[[lim]][whichlog])
  	lims[[lim]] = lims[[lim]][whichfit]
  }

  isnlopt = algo.func %in% c("NEWUOA","BOBYQA", "NM")
	if(isnlopt)
	{
		if(identical(algo.func,"NM")) algo.func = "NELDERMEAD"
		else if(identical(algo.func,"NEWUOA")) algo.func = "NEWUOA_BOUND"
		rv = list(rval=nloptr(x0 = init, eval_f=function(...) {return(-profitLikeModel(...)$LP)},
			opts = list(algorithm=paste0("NLOPT_LN_",algo.func),maxtime=maxwalltime*60,
			maxeval=maxeval,xtol_rel=0,ftol_abs=ftol_rel), lb=lims$lower, ub = lims$upper,
			Data=Data))
		rv$value = -rv$rval$objective
		rv$par = rv$rval$solution
		rv$rval$eval_f = get("eval_f",envir = environment(rv$rval$eval_f))
		environment(rv$rval$eval_f) = .GlobalEnv
		rv$rval$nloptr_environment = .GlobalEnv
		names(rv$par) = names(init)
	} else if(algo.func == "CMA") {
		Data$algo.func = "CMA"
		rv = list(rval=cmaeshpc(init, profitLikeModel, Data=Data, control=list(
    	maxit=maxiter, fnscale=-1.0, sigma=cmasigma, diag.sigma=TRUE, diag.eigen=TRUE, diag.pop=TRUE,
    	diag.value=TRUE, maxwalltime=maxwalltime, trace=TRUE, stopfitness = 0, stop.tolx=tolerance*cmasigma,
    	nthreadseval=cmathreads),	lower = lims$lower, upper=lims$upper))
		rv$value = return$rval$value
		rv$par = return$rval$par
	} else {
		Data$algo.func="LA"
		rv=list(rval=LaplaceApproximation(profitLikeModel, parm=Data$init, Data=Data,
			Iterations=maxiter, Method=algo.func, CovEst = "Identity",
			MaxWalltime = maxwalltime, sir=FALSE, CheckDataMatrixRanks = FALSE))
	}
  rv$overtime = exceededmaxtime(timeinmins(),tfinal)
	return(rv)
}

mlFitTest <- function(Data, init=Data$init, maxeval=1000, baseiter=100, doCMA=TRUE, cmaiter=baseiter,
	cmasigmafromcovar=FALSE, maxwalltime=Inf)
{
  # Setup limits for CMA, which for now needs to know them since it doesn't take back the modified params from the fit function
  # We could modify CMA to do that...
  lims = unlist(Data$intervals)
  whichfit = which(unlist(Data$tofit))
  whichlog = which(unlist(Data$tolog))
  lims = list(
  	lower = lims[endsWith(names(lims),"lim1")],
  	upper = lims[endsWith(names(lims),"lim2")]
  )
  for(lim in names(lims))
  {
  	lims[[lim]][whichlog] = log10(lims[[lim]][whichlog])
  	lims[[lim]] = lims[[lim]][whichfit]
  }
	inputinit = init

  #Data$verbose = TRUE
	# OK: HAR, HJ, LM, NM, SPG; bad: AGA, BFGS, CG, DFP, SR1; unuseable: BHHH, CG, SGD; both: NR, PSO, SOMA, TR
	algos = c("NEWUOA","BOBYQA","HAR","LM","NM","SPG")
	testiter=50
	res = list()
	for(algo in algos)
	{
		isnewuoa = algo == "NEWUOA"
		isbobyqa = algo == "BOBYQA"
		if(isnewuoa || isbobyqa)
		{
			Data$algo.func="CMA"
			mins = 0
			init = inputinit
			nloop = 10
			res[[algo]] = list(Deviance=numeric(nloop))
			t0 = proc.time()["elapsed"]
			for(i in 1:nloop)
			{
				if(isbobyqa)
				{
					nloptrv = bobyqa(init, lower=lims$lower, upper=lims$upper,
						fn = function(...) {return(-2*profitLikeModel(Data=Data,...))},
						control = list(maxeval=testiter))
				} else {
					nloptrv = newuoa(init,
						fn = function(...) {return(-2*profitLikeModel(Data=Data,...))},
						control = list(maxeval=testiter))
				}
				t1 = proc.time()["elapsed"]
				mins = mins + (t1-t0)/60
				res[[algo]]$Deviance[i] = nloptrv$value
				init = nloptrv$par
				t0 = proc.time()["elapsed"]
			}
			res[[algo]] = c(res[[algo]], list(Minutes=mins, Par=init))
		} else {
			Data$algo.func="LA"
			covest = "Identity"
			if(cmasigmafromcovar && algo == "LM") covest = "Hessian"
  		LAfit=LaplaceApproximation(profitLikeModel, parm=inputinit, Data=Data,
				Iterations=testiter*(1+5*(algo == "HAR") - 0.5*(algo=="SPG") - 0.25*(algo=="LM")),
				Method=algo, CovEst = covest && (algo == "HAR"), sir=FALSE, CheckDataMatrixRanks = FALSE)
  		res[[algo]] = LAfit
		}
	}
  if(doCMA)
  {
  	cmasigma = formals(Data$priors)$sigmas[whichfit]/2
  	Data$algo.func = "CMA"
  	if(!all(is.finite(cmasigma)) || cmasigmafromcovar)
  	{
			hess = hessian(f = function(...) {return(profitLikeModel(Data=Data,...))}, x0=inputinit, h = 1e-3)
			cmasigma = diag(inv(-hess))
			for(algo in algos) if(!is.null(res[[algo]]$Covar)) print(diag(res[[algo]]$Covar))
			print(cmasigma)
			cmasigma = cmasigma/sqrt(sum(cmasigma^2))
  	}

  	cmafit = cmaeshpc(inputinit, profitLikeModel, Data=Data, control=list(
    	maxit=cmaiter, fnscale=-1.0, sigma=cmasigma, diag.sigma=TRUE, diag.eigen=TRUE, diag.pop=TRUE,
    	diag.value=TRUE, maxwalltime=maxwalltime, trace=TRUE, stopfitness = 0, stop.tolx=5e-2*cmasigma),
    	lower = lims$lower, upper=lims$upper)
  }
	xlims = c(0,max(unlist(lapply(res,function(x) {return(x$Minutes)}))))
	ylims = range(unlist(lapply(res,function(x) {return(range(x$Deviance))})))
	if(doCMA)
	{
		cmaiter = dim(cmafit$diagnostic$value)[1]
		x = (1:cmaiter)/cmaiter*cmafit$diagnostic$walltime
		y = 2*rowMins(cmafit$diagnostic$value)
		ylims = range(c(ylims,-2*cmafit$value))
	} else {
		x = y = c()
	}
	Data$algo.func = "LD"
	initdev = profitLikeModel(inputinit,Data)$Dev
	if(initdev > ylims[2]) ylims[1] = initdev
	ylims = rev(ylims)
	magplot(x, y, lty=1, lwd=3, col="red", type="l",
		xlim = xlims, ylim = ylims, xaxs="i",yaxs="i")
	for(a in 1:length(res))
	{
		fit = res[[a]]
		algo = names(res)[a]
		ndev = length(fit$Deviance)
		if(algo == "NEWUOA" || algo == "BOBYQA") iter = ndev
		else iter = min(fit$Iterations,ndev)
		res[[a]]$lty = a
		res[[a]]$lwd = 1+(a%%2)
		lines((1:iter)/iter*fit$Minutes,fit$Deviance, lty=res[[a]]$lty, lwd=res[[a]]$lwd)
	}
	legs = algos
	ltys = unlist(lapply(res, function(x) { return(x$lty)}))
	lwds = unlist(lapply(res, function(x) { return(x$lwd)}))
	cols = rep("black",length(algos))
	if(doCMA)
	{
		legs = c("CMA",legs)
		cols = c("red",cols)
		ltys = c(1, ltys)
		lwds = c(3,lwds)
	}
	legend("topright", legend = legs, lty = ltys, lwd = lwds, col=cols)
	# TODO: Return something...
	# Actually, rewrite so that CMA is included in the list
}

convergeFit <- function(Data, init=Data$init, initalg = "CHARM", initspecs = list(alpha.star=0.44),
	finalalg=initalg, finalspecs=initspecs, inititer = 500, finaliter=inititer*4, maxruns = Inf,
	inititerharm=inititer, inititermcmc=inititer, maxeval=inititer, ftol_rel=0,
	domlfit = TRUE, domcmc = TRUE, docma = FALSE, cmasigma = NULL, cmasigmamult = c(2,1,0.5),
	cmatolfrac =0.05, cmaiter=inititer,	cmaresetmaxruns = TRUE, maxwalltime=Inf)
{
	tfinal = timeinmins() + maxwalltime
	overtime = FALSE
	bestLP = profitLikeModel(init,Data,makeplots = F)$LP
	converged = FALSE
	run = 0
	if(docma)
	{
		stopifnot(!domcmc)
		stopifnot(is.numeric(cmasigma) && (length(cmasigma) == length(init)))
		lims = unlist(Data$intervals)
	  whichfit = which(unlist(Data$tofit))
	  whichlog = which(unlist(Data$tolog))
	  lims = list(
	  	lower = lims[endsWith(names(lims),"lim1")],
	  	upper = lims[endsWith(names(lims),"lim2")]
	  )
	  for(lim in names(lims))
	  {
	  	lims[[lim]][whichlog] = log10(lims[[lim]][whichlog])
	  	lims[[lim]] = lims[[lim]][whichfit]
	  }
	}
	if(domlfit)
	{
		algofuncs = c("NEWUOA","BOBYQA","NM")
		fits = list()
		for(func in algofuncs)
		{
			fits[[func]] = mlFit(Data, init = init, maxeval=maxeval, algo.func=func,
				maxwalltime = tfinal-timeinmins())
			if(fits[[func]]$overtime) break
		}
		best = which.max(unlist(lapply(fits, function(x) { x$value })))
		fit = fits[[best]]
		init = fit$par
		bestLP = fit$value
		mlfunc = algofuncs[best]
		print(sprintf(paste0("Got bestfunc: ",mlfunc,"; value: %.6e; par:"),bestLP))
		print(init)
	}
	overtime = exceededmaxtime(timeinmins(),tfinal)
	while(!converged && !overtime)
	{
		if(domlfit)
		{
			fit = mlFit(Data, init = init, algo.func = mlfunc, maxeval = maxeval,
				ftol_rel = ftol_rel, maxwalltime = tfinal-timeinmins())
			newLP = fit$value
			if(newLP > bestLP)
			{
				# Re-apply intervals
				init = profitRemakeModellist(fit$par, Data=Data)$parm
			}
			print(paste0("mlfit from LP=",sprintf("%.5e",bestLP),
				" to LP=",sprintf("%.5e",newLP)," deltaLP=",sprintf("%.3e",newLP-bestLP),"; parm:"))
			print(fit$par)
		}
		else
		{
			fit = mlFit(Data=Data, init=init, maxiter=inititer, algo.func="HAR",
				maxwalltime=tfinal-timeinmins())
			newLP = fit$value
			if(newLP > bestLP) init = fit$par
		}
		run = run + 1
		converged = (newLP - bestLP) < exp(1) || (run >= maxruns)
		if(fit$overtime) break
		# Double check convergence if not using NM
		if(converged && !(domlfit && mlfunc == "NM"))
		{
			converged = FALSE
			domlfit = TRUE
			mlfunc = "NM"
		}
		if(converged)
		{
			if(newLP > bestLP) bestLP = newLP
			print(sprintf(paste0("MLFit converged at value: %.6e; par:"),bestLP))
			print(init)
			if(docma)
			{
				algofunc = Data$algo.func
				Data$algo.func = "CMA"

				for(mult in cmasigmamult)
				{
					# TODO: replace with MLFit call
					fit = cmaeshpc(init, profitLikeModel, Data=Data, control=list(
			    	maxit=cmaiter, fnscale=-1.0, sigma=mult*cmasigma,
			    	diag.sigma=TRUE, diag.eigen=TRUE, diag.pop=TRUE, diag.value=TRUE,
			    	maxwalltime=tfinal-timeinmins(),	trace=TRUE, stopfitness = Inf,
			    	stop.tolx=cmatolfrac*mult*cmasigma,
			    	lower = lims$lower, upper=lims$upper))
					if((fit$value - bestLP) > exp(1))
					{
						converged = FALSE
						fit$par = profitRemakeModellist(fit$par, Data=Data)$parm
						init = fit$par
						bestLP = fit$value
						print(sprintf(paste0("CMA converged at value: %.6e; par:"),bestLP))
						print(init)
						if(cmaresetmaxruns) run = 0
					} else docma = FALSE
					overtime = exceededmaxtime(timeinmins(),tfinal)
					if(overtime) break
				}
				Data$algo.func = algofunc
			}
			if(domcmc && !overtime)
			{
				print(paste0("LD convergence starting with params: c(",paste0(sprintf("%.6f",init),collapse=","),")"))
				# See if HARM can find a better solution
				fit = LaplacesDemon(profitLikeModel, Initial.Values=init, Data=Data, Iterations=inititerharm,
					Algorithm="HARM", Thinning=1, Specs=list(alpha.star=0.234, B=NULL),
					CheckDataMatrixRanks=FALSE, MaxWalltime=tfinal-timeinmins())
				best = which.max(fit$Monitor[,"LP"])
				newLP = fit$Monitor[best,"LP"]
				converged = (newLP - bestLP) < exp(1)
				if(newLP > bestLP)
				{
					init = fit$Posterior1[best,]
					bestLP = newLP
				}
				overtime = exceededmaxtime(timeinmins(),tfinal)
				if(overtime) break
				# If it's still converged, test it out with the final algorithm; otherwise try ML/CMA again
				if(converged)
				{
					print(paste0("LD MCMC starting with params: c(",paste0(sprintf("%.6f",init),collapse=","),")"))
					if(finalalg != "HARM")
					{
						fit = LaplacesDemon(profitLikeModel, Initial.Values=init, Data=Data, Iterations=inititermcmc,
							Algorithm=finalalg, Thinning=1, Specs=finalspecs, CheckDataMatrixRanks=FALSE,
							MaxWalltime=tfinal-timeinmins())
						best = which.max(fit$Monitor[,"LP"])
						newLP = fit$Monitor[best,"LP"]
					}
					converged = (newLP - bestLP) < exp(1)
					if(newLP > bestLP)
					{
						init = fit$Posterior1[best,]
						bestLP = newLP
					}
					overtime = exceededmaxtime(timeinmins(),tfinal)
					if(overtime) break
					# If a test run of the final algorithm is still converged, go for it
					# Usually we don't want to use the final algorithm for convergence if it's e.g. CHARM vs HARM
					if(converged)
					{
						fit = LaplacesDemon(profitLikeModel, Initial.Values=init, Data=Data, Iterations=finaliter,
							Algorithm=finalalg, Thinning=1, Specs=finalspecs, CheckDataMatrixRanks=FALSE,
							MaxWalltime=tfinal-timeinmins())
						best = which.max(fit$Monitor[,"LP"])
						newLP = fit$Monitor[best,"LP"]
						# Need to rerun if the final MCMC run wasn't converged (i.e. still needed burn-in)
						converged = (newLP - bestLP) < exp(1)
						if(newLP > bestLP)
						{
							init = fit$Posterior1[best,]
							bestLP = newLP
						}
						overtime = exceededmaxtime(timeinmins(),tfinal)
						if(overtime) break
					}
					print(paste0("LD MCMC finished (converged=",converged,") with params: c(",paste0(sprintf("%.6f",init),collapse=","),")"))
				}
			}
		}
		if(newLP > bestLP) bestLP = newLP
		overtime = exceededmaxtime(timeinmins(),tfinal)
	}
	if(!exists("fit"))
	{
		fit = list(par=init)
		Data$algo.func="LD"
		fit$value = profitLikeModel(fit$par,Data=Data)$LP
	}
	fit$converged = converged
	fit$overtime = overtime
	if(is.demonoid(fit) || is.laplace(fit))
	{
		fit$par = init
		fit$value = bestLP
	}
	fit$par = profitRemakeModellist(fit$par, Data=Data)$parm
	return(fit)
}

# finesample makes integration slightly more accurate, though it isn't
# entirely necessary because subpixel integration is adaptive by default
# returnfine should be true if you want to finesample the PSF for convolving
# a finesampled model
# cropfine will crop the finesampled psf so dim(psf)-1 is divisible by pixstep,
# and so the model padding is guaranteed to be divisible by pixstep as well
# this helps ensure an efficient image size for FFT convolution, but is
# unnecessary for brute force convolution
makeProfitPSFAutoSize <- function(modellist, finesample=1L, minpsfsum=0.999,
	pixstep=10,returnfine=FALSE,cropfine=returnfine)
{
	stopifnot("moffat" %in% names(modellist))
	stopifnot(minpsfsum < 1 & minpsfsum > 0)
	stopifnot(is.integer(finesample) && (finesample >= 1L))
	# Force returning an odd PSF
	if(returnfine) stopifnot((finesample %% 2) == 1 && (finesample > 1L))
	if(cropfine) stopifnot(returnfine)
	psfdim = pixstep*round(10*modellist$moffat$fwhm/pixstep)
	npsfs = length(psfdim)
	# Make sure the dimensions are odd, so that there's one extra pixel in the middle
	psfdim = psfdim + ((psfdim %% 2) == 0)
	if(npsfs > 1) psfdim = rep(max(psfdim),npsfs)
	modellist$moffat = c(list(xcen=psfdim/2,ycen=psfdim/2),modellist$moffat)
	psfim = profitMakeModel(modellist, dim = rep(psfdim,2),finesample = finesample)$z
	sumpsf = sum(psfim)
	hpixstep = pixstep/2
	while(sumpsf < minpsfsum)
	{
		psfdim = psfdim + pixstep
		modellist$moffat[["xcen"]] = modellist$moffat[["xcen"]] + hpixstep
		modellist$moffat[["ycen"]] = modellist$moffat[["ycen"]] + hpixstep
		psfim = profitMakeModel(modellist, dim = rep(psfdim,2),finesample = finesample)$z
		sumpsf = sum(psfim)
		stopifnot(sumpsf < (2-minpsfsum))
	}
	if(returnfine && (finesample>1))
	{
		psfim = profitMakeModel(modellist, dim = rep(psfdim,2),
			finesample = finesample, returnfine=TRUE)$z
		if(cropfine)
		{
			finecrop = (finesample-1)/2
			psfim = psfim[finecrop+(1:(dim(psfim)[1] - finecrop*2)),
				finecrop+(1:(dim(psfim)[2] - finecrop*2))]
		}
	}
	return(psfim)
}

# TOOD: Finish this
# The idea is to make a segim from the model fits
# TBD as to how it should be structured exactly
makeProfitComponentSegim <- function(parm, Data)
{
	modellist = profitRemakeModellist(parm, Data=Data)

	for(comps in names(modellist))
	{

	}
}

makeProfitSingleSourceLists <- function(profile = "sersic", x, y,
	mag, size, slope, ang, axrat, box=0, extra=list(),
	xmin, xmax, ymin=xmin, ymax=xmax,
	magmin = mag-5, magmax=mag+5, sizemin, sizemax,
	slopemin, slopemax,axratmin=0.02, axratmax=1,
	extramin=list(), extramax=list(),
	xsd, ysd, magsd, addsky=FALSE, fitsky=FALSE, skybg=0, skybgsd=0,
	sigmas=NULL, fitx = TRUE, fity=fitx, fitmag=TRUE, fitsize=TRUE,
	fitslope=TRUE, fitang=TRUE, fitaxrat=TRUE, fitbox=FALSE, fitextra=list(),
	unlog=FALSE)
{
	lists = profitMakeLists(profile=profile, xcen = x, ycen=y, mag = mag,
		size = size, slope = slope, ang = ang, axrat = axrat, box=box, extra = extra,
		xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax,
		magmin = magmin, magmax=magmax, sizemin = sizemin, sizemax=sizemax,
		slopemin = slopemin, slopemax=slopemax, axratmin=axratmin,
		axratmax=axratmax, extramin=extramin, extramax=extramax,
		addsky = addsky, skybg = skybg, fitsky = fitsky,
		xsd=xsd, ysd=ysd, magsd=magsd, fitx = fitx, fity=fity,
		fitmag=fitmag, fitsize=fitsize, fitslope=fitslope, fitang = fitang,
		fitaxrat = fitaxrat, fitbox = fitbox, fitextra = fitextra, unlog = unlog
	)
	allsiginf = is.null(sigmas)
	if(allsiginf) sigmas = rep(Inf,length(unlist(lists$modellist)))
	else stopifnot(length(sigmas) == length(unlist(lists$modellist)))
	sigmas = relist(sigmas,skeleton = lists$modellist)
	if(addsky && allsiginf) sigmas$sky$bg = skybgsd

	lists$priors = profitMakePriors(lists$modellist, sigmas, lists$tolog, means=lists$modellist,
		allowflat = TRUE, tofit = lists$tofit)
	environment(lists$priors) = .GlobalEnv
	return(list(lists=lists, sigmas=sigmas))
}

# A big messy code to do multi-component fits to a galaxy, (vaguely) intelligently
# adding components. Consider it an experimental first pass.

# Yes, I'm aware there's lots of repetition and copypasta here, and that's bad.
# It all needs to be refactored at some point. Probably what we need is some
# kind of "source" class that can be an arbitrary combination of components
# linked by a common centre. Perhaps components should be named so it's
# easy to switch the disk from Sersic to broken exponential, bar from
# Sersic to Ferrers, etc.

profitFitGalaxyComponents <- function(image, sigma, psfim,
	segim=matrix(0,dim(image)[1],dim(image)[2]), maxcomp=2,
	segobjs = c(0), gain_eff=1, chisqrtarg=1, chisqrminratio=1.05,
	init=list(x=dim(image)[1]/2,y=dim(image)[2]/2,
		size=sqrt(prod(dim(image))), mag=NaN,
		slope=2.5, ang=0, axrat=1), skylevel = 0, skyRMS = median(sigma,na.rm = TRUE),
	nthreads=1, maxiter=1e3, maxruns=Inf, domcmc=TRUE, finesample=1L,
	psffine=finesample>1L, autocrop=FALSE, adddiskfirst=TRUE, dobrokenexp=FALSE,
	galfit=NULL, maxwalltime=Inf, openclenvs = profitGetOpenCLEnvs(make.envs=TRUE),
	precisions=c("double"))
{
	tfinal = timeinmins() + maxwalltime
	dimimg = dim(image)

	if(any(segim != 0))
	{
		region = matrix(0,dim(image)[1],dim(image)[2])
		for(obj in segobjs) region[segim == obj] = 1
	}
	else region=1+segim

	if(!is.finite(init$mag)) init$mag = -2.5*log10(sum(image[region],na.rm = TRUE))
	maginit = init$mag - 2.5*log10(gain_eff)
	siglists = makeProfitSingleSourceLists(profile = "sersic",
		x=init$x, y=init$y, mag=maginit, size=init$size, slope = init$slope,
		ang = init$ang, axrat = init$axrat, xmin = 0,
		xmax=dimimg[2], ymax=dimimg[2],
		magmin = maginit-5, magmax=maginit+5,
		sizemin = 0.01, sizemax=5*sqrt(prod(dimimg))/2,
		slopemin = 0.1, slopemax=20, axratmin=0.02, axratmax=1, addsky = TRUE,
		skybg = skylevel*gain_eff, fitsky = FALSE, xsd=2, ysd=2, magsd=0.5,
		skybgsd = skyRMS*gain_eff
	)
	lists = siglists$lists

	#lists$intervals$sky = list(bg=list(lim=skylevel+skyRMS*c(-5,5)))
	psflo = psfim
	if(finesample > 1L)
	{
		if((finesample %% 2) == 1)
		{
			finepad = (finesample-1)/2
			stopifnot((finepad %% 1) == 0)
			psfdim = dim(psfim)
			psflo = matrix(0, psfdim[1]+2*finepad, psfdim[2]+2*finepad)
			psflo[(1:psfdim[1])+finepad,(1:psfdim[2])+finepad] = psfim
		}
		psflo = profitDownsample(psfim,finesample)
	}

	convmethods = "fft"
	if(prod(dim(image),dim(psflo)) < 1e8) convmethods = c(convmethods,"brute")
	if(profitHasOpenCL()) convmethods = c(convmethods,"opencl")

	Data=profitSetupData(image, sigma=sigma, region=region,
		modellist=lists$modellist, tofit=lists$tofit,
		tolog=lists$tolog, intervals=lists$intervals,
		psf=psflo, finesample = 1L,
		priors=lists$priors, algo.func='LD',like.func="norm",
		nbenchmark = 3L, benchconvmethods = convmethods,
		benchintmethods = profitAvailableIntegrators(),
		benchprecisions = precisions,
		verbose=FALSE, omp_threads = nthreads,
		benchopenclenvs = openclenvs
	)
	Data$gain = gain_eff
	Data$mon.names = c(Data$mon.names, chisq="chisq")
	#return(list(Data=Data))

	if(is.null(galfit) || !is.list(galfit) || length(galfit) == 0) galfit = list(data=list())
	# TBD: Will implement fitting a/the mock at some point
	fittypes = list(data=FALSE,mock=TRUE)
	LDFit = c()
	ncomp = 1

	# Do a quick single Sersic fit to determine if the image should be cropped (for autocrop)
	modelname = "ser_frpa_init"
	isnull = is.null(galfit$data[[modelname]])
	complete = !isnull && !is.null(galfit$data[[modelname]]$MLDone) && galfit$data[[modelname]]$MLDone
	if(!isnull)
	{
		mlfit = galfit$data[[modelname]]$MLFit
		Data$init = mlfit$par
	}
	if(!complete)
	{
		mlfit = convergeFit(Data, inititer = maxiter/4,	domcmc = FALSE, maxruns = maxruns,
			maxwalltime = tfinal-timeinmins())
		if(!is.null(Data$convopt$fft$fftwplan)) Data$convopt$fft$fftwplan = NULL
		galfit$data$ser_frpa_init = list(Data=Data,MLFit=mlfit,MLDone=!mlfit$overtime,DoMCMC=FALSE)
		if(mlfit$overtime) return(galfit)
	}
	init = mlfit$par
	names(init) = names(Data$init)

	# Autocrop a big image to smaller one for efficiency, using single-Sersic fit to estimate
	# the minimum required cutout size
	if(autocrop)
	{
		Data$algo.func = ""
		model = profitLikeModel(init, Data)$model$z
		dimmod = dimimg
		summod = 10^(-0.4*(init["sersic.mag"]-Data$magzero))
		tocrop = 0
		mindim = dim(psfim)/finesample
		mindim[mindim<100] = 100
		while(all(dimmod>=mindim) && (sum(model[tocrop+(1:dimmod[1]),tocrop+(1:dimmod[2])])/summod)>0.98)
		{
			tocrop = tocrop+10
			dimmod = dimmod-20
		}
		toolow = dimmod < mindim
		if(any(toolow))
		{
			dimmod[toolow] = mindim[toolow]
		}
		dimmod = dimmod + ((dimmod %% 2) == 1)
		tocrop = (dimimg - dimmod)/2
		stopifnot(all(tocrop >= 0))
		if(any(tocrop > 0))
		{
			xy = paste0(paste0("sersic.",c("x","y")),"cen")
			init[xy] = init[xy] - tocrop
			cropx = tocrop + (1:dimmod[1])
			cropy = tocrop + (1:dimmod[2])
			image = image[cropx,cropy]
			sigma = sigma[cropx,cropy]
			region = region[cropx,cropy]
			Data=profitSetupData(image, sigma=sigma, region=region,
				modellist=profitRemakeModellist(init,Data=Data)$modellist,
				tofit=lists$tofit, tolog=lists$tolog, intervals=lists$intervals,
				psf=psflo, finesample = 1L,
				priors=lists$priors, algo.func='LD',like.func="norm",
				nbenchmark = 3L, benchconvmethods = convmethods,
				benchintmethods = profitAvailableIntegrators(),
				benchprecisions = precisions,
				verbose=FALSE, omp_threads = nthreads,
				benchopenclenvs = openclenvs
			)
			init = Data$init
			Data$mon.names = c(Data$mon.names, chisq="chisq")
		} else Data$algo.func = "LD"
	}
	ndof = sum(Data$region)
	chisqredgal = profitLikeModel(init, Data)$Monitor["chisq"]/ndof
	# If the quick Sersic is good enough already, converge and flag it to do MCMC
	if(chisqredgal/chisqrtarg <= chisqrminratio)
	{
		modelname = "ser_frpa_1"
		isnull = is.null(galfit$data[[modelname]])
		complete = !isnull && !is.null(galfit$data[[modelname]]$MLDone) && galfit$data[[modelname]]$MLDone
		if(!isnull)
		{
			mlfit = galfit$data[[modelname]]$MLFit
			Data$init = mlfit$par
		}
		if(!complete)
		{
			mlfit = convergeFit(Data, init = mlfit$par, inititer = maxiter/4,
				domcmc = TRUE,	maxruns = 50, maxwalltime = tfinal-timeinmins())
			if(!is.null(Data$convopt$fft$fftwplan)) Data$convopt$fft$fftwplan = NULL
			galfit$data[[modelname]] = list(Data=Data,MLFit=mlfit,Complete=!overtime,DoMCMC=TRUE)
			if(mlfit$overtime) return(galfit)
		}
	}
	initpars = c("mag","re","axrat","ang","nser")

	# Do at least a B/D decomposition (even if single Sersic is ok) and add more components if not
	while(((chisqredgal/chisqrtarg > chisqrminratio) || ncomp < 2) && ncomp < maxcomp)
	{
		ncomp = ncomp + 1
		# Try adding an exponential disk with fixed position angle first
		adddisk = ncomp == 2 && (adddiskfirst || init["sersic.nser"] < 0.8)

		if(adddisk)
		{
			modelname = paste0("ser_fipa_",ncomp)
		} else {
			modelname = paste0("ser_frpa_",ncomp)
		}

		isnull = is.null(galfit$data[[modelname]])
		if(isnull)
		{
			initnames = names(Data$init)
			if(is.null(initnames)) initnames = names(init)
			stopifnot(is.character(initnames))
			inits = list()
			for(initpar in initpars)
			{
				inits[[initpar]] = unname(init[which(startsWith(initnames,paste0("sersic.",initpar)))])
			}
			if(adddisk)
			{
				fluxfrac = 0.95 - min(0.85,max(init["sersic.nser"],0))
				if(fluxfrac > 0.5) inits$re = c(inits$re-0.2, inits$re+0.1)
				else inits$re = c(inits$re, inits$re-0.2)
				inits$axrat = c(-0.05, inits$axrat)
			} else {
				fluxfrac=0.1
				inits$re = c(inits$re, mean(inits$re) - (ncomp==2))
				inits$axrat = c(inits$axrat, -0.05)
			}
			inits$mag = c(inits$mag - 2.5*log10(1-fluxfrac), -2.5*log10(fluxfrac*sum(10^(inits$mag/-2.5))))
			inits$axrat[inits$axrat > 0] = -0.01
			inits$nser = c(0.4+!adddisk*inits$nser, 0+0.3*!adddisk)
			inits$nser[inits$nser > 0.6] = 0.6
			inits$nser[inits$nser < -0.6] = -0.6
			inits$ang = c(inits$ang, inits$ang[ncomp-1]+90*!adddisk) %% 360
			for(initpar in c("re","axrat","nser")) inits[[initpar]] = 10^inits[[initpar]]

			postfix = ""
			if(ncomp > 2) postfix = "1"
			alltrue = rep(TRUE, ncomp)
			siglists = makeProfitSingleSourceLists(profile = "sersic",
				x=unname(rep(init[paste0("sersic.xcen",postfix)],ncomp)),
				y=unname(rep(init[paste0("sersic.ycen",postfix)],ncomp)),
				mag = inits$mag, size = inits$re, slope = inits$nser,
				ang = inits$ang, axrat = inits$axrat, box = rep(0,ncomp),
				xmin = 0, xmax=dimimg[1], ymax=dimimg[2],
				sizemin = max(0.02,0.1/finesample), sizemax = 5*sqrt(prod(dimimg))/2,
				slopemin = 0.1, slopemax=20, axratmin=0.02, axratmax=1,
				magmin = min(inits$mag)-25, magmax=max(inits$mag)+25,
				fitx = c(TRUE, rep(NA,ncomp-1)), fitmag = alltrue,
				fitsize=alltrue, fitslope = c(rep(TRUE,ncomp-1),!adddisk),
				fitang = c(rep(TRUE,ncomp-1),!adddisk), fitaxrat = alltrue,
				fitbox=!alltrue,  fitsky = FALSE, addsky = TRUE,
				skybg = skylevel*gain_eff, skybgsd = skyRMS*gain_eff,
				xsd=2, ysd=2, magsd=0.5
			)
			lists = siglists$lists
			rm("Data")
			gc()
			Data=profitSetupData(image, sigma=sigma, region=region,
				modellist=lists$modellist, tofit=lists$tofit,
				tolog=lists$tolog, intervals=lists$intervals, psf=psfim,
				priors=lists$priors, algo.func='LD',like.func="norm",
				nbenchmark = 3L, benchconvmethods = convmethods,
				benchintmethods = profitAvailableIntegrators(),
				benchprecisions = precisions,
				verbose=FALSE, omp_threads = nthreads,
				finesample=finesample, psffinesampled = TRUE,
				benchopenclenvs = openclenvs)
			Data$parm.names = names(Data$init)
			Data$mon.names = c(Data$mon.names, chisq="chisq")
			Data$gain = gain_eff
			if(finesample > 1L)
			{
				Data$lo=profitSetupData(image, sigma=sigma, region=region,
					modellist=lists$modellist, tofit=lists$tofit,
					tolog=lists$tolog, intervals=lists$intervals, psf=psflo,
					priors=lists$priors, algo.func='LD',like.func="norm",
					nbenchmark = 3L, benchconvmethods = convmethods,
					benchintmethods = profitAvailableIntegrators(),
					benchprecisions = precisions,
					verbose=FALSE, omp_threads = nthreads,
					finesample=1L,
					benchopenclenvs = openclenvs)
				Data$lo$parm.names = names(Data$lo$init)
				Data$lo$mon.names = c(Data$lo$mon.names, chisq="chisq")
				Data$lo$gain = gain_eff
			}
		} else {
			Data = galfit$data[[modelname]]$Data
		}
		#Data$xyoffsets = xyoffsets
		if(adddisk)
		{
			complete = !isnull && !is.null(galfit$data[[modelname]]$MLDone) && galfit$data[[modelname]]$MLDone
			if(!isnull)
			{
				mlfit = galfit$data[[modelname]]$MLFit
				Data$init = mlfit$par
			}
			if(!complete)
			{
				if(!isnull)
				{
					benchmarks = profitDataBenchmark(
						modellist = Data$modellist, calcregion = Data$calcregion,
						imgdim = dim(Data$image), finesample = Data$finesample, psf=Data$psf, fitpsf = FALSE,
						omp_threads = nthreads, nbenchmark = 3L, benchconvmethods = convmethods,
						benchintmethods = profitAvailableIntegrators(), benchprecisions = precisions,
						benchopenclenvs = openclenvs)
					Data = profitDataSetOptionsFromBenchmarks(Data, benchmarks)
					if(finesample > 1L)
					{
						benchmarks = profitDataBenchmark(
							modellist = Data$lo$modellist, calcregion = Data$lo$calcregion,
							imgdim = dim(Data$lo$image), finesample = Data$lo$finesample, psf=Data$lo$psf, fitpsf = FALSE,
							omp_threads = nthreads, nbenchmark = 3L, benchconvmethods = convmethods,
							benchintmethods = profitAvailableIntegrators(), benchprecisions = precisions,
							benchopenclenvs = openclenvs)
						Data$lo = profitDataSetOptionsFromBenchmarks(Data$lo, benchmarks)
					}
				}
				if(finesample > 1L)
				{
					mlfitf = convergeFit(Data$lo, inititer = maxiter/4, domcmc = FALSE,
															maxruns = maxruns, maxwalltime = tfinal-timeinmins())
					Data$init = mlfitf$par
				}
				overtime = exceededmaxtime(timeinmins(),tfinal)
				if(!overtime)
				{
					mlfit = convergeFit(Data, inititer = maxiter/4, domcmc = FALSE, maxruns = maxruns,
						maxwalltime = tfinal-timeinmins())
					# Try reversing component mags again
					if(!mlfit$overtime)
					{
						whichmag = which(startsWith(names(Data$init),"sersic.mag"))
						init = mlfit$par
						init[whichmag] = rev(init[whichmag])
						if(finesample > 1L)
						{
							initf = mlfitf$par
							initf[whichmag] = rev(initf[whichmag])
							mlfitnew = convergeFit(Data$lo, init=mlfit$par, inititer = maxiter/4, domcmc = FALSE,
								maxruns = maxruns, maxwalltime = tfinal-timeinmins())
							if((mlfitnew$value > mlfitf$value) && !mlfitnew$overtime)
							{
								mlfitnew = convergeFit(Data, init = mlfitnew$par, inititer = maxiter/4, domcmc = FALSE,
									maxruns = maxruns, maxwalltime = tfinal-timeinmins())
							}
						} else {
							mlfitnew = convergeFit(Data, init=init, inititer = maxiter/4, domcmc = FALSE,
								maxruns = maxruns, maxwalltime = tfinal-timeinmins())
						}
						if(mlfitnew$value > mlfit$value)
						{
							mlfit = mlfitnew
							Data$init = init
						}
						overtime = exceededmaxtime(timeinmins(),tfinal)
					}
				}
				if(!is.null(Data$convopt$fft$fftwplan)) Data$convopt$fft$fftwplan = NULL
				galfit$data[[modelname]] = list(Data=Data,MLFit=mlfit,MLDone=!overtime,DoMCMC=FALSE)
				if(overtime) return(galfit)
			}
		}
		if(domcmc) maxrunsml = Inf
		else maxrunsml=maxruns
		# Try a broken exponential
		if(adddisk && dobrokenexp)
		{
			bexpname = "serbexp_fipa"
			complete = !is.null(galfit$data[[bexpname]]) &&
				!is.null(galfit$data[[bexpname]]$MLDone) &&
				galfit$data[[bexpname]]$MLDone
			if(complete)
			{
				newData = galfit$data[[bexpname]]$Data
			}
			else
			{
				if(is.null(galfit$data[[bexpname]]$Data))
				{
					best = galfit$data$ser_fipa_2$MLFit$par
				}
				if(!is.null(galfit$data[[bexpname]]$LDFit))
				{
					best = getLDFitBest(galfit$data[[bexpname]])
				} else if(!is.null(galfit$data[[bexpname]]$MLFit)) {
					best = galfit$data[[bexpname]]$MLFit$par
					names(best) = names(galfit$data[[bexpname]]$Data$init)
				}
				if(!is.null(galfit$data[[bexpname]]$Data))
				{
					newData = galfit$data[[bexpname]]$Data
					newData$init = best
					benchmarks = profitDataBenchmark(
						modellist = newData$modellist, calcregion = newData$calcregion,
						imgdim = dim(newData$image), finesample = newData$finesample, psf = newData$psf, fitpsf = newData$fitpsf,
						nbenchmark = 1L, benchconvmethods = profitAvailableConvolvers(),
						benchintmethods = profitAvailableIntegrators(), benchopenclenvs = openclenvs,
						printbenchmark = TRUE)
					newData = profitDataSetOptionsFromBenchmarks(newData, benchmarks)
				} else {
					names(best) = unlist(lapply(names(best), function(x) { return(substr(x,8,nchar(x)))}))
					mag = -2.5*log10(sum(10^(-0.4*best["mag1"])+10^(-0.4*best["mag2"])))

					siglistsbexpb = makeProfitSingleSourceLists(
						profile = "sersic",
						x=unname(best["xcen1"]), y=unname(best["ycen1"]), mag = unname(best["mag1"]),
						size = unname(best["re1"]), slope = unname(best["nser1"]), ang=unname(best["ang1"]),
						axrat = unname(best["axrat1"]), box=0, xmin = 0, xmax=dimimg[1], ymax=dimimg[2],
						sizemin = max(0.02,0.1/finesample), sizemax = 5*sqrt(prod(dimimg))/2,
						slopemin = 0.1, slopemax=20, axratmin=0.02, axratmax=1,
						magmin = mag-25, magmax=mag+25,
						fitx = TRUE, fitmag = TRUE, fitsize=TRUE, fitslope = TRUE,
						fitang = TRUE, fitaxrat = TRUE, fitbox=FALSE,
						fitsky = FALSE, addsky = TRUE,
						skybg = skylevel*gain_eff, skybgsd = skyRMS*gain_eff,
						xsd=2, ysd=2, magsd=0.5, unlog=TRUE
					)
					siglistsbexpd = makeProfitSingleSourceLists(
						profile = "brokenexp",
						x=unname(best["xcen1"]), y=unname(best["ycen1"]), mag = unname(best["mag2"]),
						size = unname(best["re2"]), slope = unname(best["re2"])-0.3, ang=unname(best["ang1"]),
						axrat = unname(best["axrat2"]), box=0, extra=list(a=1,rb=unname(best["re2"])-0.3),
						xmin = 0, xmax=dimimg[1], ymax=dimimg[2],
						sizemin = max(0.02,0.1/finesample), sizemax = 5*sqrt(prod(dimimg))/2,
						slopemin = 0.1, slopemax=20, axratmin=0.02, axratmax=1,
						magmin = mag-25, magmax=mag+25, extramin = list(a=0,rb=1e-2), extramax=list(a=100,rb=dimimg[1]),
						fitx = FALSE, fitmag = TRUE, fitsize=TRUE, fitslope = TRUE,
						fitang = FALSE, fitaxrat = TRUE, fitbox=FALSE, fitsky = FALSE,
						fitextra = list(a=TRUE,rb=TRUE), addsky = TRUE,
						skybg = skylevel*gain_eff, skybgsd = skyRMS*gain_eff,
						xsd=2, ysd=2, magsd=0.5, unlog=TRUE
					)
					lists = profitCombineLists(list(siglistsbexpb$lists,siglistsbexpd$lists))

					constraints <- function(modellist) {
						modellist$brokenexp$xcen = modellist$sersic$xcen
						modellist$brokenexp$ycen = modellist$sersic$ycen
						modellist$brokenexp$ang = modellist$sersic$ang
						return(modellist)
					}
					environment(constraints) = .GlobalEnv

					# TODO: there should be a better way of doing this...
					# It's a verbose way of re-applying all of the relevant
					# limits and constraints to the input params so that
					# profitSetupData doesn't barf on axrat>1, say
					whichfit = which(unlist(lists$tofit))
					best = unlist(lists$modellist)[whichfit]
					whichlog = which(unlist(lists$tolog)[whichfit])
					best[whichlog] = log10(best[whichlog])

					remake = profitRemakeModellist(parm=best,
																				 modellist=lists$modellist, tofit = lists$tofit, tolog = lists$tolog,
																				 intervals = lists$intervals, constraints = constraints)
					lists$modellist = remake$modellist

					newData=profitSetupData(
						image, sigma=sigma, region=region,
						modellist=lists$modellist, tofit=lists$tofit,
						tolog=lists$tolog, intervals=lists$intervals, psf=psfim,
						priors=lists$priors, algo.func='LD',like.func="norm",
						nbenchmark = 3L, benchconvmethods = convmethods,
						benchintmethods = profitAvailableIntegrators(),
						benchprecisions = precisions,
						verbose=FALSE, omp_threads = nthreads,
						finesample=finesample, psffinesampled = TRUE,
						constraints = constraints, benchopenclenvs = openclenvs)
					newData$tolog$brokenexp$rb = TRUE
					newData$mon.names = c(newData$mon.names, chisq="chisq")
					newData$gain = gain_eff
				}

				if(finesample > 1L)
				{
					newData$lo=profitSetupData(
						image, sigma=sigma, region=region,
						modellist=newData$modellist, tofit=newData$tofit,
						tolog=newData$tolog, intervals=newData$intervals, psf=psflo,
						priors=newData$priors, algo.func='LD',like.func="norm",
						nbenchmark = 3L, benchconvmethods = convmethods,
						benchintmethods = profitAvailableIntegrators(),
						benchprecisions = precisions,
						verbose=FALSE, omp_threads = nthreads,
						finesample=1L,
						constraints = newData$constraints, benchopenclenvs = openclenvs)
					newData$lo$mon.names = c(newData$mon.names, chisq="chisq")
					newData$lo$gain = gain_eff

					mlfit = convergeFit(newData$lo, init=newData$init, inititer = maxiter/4, maxruns=maxrunsml,
						domcmc = FALSE, docma = FALSE, cmasigma = NULL, cmaiter = 200,
						cmaresetmaxruns = TRUE, cmasigmamult = 1, cmatolfrac = 0.01,
						maxwalltime = tfinal-timeinmins())
					newinit = mlfit$par
				} else {
					newinit = newData$init
				}

				mlfit = convergeFit(newData, init = newinit, inititer = maxiter/4, maxruns=maxrunsml,
					domcmc = FALSE, docma = FALSE, cmasigma = NULL, cmaiter = 200,
					cmaresetmaxruns = TRUE, cmasigmamult = 1, cmatolfrac = 0.01,
					maxwalltime = tfinal-timeinmins())
				if(!is.null(newData$convopt$fft$fftwplan)) newData$convopt$fft$fftwplan = NULL
				galfit$data[[bexpname]] = list(Data=newData,MLFit=mlfit,MLDone=!mlfit$overtime,DoMCMC=TRUE)
				if(mlfit$overtime) return(galfit)
			}
			fipabexp = bexpname
			bexpname = "serbexp_frpa"
			complete = !is.null(galfit$data[[bexpname]]) &&
				!is.null(galfit$data[[bexpname]]$MLDone) &&
				galfit$data[[bexpname]]$MLDone
			if(!complete)
			{
				newData = galfit$data[[fipabexp]]$Data
				newData$constraints <- function(modellist) {
					modellist$brokenexp$xcen = modellist$sersic$xcen
					modellist$brokenexp$ycen = modellist$sersic$ycen
					return(modellist)
				}
				environment(newData$constraints) = .GlobalEnv
				newData$tofit$brokenexp$ang = TRUE
				if(!is.null(galfit$data[[fipabexp]]$LDFit)) best = getLDFitBest(galfit$data[[fipabexp]])
				else {
					best = galfit$data[[fipabexp]]$MLFit$par
					names(best) = names(galfit$data[[fipabexp]]$Data$init)
				}
				newData$init = unlist(newData$modellist)[which(unlist(newData$tofit))]
				newData$init[names(best)] = best
				newData$parm.names = names(newData$init)

				if(finesample > 1L)
				{
					newData$lo$constraints = newData$constraints
					newData$lo$tofit$brokenexp$ang = TRUE
					newData$lo$init = newData$init
					mlfit = convergeFit(newData$lo, inititer = maxiter/4, maxruns=maxrunsml,
						domcmc = FALSE, maxwalltime = tfinal-timeinmins())
					newData$init = mlfit$par
				}

				overtime = exceededmaxtime(timeinmins(), tfinal)
				{
					mlfit = convergeFit(newData, inititer = maxiter/4, maxruns=maxrunsml,
						domcmc = FALSE, maxwalltime = tfinal-timeinmins())
					# Try reversing component mags again
					if(!mlfit$overtime)
					{
						init = mlfit$par
						names(init) = names(newData$init)
						whichmag = c("sersic.mag","brokenexp.mag")
						init[whichmag] = rev(init[whichmag])
						mlfitnew = convergeFit(newData, init = init, inititer = maxiter/4,
							maxruns=maxrunsml, domcmc = FALSE, maxwalltime = tfinal-timeinmins())
						if(mlfitnew$value > mlfit$value)
						{
							mlfit = mlfitnew
							newData$init = init
						}
						mlfit$overtime = mlfitnew$overtime
					}
				}
				if(!is.null(newData$convopt$fft$fftwplan)) newData$convopt$fft$fftwplan = NULL
				galfit$data[[bexpname]] = list(Data=newData,MLFit=mlfit,MLDone=!mlfit$overtime,DoMCMC=TRUE)
				if(mlfit$overtime) return(galfit)
			}
		}
		if(adddisk)
		{
			modelname = paste0("ser_frpa_",ncomp)
			isnull = is.null(galfit$data[[modelname]])
			if(isnull)
			{
				init = galfit$data[[paste0("ser_fipa_",ncomp)]]$MLFit$par
				# Now let the position angle and Sersic index go free and update init
				for(par in c("nser","ang"))
				{
					Data$tofit$sersic[[par]][ncomp] = TRUE
					modpar = paste0("sersic.",par)
					last = which(startsWith(names(init),modpar))[ncomp-1]
					val = Data$modellist$sersic[[par]][ncomp]
					if(Data$tolog$sersic[[par]][ncomp-1]) val = log10(val)
					init = c(init[1:last],val,init[(1+last):length(init)])
					names(init)[last+1] = paste0(modpar,ncomp)
				}
				Data$init = init
				Data$parm.names = names(Data$init)
				for(var in c("tofit","init","parm.names")) Data$lo[[var]] = Data[[var]]
			}
		}
		complete = !isnull && !is.null(galfit$data[[modelname]]$MLDone) && galfit$data[[modelname]]$MLDone
		if(!isnull)
		{
			Data = galfit$data[[modelname]]$Data
			mlfit = galfit$data[[modelname]]$MLFit
			Data$init = mlfit$par
		} else if(!adddisk) init = Data$init
		if(!isnull || identical(Data$convopt$convolver,new("externalptr")))
		{
			benchmarks = profitDataBenchmark(
				modellist = Data$modellist, calcregion = Data$calcregion,
				imgdim = dim(Data$image), finesample = Data$finesample, psf=Data$psf, fitpsf = FALSE,
				omp_threads = nthreads, nbenchmark = 3L, benchconvmethods = convmethods,
				benchintmethods = profitAvailableIntegrators(), benchprecisions = precisions,
				benchopenclenvs = openclenvs)
			Data = profitDataSetOptionsFromBenchmarks(Data, benchmarks)
			if(finesample > 1L)
			{
				benchmarks = profitDataBenchmark(
					modellist = Data$lo$modellist, calcregion = Data$lo$calcregion,
					imgdim = dim(Data$lo$image), finesample = Data$lo$finesample, psf=Data$lo$psf, fitpsf = FALSE,
					omp_threads = nthreads, nbenchmark = 3L, benchconvmethods = convmethods,
					benchintmethods = profitAvailableIntegrators(), benchprecisions = precisions,
					benchopenclenvs = openclenvs)
				Data$lo = profitDataSetOptionsFromBenchmarks(Data$lo, benchmarks)
			}
		}
		if(!complete)
		{
			cmasigma = c(rep(0.5,2),rep(0.1,3*ncomp), rep(2,ncomp), rep(0.1,ncomp))
			if(domcmc) maxrunsml = Inf
			else maxrunsml=maxruns
			if(finesample > 1L)
			{
				mlfit = convergeFit(Data$lo, inititer = maxiter/4, domcmc = FALSE,
					maxruns = maxruns, maxwalltime = tfinal-timeinmins())
				Data$init = mlfit$par
			}
			mlfit = convergeFit(Data, inititer = maxiter/4, maxruns=maxrunsml,
				domcmc = FALSE, docma = TRUE, cmasigma = cmasigma, cmaiter = 200,
				cmaresetmaxruns = TRUE, cmasigmamult = 1, cmatolfrac = 0.01,
				maxwalltime = tfinal-timeinmins())
			if(!is.null(Data$convopt$fft$fftwplan)) Data$convopt$fft$fftwplan = NULL
			galfit$data[[modelname]] = list(Data=Data,MLFit=mlfit,MLDone=!mlfit$overtime,DoMCMC=TRUE)
			if(mlfit$overtime) return(galfit)
		} else {
			mlfit = galfit$data[[modelname]]$MLFit
		}
		init = profitRemakeModellist(mlfit$par, Data=Data)$parm
		if(identical(Data$convopt$convolver,new("externalptr"))) Data$convopt$convolver = NULL
		chisqredgal = profitLikeModel(init, Data)$Monitor["chisq"]/ndof
	}
	if(domcmc)
	{
		for(fittype in names(fittypes)[1])
		{
			for(modtype in names(galfit$data))
			{
				stopifnot(!is.null(galfit$data[[modtype]]$DoMCMC))
				if(galfit$data[[modtype]]$DoMCMC)
				{
					complete = !is.null(galfit$data[[modtype]]$MCMCDone) && galfit$data[[modtype]]$MCMCDone
					if(!complete)
					{
						Data = galfit$data[[modtype]]$Data
						Data$init = galfit$data[[modtype]]$MLFit$par
						# Make a mock of the best fit. TODO: Re-implement this if it seems compelling enough
						if(fittypes[[fittype]])
						{
							bestmodel = profitMakeModel(
								profitRemakeModellist(Data$init,Data = Data)$modellist,
								magzero = Data$magzero,	psf = Data$psf,
								dim = dim(Data$image), finesample = Data$finesample
							)
							# Todo: fix this at some point; but more likely just give up
							# for images with variable exposure times
							skygain = gain_eff
							Data$image = profitPoissonMonteCarlo(bestmodel$z + skylevel*skygain)
							Data$sigma = sqrt(Data$image)
							Data$image = Data$image - skylevel*skygain
						}
						benchmarks = profitDataBenchmark(
							modellist = Data$modellist, calcregion = Data$calcregion,
							imgdim = dim(Data$image), finesample = Data$finesample, psf=Data$psf, fitpsf = FALSE,
							omp_threads = nthreads, nbenchmark = 3L, benchconvmethods = convmethods,
							benchintmethods = profitAvailableIntegrators(), benchprecisions = precisions,
							benchopenclenvs = openclenvs)
						Data = profitDataSetOptionsFromBenchmarks(Data, benchmarks)
						parnames = names(Data$init)
						if(!identical(Data$parm.names,parnames))
						{
							stopifnot(is.character(parnames) && identical(length(parnames),length(Data$init)))
							Data$parm.names = parnames
						}
						LDFit = convergeFit(Data, init=Data$init, inititer = maxiter/2,
							finaliter = maxiter, domlfit = TRUE, domcmc = TRUE,
							inititermcmc = maxiter/5, maxwalltime = tfinal-timeinmins())
						overtime = exceededmaxtime(timeinmins(),tfinal)
						galfit[[fittype]][[modtype]]$LDFit = LDFit
						galfit[[fittype]][[modtype]]$MCMCDone = !overtime
						if(overtime) return(galfit)
					}
				}
			}
		}
	}
	return(galfit)
}

getLDFitBest <- function(fitdata=NULL, Data = fitdata$Data, Fit = fitdata$LDFit)
{
	return(which.max(Fit$Monitor[,"LP"]))
}

getLDFitParm <- function(fitdata=NULL, Data = fitdata$Data, Fit = fitdata$LDFit,
	which = getLDFitBest(fitdata,Data,Fit))
{
	return(Fit$Posterior1[which,1:Fit$Parameters])
}

getLDFitMon <- function(fitdata=NULL, Data = fitdata$Data, Fit = fitdata$LDFit,
	which = getLDFitBest(fitdata,Data,Fit))
{
	return(Fit$Monitor[which,])
}

profitLikeModelLDFit <- function(fitdata=NULL, Fit=fitdata$LDFit, Data=fitdata$Data,
	parm = getLDFitParm(fitdata = fitdata, Fit=Fit, Data=Data), ...)
{
	return(profitLikeModel(parm,Data, ...))
}

profitBenchmarkCutoutFFT <- function(finesample = 3L)
{
	neval = c(900,2*5*7*11,2*5*7*13,4*5*5*7,4*5*7*7,1000,2*5*7*17,8*3*5*7,2*3*5*5*7,32*5*7,1100)
	testevals = c(neval-2,neval-1,neval,neval+1,neval+2)
	psf = profitMakeGaussianPSF(fwhm=3*finesample,dim = c(31,31)*finesample)
	times = numeric(length(testevals))
	for(i in 1:length(testevals))
	{
		n = testevals[i]*3
		tmp = matrix(rnorm(n^2),n,n)
		rv = profitBenchmarkConv(tmp, psf, nbench = 10L, methods = c("FFTWconv"))$best$time
		print(rv)
		times[i] = rv
	}
	return(list(x=testevals,y=times))
}

readkey <- function(str="")
{
    cat(paste0(str,"Press [q + Enter] to quit, [Enter] to continue."))
    line <- readline()
    if(tolower(line) == "q")
    {
    	opt <- options(show.error.messages = FALSE)
  		on.exit(options(opt))
  		stop()
    }
}
