cutoutSourceWCS <- function(image, header, ra, dec, box,
	extras=list())
{
	stopifnot(!(("image" %in% names(extras)) || ("xy" %in% names(extras))))
	rv = list(maps=list(
		image = magcutoutWCS(image, header=header,
			loc=c(ra,dec), loc.type =  c("coord","image"),
			box=box)))
	rv$xy = radec2xy(ra,dec,rv$maps$image$header)
	for(extra in names(extras))
	{
		rv$maps[[extra]] = magcutoutWCS(extras[[extra]],
			loc=rv$maps$image$loc.orig, box=box,
			loc.type = rep("image",2))$image
	}
	return(rv)
}

# Nominally, the uncertainty on a median should be 1/(4*N*f(m)^2)
# But this would typically give below one count, e.g.:
# 1/(4*prod(dim(sigma))*dnorm(0,sd=median(sigma)*gain)^2)
# N = 1000^2, sd=100 (1e4 sky photons) gives skysd=0.0157 counts!
fitPointSources <- function(image, sigma, header, mask,
	segim, segstats, fitobj, gain=1, pixscale=1, maxcomp=1,
	skybgsd=sqrt(median(sigma,na.rm = TRUE)*gain),
	mincutoutsize=50, maxiter=1e3, nthreads=1L, domcmc=TRUE,
	fitsky=FALSE, fits=list())
{
	objn=1
	nobj = length(fitobj)
	for(obj in fitobj)
	{
		row = which(segstats$segID == obj)
		ra=segstats[row,"RAcen"]
		dec=segstats[row,"Deccen"]
		mag=-2.5*log10(segstats[row,"flux"])
		fwhm=2*segstats[row,"R50"]/pixscale
		con=max(1.2,1/sqrt(segstats[row,"con"]))
		axrat=segstats[row,"axrat"]
		ang=segstats[row,"ang"]

		cutoutsize = 5*fwhm
		if(cutoutsize < mincutoutsize) cutoutsize = mincutoutsize
		cutoutsize = c(cutoutsize,cutoutsize)

		cutouts = cutoutSourceWCS(image, header, ra, dec, box=cutoutsize,
			extras=list(segim=segim,sigma=sigma))
		cutouts$maps$segim[!is.finite(cutouts$maps$segim)] = Inf
		xy = cutouts$xy

		maginit = mag -2.5*log10(gain)
		lists = profitMakeLists(xcen = xy[1], ycen=xy[2], mag = maginit,
			size = fwhm, slope = con, ang = ang, axrat = axrat, profile="moffat",
			xmin = 0, xmax=cutoutsize[1], magmin = maginit-5, magmax=maginit+5,
			sizemin = 0.1,sizemax=cutoutsize[1], slopemin = 1.01,slopemax=10,
			axratmin=0.1,axratmax=1, addsky = TRUE, skybg = 0, fitsky = fitsky)

		sigmas = Inf+0*unlist(lists$modellist)
		sigmas[length(sigmas)] = skybgsd
		sigmas = relist(sigmas,skeleton = lists$modellist)

		priors = profitMakePriors(lists$modellist, sigmas, lists$tolog, means=lists$modellist, allowflat = TRUE)

		#lists$intervals$sky = list(bg=list(lim=skylevel+skyRMS*c(-5,5)))

		Data=profitSetupData(cutouts$maps$image$image*gain,
			sigma=cutouts$maps$sigma*gain, priors=priors,
			region=(cutouts$maps$segim == obj | cutouts$maps$segim <= 0),
			modellist=lists$modellist, tofit=lists$tofit, tolog=lists$tolog,
			intervals=lists$intervals, algo.func='LD',like.func="norm",
			verbose=FALSE,omp_threads = nthreads)
		Data$gain = gain
		Data$header = cutouts$maps$image$header
		objname = as.character(obj)

		print(sprintf("Fitting object %i/%i (id=%s), init=",objn,nobj,objname))
		print(Data$init)

		if(is.null(fits[[objname]])) fits[[objname]] = list()
		fittypes = list(data=FALSE,mock=TRUE)
		LDFit = c()
		for(fittype in names(fittypes)[1])
		{
			if(is.null(fits[[objname]][[fittype]])) fits[[objname]][[fittype]] = list()
			if(is.null(fits[[objname]][[fittype]]$moffat1))
			{
				fits[[objname]][[fittype]]$moffat1=list(Data=Data)
			}
			if(is.null(fits[[objname]][[fittype]]$moffat1$LDFit))
			{
				# TODO: fix this
				if(fittypes[[fittype]])
				{
					best = getLDFitBest(LDFit)
					Data$init = LDFit$Posterior1[best,]
					bestmodel = profitMakeModel(
						profitRemakeModellist(Data$init,Data = Data)$modellist,
						magzero = Data$magzero,	psf = Data$psf,
						dim = dim(Data$image), finesample = Data$finesample
					)
					Data$image = profitPoissonMonteCarlo(bestmodel$z + skyflux*skygain)
					Data$sigma = sqrt(Data$image)
					Data$image = Data$image - skyflux*skygain
				}
				fits[[objname]][[fittype]]$moffat1$LDFit = convergeFit(
					Data, inititer = maxiter/2, finaliter=maxiter, domcmc=domcmc)
			}
			if(maxcomp > 1)
			{
				if(is.null(fits[[objname]][[fittype]]$moffat2)) fits[[objname]][[fittype]]$moffat2 = list()
				if(is.null(fits[[objname]][[fittype]]$moffat2$LDFit))
				{
					init = fits[[objname]][[fittype]]$moffat1$LDFit$par
					bt = c(TRUE,TRUE)
					lists = profitMakeLists(xcen=rep(unname(init["moffat.xcen"]),2),
						ycen=rep(unname(init["moffat.ycen"]),2),
						mag=rep(unname(init["moffat.mag"]-2.5*log10(0.5)),2),
						size=rep(unname(init["moffat.fwhm"]),2) + c(0.05,-0.05),
						slope = rep(unname(init["moffat.con"]),2) + c(0,0.3),
						ang = rep(unname(init["moffat.ang"]),2),
						axrat = rep(unname(init["moffat.axrat"]),2), box = c(0,0),
						profile="moffat",	unlog=TRUE, addsky = TRUE,
						skybg = unname(init["sky.bg"]), fitx = c(TRUE,NA),
						fitmag = c(TRUE,TRUE), fitsize = bt, fitslope=bt, fitang=bt,
						fitaxrat=bt, fitbox=!bt, fitsky=TRUE)
					sigmas = Inf+0*unlist(lists$modellist)
					sigmas[length(sigmas)] = skybgsd
					sigmas = relist(sigmas,skeleton = lists$modellist)
					priors = profitMakePriors(lists$modellist, sigmas, lists$tolog, means=lists$modellist, allowflat = TRUE)
					environment(priors) = new.env()
					Data=profitSetupData(cutouts$maps$image$image*gain,
						sigma=cutouts$maps$sigma*gain, priors=priors,
						segim=(cutouts$maps$segim == obj | cutouts$maps$segim <= 0),
						modellist=lists$modellist, tofit=lists$tofit, tolog=lists$tolog,
						intervals=lists$intervals, algo.func='LD',like.func="norm",
						verbose=FALSE,omp_threads = nthreads)
					Data$gain = gain
					Data$header = cutouts$maps$image$header
					LDFit = convergeFit(Data, inititer = maxiter/2, finaliter=maxiter, domcmc=domcmc)
					fits[[objname]][[fittype]]$moffat2=list(Data=Data,LDFit=LDFit)
				}
			}
		}
		objn = objn+1
	}
	return(fits)
}

# TODO: Use the same kind of return data structure as
# the galaxy fitting procedure, overtime checks, etc.
fitPSF <- function(image, sigma, hdr, mask, segim, psfits,
	pointsources, means, sds, skypix = segim <= 0, init = means[1,],
	maxcomp=1, skybg=0, gain_eff=1, finesample=1L,nthreads=1L,
	inititer=500, finaliter=2000, docma=TRUE, fitboxmulti=FALSE,
	fits = list(), plot=FALSE,	verbose=FALSE)
{
	if(maxcomp > 2) maxcomp = 2
	posmask = skypix
	for(id in pointsources) posmask = posmask | (segim == id)
	mask = mask | !posmask
	magimage(image*!mask)
	allt = !logical(length(pointsources))
	#init = means[]
	psflists = profitMakeLists(xcen=0,ycen=0,mag=0,	size=unname(init["moffat.fwhm"]),
		slope = unname(init["moffat.con"]), ang = unname(init["moffat.ang"]),
		axrat = unname(init["moffat.axrat"]), ispsf = TRUE, profile="moffat",
		boxmin = -0.2, boxmax=0.2, fitmag = FALSE, fitbox = TRUE, unlog=TRUE)
	constraints = NULL
	fitname = paste0("moffat",maxcomp)
	if(is.null(fits[[fitname]]$Data))
	{
		for(fitps in c(FALSE,TRUE))
		{
			pslists = profitMakeLists(xcen = means[,"moffat.xcen"],
				ycen=means[,"moffat.ycen"], mag=means[,"moffat.mag"]-2.5*log10(gain_eff),
				profile="pointsource", fitx=fitps&allt, fity = fitps&allt, fitmag=fitps&allt,
				xsd = 3*sds[,"moffat.xcen"], ysd = 3*sds[,"moffat.ycen"],
				magsd = 3*sds[,"moffat.mag"], unlog=TRUE,
				skybg = skybg*gain_eff, addsky=TRUE, fitsky = TRUE)
			lists = profitCombineLists(list(psf=psflists, pointsource=pslists))

			Data = profitSetupData(image=image*gain_eff,sigma = sigma*gain_eff, mask=mask,
				modellist = lists$modellist, tofit = lists$tofit, tolog=lists$tolog, intervals = lists$intervals,
				priors=lists$priors, constraints = constraints, algo.func="LD", like.func="norm",
				verbose=(verbose && !fitps), omp_threads=nthreads)
			Data$gain = gain_eff
			if(!fitps) init = Data$init
			else if(maxcomp == 1) Data$init[names(init)] = init

			if(plot) {
				rv = profitLikeModel(Data$init, Data, makeplots = T, plotchisq = T)
				print(rv$Monitor)
			}

			if(!fitps)
			{
				cmasigma = formals(lists$priors)$sigmas[formals(lists$priors)$tofit]/3
				extrapar = c("fwhm","con","ang","axrat")
				cmasigma[c(paste0("psf.moffat.",c(extrapar,"box")),"sky.bg")] = c(median(sds[,c(paste0("moffat.",extrapar))]),
					0.1,median(sds[,"sky.bg"]))
			} else {
				moffatpars = extrapar
				if(!((maxcomp > 1) && !fitboxmulti)) moffatpars = c(extrapar,"box")
				cmasigma = c(rep(0.1,each=maxcomp-1),rep(cmasigma[paste0("psf.moffat.",
						moffatpars)],each=maxcomp),
					as.vector(sds[,paste0("moffat.",c("xcen","ycen","mag"))]),
					cmasigma["sky.bg"])
			}
			if(fitps && (!is.null(fits$par) && (length(fits$par)) == length(Data$init))) init = fits$par
			else
			{
				fit = convergeFit(Data, init=Data$init, inititer = inititer, finaliter=finaliter,
					maxeval=inititer, ftol_rel = (1e-6)*(1-fitps), domcmc = FALSE, domlfit = TRUE,
					docma=docma, cmasigma = cmasigma, cmasigmamult = c(1), cmaiter = 100)
				init = fit$par
			}
			if(is.null(names(init))) names(init) = names(Data$init)
			if(maxcomp > 1 && !fitps)
			{
				if(is.null(fits$moffat1)) fits$moffat1 = list(Data=Data, MLFit = fit, pointsources=pointsources)
				psflists = profitMakeLists(xcen=c(0,NA),ycen=c(0,NA), mag=rep(-2.5*log10(0.5),2),
					size=rep(unname(init["psf.moffat.fwhm"]),2) + c(0,-0.1),
					slope = rep(unname(init["psf.moffat.con"]),2) + c(0,0.3), ang = rep(unname(init["psf.moffat.ang"]),2),
					axrat = rep(unname(init["psf.moffat.axrat"]),2), box=rep(0,2), #box = rep(unname(init["psf.moffat.box"]),2),
					ispsf = TRUE, profile="moffat", fitmag = c(TRUE,FALSE), fitbox = rep(fitboxmulti,2),
					magmin = -2.5*log10(0.9999), magmax = -2.5*log10(1-0.9999), boxmin=-0.2, boxmax=0.2, unlog=TRUE)
				constraints <- function(modellist) {
					modellist$psf$moffat$mag[2] = -2.5*log10(1-10^(-0.4*modellist$psf$moffat$mag[1]))
					return=modellist
				}
			}
		}
		Data$init = init
		fits[[fitname]] = list(Data=Data,pointsources=pointsources)
	} else {
		for(var in c("xcen","ycen","mag")) fits[[fitname]]$Data$modellist[[var]] =
			fits[[fitname]]$Data$modellist[[var]] | TRUE
	}
	fits[[fitname]]$Data$algo.func = "LD"
	fits[[fitname]]$LDFit = convergeFit(fits[[fitname]]$Data, inititer = inititer,
		finaliter=finaliter, initalg = "HARM", initspecs = list(alpha.star=0.234, B=NULL))
	return(fits)
}

fitGalaxies <- function(process, psffit, fitobj, band, header,
	fits=list(), maxcomp=2, gain=1, mincutoutsize=36, finesample=1L,
	nthreads=1L, maxiter=2e3L, domcmc = (length(fitobj) == 1),
	pixscale=1, autocropband = band, dobrokenexp=FALSE, test=FALSE)
{
	stopifnot(band %in% names(process$single))
	image = psffit$Data$image
	mask = process$multi$mask
	sigma = psffit$Data$sigma
	segim = psffit$Data$segim
	# or perhaps:
	segim = process$multi$proc$segim
	# For individual galaxy fits
	segim[mask] = Inf
	segstats = process$single[[band]]$stats

	ndata = sum(psffit$Data$region)
	psfhaschisq = "chisq" %in% psffit$Data$mon.names
	if(psfhaschisq)
	{
		bestpsf = getLDFitBest(psffit)
		chisqredpsf = getLDFitMon(psffit, which=bestpsf)["chisq"]/ndata
		bestpsf = getLDFitParm(psffit, which=bestpsf)
	} else {
		psffit$Data$mon.names = c(psffit$Data$mon.names,"chisq")
		bestpsf = getLDFitParm(psffit, Fit = psffit$LDFit)
		chisqredpsf = profitLikeModel(bestpsf, psffit$Data)$Monitor["chisq"]/ndata
	}

	#psffwhm = bestpsf["psf.moffat.fwhm"]
	# Make an image of the PSF which should be big enough for fitting purposes
	if("psf.moffat.fwhm" %in% colnames(psffit$LDFit$Posterior1))
		biggestpsf = getLDFitParm(psffit, which=which.max(psffit$LDFit$Posterior1[,"psf.moffat.fwhm"]))
	else # too much effort to find the biggest one in multi-component PSFs
		biggestpsf = bestpsf
	psfmodel = profitRemakeModellist(biggestpsf, Data = psffit$Data)$modellist$psf
	psfim = makeProfitPSFAutoSize(psfmodel,minpsfsum = 0.998, finesample=finesample,
		returnfine = TRUE)
	psfim = psfim/sum(psfim)
	psfdim = dim(psfim)
	if(mincutoutsize < max(psfdim)/finesample) mincutoutsize = round(max(psfdim)/finesample)
	if(FALSE)
	{
		ncomps = length(psfmodel$moffat$mag)
		psfmodel$moffat$xcen = rep(psfdim[1]*finesample/2,ncomps)
		psfmodel$moffat$ycen = rep(psfdim[2]*finesample/2,ncomps)
		psfim = profitMakeModel(psfmodel,finesample = finesample,dim = psfdim*finesample)$z
	}

	for(obj in as.character(fitobj))
	{
		#if(is.null(fits[[obj]]))
		{
			ra=segstats[obj,"RAcen"]
			dec=segstats[obj,"Deccen"]
			mag=segstats[obj,"mag"]
			re=segstats[obj,"R50"]/pixscale
			axrat=segstats[obj,"axrat"]
			ang=segstats[obj,"ang"]

			# ensure cutout size is divisible by 10
			cutoutsize = 2*5*ceiling(max(c(3,1.5)*process$single[[autocropband]]$stats[obj,
				c("R50","R90")]/pixscale)/5)
			if(cutoutsize < mincutoutsize) cutoutsize = mincutoutsize
			cutoutsize = c(cutoutsize,cutoutsize)

			cutouts = cutoutSourceWCS(image, header, ra, dec, box=cutoutsize,
				extras=list(segim=segim,sigma=sigma))
			xy = cutouts$xy

			init = list(x=unname(xy[,"x"]), y=unname(xy[,"y"]), size=re, mag=mag,
				slope=2, ang=ang, axrat=axrat, box=0)
			segobjs = as.integer(c(0, obj))

			galfit = profitFitGalaxyComponents(image=cutouts$maps$image$image,
				sigma=cutouts$maps$sigma, psfim=psfim, segim=cutouts$maps$segim,
				segobjs = segobjs, init=init, gain_eff=gain, skylevel = 0, maxcomp=maxcomp,
				chisqrtarg=chisqredpsf, chisqrminratio = chisqredpsf+0.02,
				nthreads=nthreads, maxiter=maxiter, domcmc = domcmc, galfit=fits[[obj]],
				finesample = finesample, dobrokenexp = dobrokenexp)
			fits[[obj]] = galfit
		}
	}
	return(fits)
}

# TODO: Finish this
makeModelSegmap <- function(process, psffit, galfits=NULL,
	finesample=1L,nthreads=1L,mincutoutsize=36, maxiter=2e3L,
	fitall=FALSE, maxcomp=2)
{
	bests = c()
	sds = c()
	galaxies = c()
	if(!is.null(galfits)) galaxies = names(galfits)
	pointsources = psffit$pointsources
	stopifnot(!(any(galaxies %in% pointsources) || any(pointsources %in% galaxies)))

	tocopy = c("tofit","tolog","intervals")
	objlists = list()
	for(item in tocopy) objlists[[item]] = psffit$Data[[item]]
	best = getLDFitParm(psffit)

	for(obj in galaxies)
	{
		if(sources[[obj]])
		{
			fitdata = psfits[[obj]]
		} else {
			fitdata = galfits[[obj]]$data$ser_frpa_3
			if(is.null(fitdata)) fitdata = galfits[[obj]]$data$ser_frpa_2
		}
		stopifnot(!is.null(fitdata))
		objlists[[obj]] = list()
		for(item in tocopy) objlists[[obj]][[item]] = fitdata[[obj]]$Data[[item]]
		best = getLDFitBest(fitdata,Fit = fitdata$LDFit)
		sds = rbind(sds,fitdata$LDFit$Summary1[1:fitdata$LDFit$Parameters,"SD"])
		modcols = c(paste0("sersic.",c(paste0(c("x","y"),"cen"),"mag")))#,"sky.bg")
		magmod = 0 #-2.5*log10(fits[[1]]$data$Data$gain$gain_eff)
		best[modcols] = best[modcols] - c(galfits[[obj]]$data$Data$xyoffsets,numeric(1)+magmod)
		bests = rbind(bests,best)
	}

	region = process$segim == 0
	pslists = list()
	for(obj in names(psfits))
	{
		pslists[[cobj]] = list()
		for(item in tocopy) pslists[[cobj]][[item]] = fits[[cobj]]$Data[[item]]
		fitdata =
		best = fits[[cobj]]$LDFit$Posterior1[which.max(fits[[cobj]]$LDFit$Monitor[,"LP"]),]
		lists[[cobj]]$modellist = profitRemakeModellist(
		)
	}
	multilists = profitCombineLists(lists)
	Data = profitSetupData(image = image, sigma=sigma, region=region,
		modellist = multilists$modellist, tofit = multilists$tofit,
		tolog = multilists$tolog, intervals = multilists$intervals,
		priors = multilists$priors
	)

	if(fitall)
	{
		psfpad = (dim(psfim)/finesample-1)/2
		cropmax = dimimg - psfpad
		cropmax = cropmax - ((cropmax-psfpad) %% 2)
		cropx = (1+psfpad[1]):cropmax[1]
		cropy = (1+psfpad[2]):cropmax[2]

		psflists = makeProfitLists(xcen=0,ycen=0,mag=0,	size=psffwhm,
			slope = unname(bestpsf["psf.moffat.con"]),
			ang = unname(bestpsf["psf.moffat.ang"]),
			axrat = unname(bestpsf["psf.moffat.axrat"]),
			ispsf = TRUE, profile="moffat",
			sizesd = psffit$Fit$Summary1["psf.moffat.fwhm","SD"],
			slopesd = psffit$Fit$Summary1["psf.moffat.con","SD"],
			angsd = psffit$Fit$Summary1["psf.moffat.ang","SD"],
			axratsd = psffit$Fit$Summary1["psf.moffat.axrat","SD"],
			boxsd = psffit$Fit$Summary1["psf.moffat.box","SD"],
			slopemin = 1.01, slopemax = Inf,
			fitsize = TRUE, fitslope=TRUE, fitaxrat = TRUE, fitang=TRUE,
			fitbox = TRUE, fitmag = FALSE,
			unlog=TRUE)

		psbest = list()
		pssds = list()
		bestvars = c("xcen","ycen","mag")
		for(var in bestvars)
		{
			varf = paste0("pointsource.",var)
			psbest[[var]] = unname(bestpsf[which(unlist(startsWith(names(bestpsf),varf)))])
			psrows = which(unlist(startsWith(rownames(psffit$Fit$Summary1),varf)))
			pssds[[var]] = unname(psffit$Fit$Summary1[psrows,"SD"])
		}

		probstars = as.numeric(psffit$starids)
		incrop = psbest[["xcen"]] > psfpad[1] & psbest[["xcen"]] < cropmax[1] &
				psbest[["ycen"]] > psfpad[2] & psbest[["ycen"]] < cropmax[2]
		for(var in bestvars)
		{
			psbest[[var]] = psbest[[var]][incrop]
			pssds[[var]] = pssds[[var]][incrop]
		}
		probstars = probstars[incrop]
		allt = !logical(length(probstars))

		segim = process$segexp$segim[cropx,cropy]
		posmask = segim <= 0
		for(id in c(galaxies,psffit$starids)) posmask = posmask | (segim == id)
		# For simultaneous fit
		mask = process$mask[cropx,cropy] | !posmask
		image = image[cropx,cropy]
		sigma = sigma[cropx,cropy]
		magimage(image*!mask)

		lists = combineProfitLists(list(psf=psflists, pointsource=pslists, sersic=sersiclists))

		Data = profitSetupData(image=image,sigma = sigma, mask=mask, psf=psfim,
			modellist = lists$modellist, tofit = lists$tofit, tolog=lists$tolog, intervals = lists$intervals,
			priors=lists$priors, algo.func="CMA", like.func="norm", verbose=FALSE, omp_threads=1,
			nbenchmarkconv = 1L, benchmarkconvmethods = c("FFTconv","FFTWconv"))
		#profitLikeModel(Data$init, Data, makeplots = T, plotchisq = T)

		mlfit = mlFit(Data, maxiter=250, algo.func="CMA", cmathreads=nthreads, tolerance=0.1)
		Data$algo.func="LD"
		Data$omp_threads = nthreads
		LDFit = convergeFit(Data, initalg="HARM", initspecs = list(alpha.star=0.234, B=NULL),
			inititer = maxiter/4L, finaliter=maxiter)

		rv$allfit = list(Data=Data,MLFit = mlfit, LDFit = LDFit)
	}

	return(rv)
}

allStarFitCompact <- function(proc, maps, bands = names(proc$single), maxcomp=1,
	pixscale=proc$multi$proc$pixscale, domcmc=TRUE, fits=list())
{
	stopifnot(all(bands %in% names(proc$single)))
	segim = proc$multi$proc$segim
	mask = proc$multi$mask
	for(band in bands)
	{
		bproc = proc$single[[band]]
		bmaps = maps[[band]]
		# Override sigma map if a new one was generated by processImage
		if(!is.null(bproc$err)) bmaps$err = bproc$err
		gain = bmaps$gain_eff
		sources = bproc$sources$compact
		fits[[band]] = fitPointSources(bmaps$img, bmaps$err, header = bmaps$hdr,
				mask=mask, segim=segim, segstats = bproc$stats, fitobj = sources,
				gain=gain, maxcomp=maxcomp, domcmc=domcmc, fits=fits[[band]],
				pixscale=pixscale)
		#plotSourceFits(psfits, procmaps, band)
	}
	return(fits)
}

allStarFitPSF <- function(procmaps, psfits, bands = names(procmaps$single), maxcomp=1,
	psselect="auto", fits = list())
{
	stopifnot(all(bands %in% names(psfits)))
	stopifnot(all(bands %in% names(procmaps$maps)))
	stopifnot(all(bands %in% names(procmaps$proc$single)))
	segim = procmaps$proc$multi$proc$segim
	mask = procmaps$proc$multi$mask
	#skypix = procmaps$proc$multi$proc$objects_redo==0
	skypix = procmaps$proc$multi$proc$segim<=0
	for(band in bands)
	{
		bproc = procmaps$proc$single[[band]]
		bmaps = procmaps$maps[[band]]
		bpsfits = psfits[[band]]
		gain = bmaps$gain_eff
		pscands = selectFittedPointsources(bpsfits, bproc$brightpsmag, bmaps$hdr,
			model="moffat1", select=psselect)
		fits[[band]] = fitPSF(bmaps$img, sigma=bproc$err, hdr=bmaps$hdr, mask=mask,
			segim=segim, psfits = bpsfits, pointsources=pscands$ids, means=pscands$means,
			sds=pscands$sds, init = pscands$means[pscands$best,], skypix=skypix, gain_eff=gain,
			skybg = 0, maxcomp=maxcomp, fits = fits[[band]])
		#finesample=1L,nthreads=1L, inititer=500, finaliter=2000)
	}
	return(fits)
}

# TODO: Finish this
allStarFitExtended <- function(procmaps, psffits,
	bands = names(procmaps$maps), maxcomp=2)
{
	stopifnot(all(bands %in% names(procmaps$maps)))
	stopifnot(all(bands %in% names(procmaps$proc$single)))
	stopifnot(all(bands %in% names(psffits)))
	fits = list()
	for(band in bands)
	{
		bpsffit = psffits[[band]]
		bmaps = procmaps$maps[[band]]
		gain = bmaps$gain_eff
		sources = procmaps$single[[band]]$sources$extended
		fits[[band]] = fitGalaxies(process = procmaps$proc, psffit = bpsffit,
			fitobj=sources, band=band, header=bmaps$hdr, maxcomp = maxcomp,
			gain = gain)
	}
	return(fits)
}

# TODO: Fix this
plotfitresiduals <- function(fits, multifit=NULL, fittype="data",
	usemultinames=!is.null(multifit), multirepad=FALSE,multiadd=c(0,0))
{
	hasmulti = !is.null(multifit)
	if(hasmulti)
	{
		multifit$Data$algo.func = ""
		multimod = profitLikeModelLDFit(multifit)$model$z
		if(multirepad)
		{
			olddim = dim(multimod)
			dimpsf = dim(multifit$Data$psf)
			psfmod = (dimpsf-1)/2
			newdim = dimpsf - 1 + dim(multimod) + multiadd
			print(newdim)
			newmod = matrix(0,newdim[1],newdim[2])
			newmod[psfmod[1] + 1:olddim[1], psfmod[2] + 1:olddim[2]] =
				multimod
			multimod = newmod
		}
	}
	if(hasmulti && usemultinames) fitnames = multifit$starids
	else fitnames = names(fits)
	for(fitname in fitnames)
	{
		fit = fits[[fitname]][[fittype]]
		lp = profitLikeModelLDFit(fit,makeplots=T,plotchisq=T)$LP
		readkey(paste0("Single fit for obj=",fitname,"; LP=",lp,". "))
		if(hasmulti)
		{
			xybl = 1-fit$Data$xyoffsets
			xytr = xybl + dim(fit$Data$image) - 1
			lp = profitLikeModelLDFit(fit,makeplots=T,plotchisq=T,
				model=list(z=multimod[xybl[1]:xytr[1],xybl[2]:xytr[2]]))$LP
			readkey(paste0("Simultaneous fit for obj=",fitname,"; LP=",lp,". "))
		}
	}
}

plotSourceFits <- function(bpsfits, procmaps, band, fullfits=list(NULL,NULL),
	profile="moffat", maxcomp=1, plotlims = list(
		mag=c(-Inf,20), size = c(-Inf,1), axrat = c(-0.1,0)))
{
	fits = bpsfits
	stopifnot(length(fullfits)==2)
	bmaps = procmaps$maps[[band]]
	bproc = procmaps$proc$single[[band]]
	fittypes = list(data=FALSE,mock=TRUE)
	fitname = paste0(profile,maxcomp)
	varnames = profitGetProfileParamNames(profile)
	magname = paste0(profile,".mag")
	xynames = paste0(profile,".",paste0(c("x","y"),"cen"))
	sizename = paste0(profile,".",varnames$size)
	for(fittype in names(fittypes)[1])
	{
		posts = list()
		for(ifit in 1:length(fits))
		{
			x = fits[[ifit]]
			fit = x[[fittype]][[fitname]]
			gaine = fit$Data$gain
			if(!is.null(gaine))
			{
				gainmag = -2.5*log10(gaine)
				fit$LDFit$Posterior1[,magname] = fit$LDFit$Posterior1[,magname]-gainmag
			}
			rv = list()
			within = TRUE
			for(name in names(plotlims))
			{
				if(name %in% names(varnames)) varname = varnames[[name]]
				else varname = name
				lims = plotlims[[name]]
				y = fit$LDFit$Posterior1[,paste0(profile,".",varname)]
				if(!any(y >= lims[1] & y <= lims[2])) within=FALSE
			}
			post=fit$LDFit$Posterior1
			radec = xy2radec(post[1,xynames[1]],post[1,xynames[2]], header = fit$Data$header)
			xyoff = radec2xy(radec[,"RA"],radec[,"Dec"],header = bmaps$hdr)-post[1,xynames]
			post[,xynames] = post[,xynames] + rep(xyoff,each=dim(post)[1])

			if(within) rv = list(post=post, mon=fit$LDFit$Monitor)
			posts[[ifit]] = rv
		}
		names(posts) = names(fits)
		for(post in names(posts))
		{
			if(length(posts[[post]]) == 0) posts[[post]] = NULL
		}
		segids = as.numeric(names(posts))

		ranges = list()
		plotdata = list()
		# https://www.r-bloggers.com/the-paul-tol-21-color-salute/
		cmap = c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
		ncmap = length(cmap)
		if(ncmap < length(posts))
		{
			cmap = colorRampPalette(cmap)(length(posts))
		}
		cmap = sample(cmap)
		objn = 1
		for(i in 1:length(posts))
		{
			if(!is.null(posts[[i]]))
			{
				post = posts[[i]]$post
				ncols = dim(post)[2]
				for(j in 1:ncols)
				{
					ranges[[j]] = range(c(post[,j],ranges[j]))
				}
				post[,paste0(profile,".ang")] = post[,paste0(profile,".ang")] %% 180
				plotdata[[paste0("obj",segids[i])]] = list(consd=0.02,
					data = post[,3:ncols], pointcol = cmap[objn], ptype=20, plotsd=FALSE
				)
				objn = objn+1
			}
		}
		print(ranges)

		best = c()
		for(post in posts) best = rbind(best,post$post[which.max(post$mon[,"LP"]),])
		colnames(best) = colnames(posts[[1]]$post)
		tolog = unlist(fits[[1]]$data[[fitname]]$Data$tolog$moffat)
		varnames = unlist(lapply(strsplit(colnames(best),paste0(profile,".")),function(x) {if(length(x)>1) return(x[[2]])}))

		print(best[1:2,])

		yeqx = c(-1,1)*.Machine$integer.max

		par(mfrow=c(2,2),mar=c(2,2,2,2),cex=par("cex"),lwd=par("lwd"))
		seg = bproc$stats
		idxs = match(segids,seg[,"segID"])
		seg = seg[idxs,]
		seg[,"mag"] = -2.5*log10(seg[,"flux"])
		magplot(seg[,"mag"], best[,magname])
		lines(yeqx,yeqx)
		seg[,"fwhm"] = sqrt(seg[,"N50"]/pi)
		if(tolog["fwhm"])
		{
			seg[,"fwhm"] = log10(seg[,"fwhm"])
			fwhmmod = best[,paste0(profile,".axrat")]
			if(tolog["axrat"]) fwhmmod = best[,sizename] + 0.5*fwhmmod
			else fwhmmod = best[,sizename]*sqrt(fwhmmod)
		}
		magplot(seg[,"fwhm"],fwhmmod)
		lines(yeqx,yeqx)
		seg[,"fwhm"] = seg[,"fwhm"] - 0.5*log10(seg[,"axrat"])
		for(var in c("axrat","con")) if(tolog[var]) seg[,var] = log10(seg[,var])
		seg[,"con"] = 1+seg[,"con"]
		magplot(seg[,"axrat"],best[,paste0(profile,".axrat")])
		lines(yeqx,yeqx)
		magplot(seg[,"ang"],best[,paste0(profile,".ang")] %% 180)
		lines(yeqx,yeqx)

		for(i in 1:length(posts))
		{
			data = seg[i,varnames]
			if("sky.bg" %in% colnames(best)) data = cbind(data,0)
			names(data) = colnames(best)
			# Why is data a list in the first place?! who knows
			stopifnot(all(is.finite(unlist(data))))
			plotdata[[paste0("obj_",i)]] = list(sample=4,samptype="nth",
				data = data[,3:length(data)], pointcol = cmap[i],pty=0,plotsd=FALSE)
		}
		for(ffi in 1:length(fullfits))
		{
			fullfit = fullfits[[ffi]]
			if(!is.null(fullfit))
			{
				gainmag = -2.5*log10(fullfit$Data$gain)
				data = fullfit$LDFit$Posterior1[1:fullfit$LDFit$Iterations,]
				psfpar = which(startsWith(colnames(data),"psf.") & !endsWith(colnames(data),"box"))
				psmag = unlist(lapply(strsplit(colnames(data),split = ".mag"), length)) == 2
				data = cbind(-gainmag -2.5*log10(rowSums(10^(data[,which(
					startsWith(colnames(data),"pointsource.") & psmag)]/-2.5))),data)
				data = data[,c(1,psfpar+1,which(colnames(data) == "sky.bg"))]
				colnames(data) = colnames(best)[3:length(colnames(best))]
				dark = (ffi-0.5)/length(fullfits)/2
				plotdata[[paste0("fullfit_",ffi)]] = list(data = data, ptype=20, consd=0.02,
					pointcol=rgb(dark,dark,dark))
			}
		}
	magtriplus(plotdata,samptype="nth",samples=5,sampallcon=TRUE)
	}
}
