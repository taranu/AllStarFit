profitMakeModelFine <- function(parm, Data,
	finesample=1L, nopsf=TRUE, openclenv=NULL)
{
	bestlist = profitRemakeModellist(parm, Data=Data)$modellist
	if(nopsf && !is.null(bestlist$psf)) bestlist$psf = NULL
	image = profitMakeModel(bestlist, dim=dim(Data$image),
		finesample = finesample, returnfine = TRUE, openclenv=openclenv)$z
	return(list(image=image,modellist=bestlist))
}

makeQuantileMaps <- function(image, modelflux=sum(image),
	quantiles=(1:9)/10, downsample = 1L, psf=NULL, plot=FALSE)
{
	stopifnot(identical(TRUE,is.sorted(quantiles)))
	dimimg = dim(image)
	sumimage = sum(image)
	sorted = sort(image, decreasing=T, index.return=T)
	cumsummed = cumsum(image[sorted$ix])
	nquantiles = length(quantiles)
	rval = list()

	quant = image
	quantic = length(cumsummed)
	stopifnot(quantic == prod(dimimg))

	quantord = rev(1:nquantiles)
	for(quanti in quantord)
	{
		quantile = quantiles[quanti]
		fluxtarget = quantile*modelflux
		if(fluxtarget < sumimage)
		{
			quanticnew = .Internal(findInterval(cumsummed[1:quantic], fluxtarget,
				FALSE, FALSE, FALSE))
			if(quanticnew < length(cumsummed)) quanticnew = quanticnew + 1
			stopifnot(quanticnew <= quantic)
			if(quantic > 1) quant[sorted$ix[(quanticnew + (quanticnew<length(cumsummed))):quantic]] = 0
			# Subtract the flux difference of the dimmest pixel to exactly reach the desired quantile
			diff = cumsummed[quanticnew]-fluxtarget
			stopifnot(diff >= 0)
			if(quanticnew == 1)
			{
				stopifnot(quant[sorted$ix[1]] >= fluxtarget)
			}
			quant[sorted$ix[quanticnew]] = quant[sorted$ix[quanticnew]] - diff
			if(quanticnew == 1) cumsummed[1] = fluxtarget
			if(plot) magimage(quant,stretch="log")
			quantic = quanticnew
			rval[[quanti]] = list(map=quant,quantile=quantile)
		} else {
			stopifnot(quanti>=2)
			quantiles = quantiles[1:(quanti-1)]
		}
	}
	names(rval) = paste0("p",sprintf("%.2f",100*quantiles))
	return(rval)
}

# Zeroes out negative pixels after FFT convolution, attempting to preserve
# total flux by keeping the sum of zeroed pixels close to zero
# Useful if you need a positive convolved image
fitfftconvpos <- function(image)
{
	if(any(image < 0))
	{
		sorted = sort.int(image, index.return = T)
		cumsums = cumsum(sorted$x)
		maxix = which(cumsums>0)[1]
		stopifnot(maxix >= 1)
		image[sorted$ix[1:maxix]] = 0
	}
	return(image)
}

testQuantileMaps <- function(gid, proj, kpcinasec, survey="kids", snap="008192",
	kinpixrat = 5/2, finesample = 1L, meanvspaxels = 8, spaxscale=0.5)
{
	name = load("/home/taranu/raid/romulus/snaps/halos/cosmo25_008192/maps/kids/0061270/profit_2.R.dat")
	stopifnot(length(name) == 1)
	fit = eval(as.name(name[[1]]))
	fit = fit[[1]]$base[[1]]$data$ser_frpa_2
	fit$Data$convopt$fft$fftwplan = NULL
	Data = setupprofitdata(basedir = file.path("~","raid","romulus","snaps","halos",
			paste0("cosmo25_",snap),"maps"), gname="", imgdir=file.path(survey,gid), dirs=c(base=""),
			imgname = paste0("image_smoothed_",proj), samipath = file.path("sami",gid))
	kinpath = file.path(kinpngpath,paste0(filename[[1]],"_kin.png"))
	measureMapQuantiles()
}

profitRescaleSersic <- function(dimimg, pars, enlarge=1, rescale=1, logrescale=FALSE)
{
	stopifnot(all(enlarge>0,rescale>0))
	if(any(enlarge != 1))
	{
		xynames = paste0("sersic.",paste0(c("x","y"),"cen1"))
		newdimimg = ceiling(enlarge*dimimg)
		pars[xynames] = pars[xynames]-dimimg/2 + newdimimg/2
	}
	else
	{
		newdimimg = dimimg
	}
	if(any(rescale != 1))
	{
		res = which(unlist(lapply(as.list(names(pars)), function(x) { return(identical(substr(x,1,nchar(x)-1),"sersic.re")) })))
		if(logrescale) pars[res] = pars[res] + log10(rescale)
		else pars[res] = rescale*pars[res]
	}
	rv = list(dimimg=newdimimg, pars=pars)
}

measureMapQuantiles <- function(galdata, fitdata, kpcinasec,
	kinpixrat = 5/2, finesample = 1L, kinfine=3L,
	meanvspaxels = 8, spaxscale=0.5, fovradinpix=15,
	returnmaps=FALSE, fileprefix=NULL, filepostfix=NULL, minkincounts=0,
	modelpars = NULL, minfluxfrac=0.95, openclenv=NULL)
{
	stopifnot(meanvspaxels >= 1)
	stopifnot(all(minfluxfrac>0,minfluxfrac <= 1))

	kinmaps = list(
		l = galdata$lmap,
		v = galdata$vmap,
		vd = galdata$vdmap
	)
	dimimg = dim(fitdata$Data$image)
	kinbox = dimimg/kinpixrat
	kinpsf = galdata$kinpsf
	# TODO - read the kin psf finesample automatically
	minkinbox = 10*ceiling(dim(kinpsf)/10/3)
	ltminkin = which(kinbox < minkinbox)
	kinbox[ltminkin] = minkinbox[ltminkin]
	for(mapname in names(kinmaps))
	{
		if(any(dim(kinmaps[[mapname]])*kinpixrat > dimimg))
		{
			kinmaps[[mapname]] = magcutout(kinmaps[[mapname]],box=kinbox)$image
		}
	}

	kindim = dim(kinmaps$v)
	stopifnot(identical(dim(kinmaps$vd),kindim))
	stopifnot(meanvspaxels <= prod(kindim))
	kindimrescale = kindim*kinpixrat*finesample

	if(!is.null(openclenv) && identical(typeof(openclenv),"externalptr") &&
		 !identical(openclenv, new("externalptr"))) fitdata$Data$openclenv = openclenv
	else if(identical(fitdata$Data$openclenv, new("externalptr"))) fitdata$Data$openclenv = NULL

	if("LDfit" %in% names(fitdata)) names(fitdata)[which(names(fitdata) == "LDfit")] = "LDFit"
	if(is.demonoid(fitdata$LDFit)) best = getLDFitParm(fitdata = fitdata)
	else {
		fit = fitdata$LDFit
		if(is.null(fit))
		{
			stopifnot(!is.null(fitdata$MLFit))
			fit = fitdata$MLFit
		}
		best = fit$par
	}
	stopifnot(is.numeric(best))
	if(is.numeric(modelpars))
	{
		stopifnot(identical(length(modelpars),length(best)))
		best=modelpars
	}
	stopifnot(length(best) == length(fitdata$Data$init))
	if(is.null(names(best))) names(best) = names(fitdata$Data$init)
	stopifnot(is.character(names(best)))
	gain = 1
	if(!is.null(fitdata$Data$gain)) gain = fitdata$Data$gain
	fitdim = dim(fitdata$Data$image)
	whichmin = which(fitdim < kindim*kinpixrat)
	if(length(whichmin) > 0)
	{
		newdim = fitdim
		newdim[whichmin] = 10*ceiling(kindim[whichmin]*kinpixrat/10)
		cennames = paste0("sersic.",c("xcen1","ycen1"))
		best[cennames[whichmin]] = best[cennames[whichmin]] + (newdim-fitdim)/2
		fitdata$LDFit$Posterior1[,cennames[whichmin]] = fitdata$LDFit$Posterior1[,cennames[whichmin]] + (newdim-fitdim)/2
		fitdata$Data$image = magcutout(galdata$image,box=newdim)$image*fitdata$Data$gain
		fitdata$Data$sigma = magcutout(galdata$errors,box=newdim)$image*fitdata$Data$gain
		#convopt = profitBenchmarkConv(image = profitUpsample(fitdata$Data$image,fitdata$Data$finesample),
		#	psf = fitdata$Data$psf)
		#fitdata$Data$convopt = convopt
		newdata = profitSetupData(image = fitdata$Data$image, psf=fitdata$Data$psf,
			psffinesampled = TRUE, finesample = fitdata$Data$finesample,
			sigma = fitdata$Data$image, modellist = fitdata$Data$modellist,
			tofit = fitdata$Data$tofit, nbenchmark = 2L)
		fitdata$Data$convopt = newdata$convopt
		fitdata$Data$calcregion = newdata$calcregion
		fitdata$Data$region = newdata$region
	}
	# Reset the lower/upper limit on re because we'll finesample later
	fitdata$Data$intervals$sersic$re = rep(list(lim=c(0,Inf)),length(fitdata$Data$intervals$sersic$re))
	# Same for x/y cen
	for(cen in paste0(c("x","y"),"cen"))
	{
		fitdata$Data$intervals$sersic[[cen]] = rep(list(lim=c(-Inf,Inf)),length(fitdata$Data$intervals$sersic[[cen]]))
	}
	finemodel = profitMakeModelFine(best, fitdata$Data, openclenv = fitdata$Data$openclenv)
	modelflux = sum(10^(-0.4*finemodel$modellist$sersic$mag))/gain
	converged = FALSE
	iter = 0
	while(!converged && (iter <= 3) && all(dimimg*finesample < 1e4))
	{
		finemodel = profitMakeModelFine(best, fitdata$Data, finesample = finesample,
			openclenv = fitdata$Data$openclenv)
		finemodel$image = finemodel$image/gain
		converged = sum(finemodel$image)/modelflux > minfluxfrac
		if(!converged)
		{
			# Enlarge the image and shift the centers to match
			prevdimimg = dimimg
			rescale = profitRescaleSersic(enlarge = 2, dimimg = dimimg, pars = best)
			best = rescale$pars
			dimimg = rescale$dimimg
			fitdata$Data$image = matrix(0,dimimg[1],dimimg[2])
		}
		iter = iter + 1
		#dimfine = dim(finemodel$image)
		#stopifnot(all(dimfine >= kindimrescale))
	}
	converged = sum(finemodel$image)/modelflux > 0.7
	#finemodel$image = magcutout(finemodel$image, box=kindimrescale)$image
	if(!converged)
	{
		stopifnot(converged)
	}

	quantm = makeQuantileMaps(finemodel$image, modelflux = modelflux)

	kinfov = ceiling(profitMakeModel(modellist = list(sersic=list(
			xcen=kindim[1]/2,ycen=kindim[2]/2, mag=0,
			nser=0.5,	re=fovradinpix, box=0)),
		dim = kindim, remax = 1, openclenv = fitdata$Data$openclenv)$z)
	kincond = as.logical(kinfov) & ((kinmaps$l*galdata$kingain) >= minkincounts)
	samipsfavg = profitMakeModel(modellist = list(
		moffat=list(xcen=kindim[1]/2,ycen=kindim[2]/2, mag=0, con=2, fwhm=4, box=0)),
		dim = kindim, openclenv = fitdata$Data$openclenv)$z
	print(c(sumpsf=sum(samipsfavg),sumpsffov=sum(samipsfavg[kincond])))

	# Make a SAMI-scale galaxy model and convolve with the SAMI PSF
	kinbest = profitRescaleSersic(dimimg = dimimg, pars = best, rescale = 1/kinpixrat, enlarge = fitdim/dimimg)
	fitdata$Data$image = matrix(0, kinbest$dimimg[1], kinbest$dimimg[2])
	kinimgfine = profitMakeModelFine(kinbest$pars, fitdata$Data, finesample = kinfine*finesample,
		openclenv = fitdata$Data$openclenv)$image/gain
	kindimfine = kindim*kinfine+dim(kinpsf)-1
	stopifnot(all(dim(kinimgfine) >= kindimfine))
	convmethods = profitAvailableConvolvers()
	kinfinedown = profitDownsample(kinimgfine,finesample)
	if(prod(c(dim(kinfinedown),dim(kinpsf))) > 1e9)
	{
		whichbrute = which(convmethods == "brute")
		if(length(whichbrute) > 0) convmethods = convmethods[-whichbrute]
	}
	bench = profitBenchmark(kinfinedown, methods = convmethods,
		psf=kinpsf,nbench = 1L, benchconvolution = TRUE)
	convolver = profitBenchmarkResultBest(bench$result)$convolver

	quantmkin = makeQuantileMaps(kinimgfine, modelflux = modelflux)

	kinconv = profitDownsample(fitfftconvpos(profitConvolve(convolver,
		kinfinedown, kernel=kinpsf)),kinfine)
	if(any(dim(kinconv) > kindim)) kinconv = magcutout(kinconv,box=kindim)$image
	print(c(sumkinconv=sum(kinconv),sumkinconvfov=sum(kinconv[kincond]))/modelflux)

	censpax = sort(kinconv, decreasing=TRUE, index.return=TRUE)$ix[1:meanvspaxels]
	meanv = mean(kinconv[censpax]*kinmaps$v[censpax]/sum(kinconv[censpax]))
	stopifnot(is.finite(meanv))
	# Recenter
	kinmaps$v = kinmaps$v - meanv
	vabs = abs(kinmaps$v)
	vsq = (kinmaps$v)^2
	vdsq = kinmaps$vd^2
	sone = sqrt(vdsq + vsq)
	goodkin = kinmaps$vd > 0

	cens = kindim/2
	spax = list(
		x = 0.5 + rep(1:kindim[1],kindim[2]),
		y = 0.5 + rep(1:kindim[2],each=kindim[1])
	)
	spax$r = sqrt((spax$x-cens[1])^2 + (spax$y-cens[2])^2)*kpcinasec*spaxscale

	rval = list()
	kpcperpix = spaxscale/kinpixrat*kpcinasec

	for(quantn in names(quantm))
	{
		quantmap=quantm[[quantn]]$map
		quantmapkin=quantmkin[[quantn]]$map

		dimfine = dim(quantmap)
		edgy = any(quantmap[1,] > 0)
		if(!edgy) edgy = any(quantmap[dimfine[1],] > 0)
		if(!edgy) edgy = any(quantmap[,1] > 0)
		if(!edgy) edgy = any(quantmap[,dimfine[1]] > 0)

		npix = sum(quantmap>0)
		quantconv = profitDownsample(fitfftconvpos(profitConvolve(convolver,
			profitDownsample(quantmapkin,finesample), kernel=kinpsf)),kinfine)
		quantconvfov = magcutout(quantconv,box=dim(kincond))$image*kincond
		sumquantconv=sum(quantconv)
		sumquantconvfov=sum(quantconvfov)
		sumquantrat = sumquantconvfov/sumquantconv
		print(c(quantconvrat=sumquantconv/modelflux,
			quantconvfovrat=sumquantconvfov/modelflux,
			sumquantrat=sumquantrat))

		stats = profoundSegimStats(quantmap, segim = 1*(quantmap>0), mask = 1*(quantmap<=0))
		lfov = quantconvfov
		lr = lfov*spax$r
		totlfov = sum(lfov)
		sumvabs = sum((vabs*lfov)[goodkin])
		sumvd = sum((kinmaps$vd*lfov)[goodkin])

		rval[[quantn]] = list(
			sigma = sumvd/totlfov,
			sone = sum(sone*lfov)/totlfov,
			vdivsigma = sumvabs/sumvd,
			vdivsone = sumvabs/sum((sone*lfov)[goodkin]),
			lambda = sum((vabs*lr)[goodkin])/sum((sone*lr)[goodkin]),
			vmax = max(vabs[goodkin & (lr>0)]),
			j = sum(vabs*lr)/totlfov,
			lquantfov=totlfov,
			lquant=sum(quantmap),
			quantratio = sumquantrat,
			quantile=quantm[[quantn]]$quantile,
			r = sqrt(npix/pi)*kpcperpix/finesample,
			ell = 1-stats[,"axrat"],
			hitedgeflag=edgy
		)
		if(returnmaps)
		{
			rval[[quantn]]$map = quantmap
			rval[[quantn]]$mapconv = quantconv
			rval[[quantn]]$mapfov = kincond
		}

		#for(map in names(kinmaps))
		if(FALSE)
		{
			maprange=c(0,150)+(map=="v")*c(-300,100)
			img = (kinmaps[[map]] - maprange[1])/(maprange[2]-maprange[1])
			img[img < 0] = 0
			img[img > 1] = 1
			#if(map=="vd") img[bad] = lmap[bad]
			kinmaps[[map]] = magcutout(img,box=newdim)$image
		}
		if(FALSE)
		{
			filename = paste0(gprefix,postfix)
			len = nchar(filename)
			if(len > 5 && substr(filename,len-4,len) == ".fits") filename = substr(filename,1,len-5)
			png(filename = , width=dim(maps$r)[1], height = dim(maps$r)[2])
		}
		#par(mar=rep(0,4))
		#image.galaxy.kin.color(l = lmap, v = kinmaps$v, vd = kinmaps$vd,
		#	stretch=kinstretch, vpow = kinvpow)
	}
	return(rval)
}

# psf = profitMakeModel(list(moffat=list(xcen=40.5,ycen=40.5,mag=0,fwhm=3*4,con=4.5,axrat=1,ang=0,box=0)),dim = rep(81,2))
# writeFITSim(psf$z/sum(psf$z),file = "~/raid/groups/gas_n4_3/inis/kids_psf_g_fine3.fits")
# psf = profitMakeModel(list(moffat=list(xcen=40.5,ycen=40.5,mag=0,fwhm=3*3.25,con=4.5,axrat=1,ang=0,box=0)),dim = rep(81,2))
# writeFITSim(psf$z/sum(psf$z),file = "~/raid/groups/gas_n4_3/inis/kids_psf_r_fine3.fits")
# psf = profitMakeModel(list(moffat=list(xcen=40.5,ycen=40.5,mag=0,fwhm=3.75*4,con=4.5,axrat=1,ang=0,box=0)),dim = rep(81,2))
# writeFITSim(psf$z/sum(psf$z),file = "~/raid/groups/gas_n4_3/inis/sdss_psf_g_fine3.fits")
# psf = profitMakeModel(list(moffat=list(xcen=40.5,ycen=40.5,mag=0,fwhm=3.5*4,con=4.5,axrat=1,ang=0,box=0)),dim = rep(81,2))
# writeFITSim(psf$z/sum(psf$z),file = "~/raid/groups/gas_n4_3/inis/sdss_psf_r_fine3.fits")
