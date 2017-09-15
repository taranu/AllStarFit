profitMakeModelFine <- function(parm, Data,
	finesample=1L, nopsf=TRUE)
{
	bestlist = profitRemakeModellist(parm, Data=Data)$modellist
	if(nopsf && !is.null(bestlist$psf)) bestlist$psf = NULL
	image = profitMakeModel(bestlist, dim=dim(Data$image),
		finesample = finesample, returnfine = TRUE)$z
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
			quant[sorted$ix[(quanticnew+1):quantic]] = 0
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
	kinpixrat = 5/2, finesample = 6L, meanvspaxels = 8, spaxscale=0.5)
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

measureMapQuantiles <- function(galdata, fitdata, kpcinasec,
	kinpixrat = 5/2, finesample = 6L,	meanvspaxels = 8, spaxscale=0.5,
	returnmaps=FALSE, fileprefix=NULL, filepostfix=NULL, minkincounts=0)
{
	stopifnot(meanvspaxels >= 1)
	kindown = 5
	alldown = 3
	stopifnot(kindown*alldown == kinpixrat*finesample)

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

	if("LDfit" %in% names(fitdata)) names(fitdata)[which(names(fitdata) == "LDfit")] = "LDFit"
	if(is.demonoid(fitdata$LDFit)) best = getLDFitParm(fitdata = fitdata)
	else {
		if(is.null(fitdata$LDFit)) {
			stopifnot(!is.null(fitdata$MLFit))
			best = fitdata$MLFit$par
		}
		else best = fitdata$LDFit$par
	}
	stopifnot(is.numeric(best))
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
			tofit = fitdata$Data$tofit, nbenchmarkconv = 2L)
		fitdata$Data$convopt = newdata$convopt
		fitdata$Data$calcregion = newdata$calcregion
		fitdata$Data$region = newdata$region

	}
	finemodel = profitMakeModelFine(best, fitdata$Data, finesample = finesample)
	dimfine = dim(finemodel$image)
	stopifnot(all(dimfine >= kindimrescale))
	finemodel$image = magcutout(finemodel$image, box=kindimrescale)$image/gain

	modelflux = sum(10^(-0.4*finemodel$modellist$sersic$mag))/gain
	quantm = makeQuantileMaps(finemodel$image, modelflux = modelflux)

	samifov = ceiling(profitMakeModel(modellist = list(
		sersic=list(xcen=kindim[1]/2,ycen=kindim[2]/2, mag=0, nser=0.5, re=15, box=0)),
		dim = kindim, remax = 1)$z)
	samicond = as.logical(samifov) & ((kinmaps$l*galdata$kingain) >= minkincounts)
	samipsfavg = profitMakeModel(modellist = list(
		moffat=list(xcen=kindim[1]/2,ycen=kindim[2]/2, mag=0, con=2, fwhm=4, box=0)),
		dim = kindim)$z
	print(c(sumpsf=sum(samipsfavg),sumpsffov=sum(samipsfavg[samicond])))

	samifine = profitDownsample(finemodel$image,kindown)
	benchconv = profitBenchmarkConv(samifine, psf=kinpsf)

	samiconv = profitDownsample(fitfftconvpos(profitConvolvePSF(
		samifine, psf=kinpsf, options=benchconv)),alldown)
	print(c(sumsamiconv=sum(samiconv),sumsamiconvfov=sum(samiconv[samicond]))/modelflux)

	censpax = sort(samiconv, decreasing=TRUE, index.return=TRUE)$ix[1:meanvspaxels]
	meanv = mean(samiconv[censpax]*kinmaps$v[censpax]/sum(samiconv[censpax]))
	# Recenter
	kinmaps$v = kinmaps$v - meanv
	vabs = abs(kinmaps$v)
	vabssq = vabs^2
	vdsq = kinmaps$vd^2
	smap = sqrt(vdsq + vabssq)

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
		quantfine = profitDownsample(quantmap,kindown)
		npix = sum(quantmap>0)/kindown^2
		quantconv = profitDownsample(fitfftconvpos(profitConvolvePSF(
			quantfine, psf=kinpsf, options=benchconv)),alldown)
		quantconvfov = quantconv*samicond
		sumquantconv=sum(quantconv)
		sumquantconvfov=sum(quantconvfov)
		sumquantrat = sumquantconvfov/sumquantconv
		print(c(quantconvrat=sumquantconv/modelflux,
			quantconvfovrat=sumquantconvfov/modelflux,
			sumquantrat=sumquantrat))

		stats = profitSegimStats(quantmap, segim = quantmap>0, mask = quantmap<=0)
		vdivsig = vabs/kinmaps$vd
		vdivsig[!is.finite(vdivsig)] = 0
		lfov = quantconvfov
		lr = lfov*spax$r
		totlfov = sum(lfov)
		rval[[quantn]] = list(
			sigma = sum(kinmaps$vd*lfov)/totlfov,
			sone = sum(smap*lfov)/totlfov,
			vdivsigma = sum(vdivsig*lfov)/totlfov,
			vdivs = sum(vabs*lfov)/totlfov,
			lambda = sum(vabs*lr)/sum(smap*lr),
			vmax = max(vabs/samiconv),
			j = sum(vabs*spax$r*lr)/totlfov,
			lquantfov=totlfov,
			lquant=sum(quantmap),
			quantratio = sumquantrat,
			quantile=quantm[[quantn]]$quantile,
			r = sqrt(npix/pi)*kpcperpix,
			ell = 1-stats[,"axrat"]
		)
		if(returnmaps)
		{
			rval[[quantn]]$map = quantmap
			rval[[quantn]]$mapconv = quantconv
			rval[[quantn]]$mapfov = samicond
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
