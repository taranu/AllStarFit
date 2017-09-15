makeThickDiskZbins <- function(nzbins=10, zdisk=1, nzbinsmax = 2*nzbins)
{
	zbins = matrix(0,nzbinsmax,4)
	dzi = 1.0/nzbins
	shrunkdz = FALSE
	dzmin=0
	dzprev=0.1*zdisk
	linearsampling = FALSE
	lineardz = 0
	zmax = 10.0*zdisk;
	linswitch = 1-1.5/nzbins
	maxzi = 1-1/nzbins/5

	zi = 0
	zindex = 1
	while(zi < maxzi && dzmin < zmax)
	{
		stopifnot(zindex <= nzbinsmax)
		if(zi > linswitch && !shrunkdz)
		{
			dzi = dzi/2
			shrunkdz = TRUE
		}

		# second last bin, don't use atanh(1) or you'll get inf
		# switch to linear sampling
		if(!linearsampling && zi > linswitch)
		{
			linearsampling = TRUE
			lineardz = dzprev
			dzi = 0;
		}
		if(linearsampling)
		{
			dzmax = dzmin + lineardz;
			zi = (tanh(dzmax/zdisk)+tanh(dzmin/zdisk))/2
			sumzi = tanh(dzmax/zdisk)-tanh(dzmin/zdisk)
			z = (dzmax+dzmin)/2
		}
		else
		{
			dzmax = zdisk*atanh(zi+dzi)
			z = zdisk*atanh(zi+dzi/2)
			sumzi = dzi
		}
		zbins[zindex,] = c(z,dzmin,dzmax,sumzi)
		zi = zi + dzi
		dzprev = dzmax - dzmin
		dzmin = dzmax
		zindex = zindex + 1
	}
	zbins = zbins[1:(zindex-1),]
	colnames(zbins) = c("mid","min","max","mass")
	nzbins = dim(zbins)[1]
	zbins[nzbins,"mass"] = 1 - sum(zbins[1:(nzbins-1),"mass"])
	return(zbins)
}

makeProfitThickDiskLists <- function(lists, disklists,
	zh, nzbins = 10, tofitzh=TRUE, tologzh=TRUE, zhinterval=c(0,Inf),
	inheritcomp = 0)
{
	stopifnot(length(lists) >= 1)
	sublists = names(lists)
	stopifnot(identical(sublists,names(disklists)))
	stopifnot(inheritcomp >= 0)

	nzbinsmax = 2*nzbins

	profile = names(disklists[[1]])
	stopifnot(length(profile) == 1)
	profile = profile[1]

	thickdiskname = paste0(profile,".thick")

	inheritpars = c()
	if(!is.null(disklists$tofit))
	{
		diskfits = disklists$tofit[[profile]]
		for(i in 1:length(diskfits))
		{
			if(is.na(diskfits[i])) inheritpars = c(inheritpars,i)
		}
	}
	ninherit = length(inheritpars)

	# Setup the new component
	for(sublist in sublists)
	{
		subprofile = names(disklists[[sublist]])
		stopifnot(length(subprofile) == 1)
		subprofile = subprofile[1]
		stopifnot(subprofile == profile)
		stopifnot(all(unlist(lapply(disklists[[sublist]][[profile]],length)) == 1))
		if(sublist == "modellist") zhval = zh
		else if(sublist == "tofit") zhval = tofitzh
		else if(sublist == "tolog") zhval = tologzh
		else if(sublist == "intervals") zhval = list(lim=zhinterval)
		else stop(paste("Unknown list type",sublist))
		if(sublist == "modellist")
		{
			nzbinsval = nzbins
			nzbinsmaxval = nzbinsmax
			indexval = length(lists$modellist[[profile]]$xcen)+1
			inheritcompval = inheritcomp
			inheritparsval = inheritpars
		} else
		{
			if(sublist == "intervals") ival = list(lim=c(-Inf,Inf))
			else ival = FALSE
			nzbinsval = ival
			nzbinsmaxval = ival
			indexval = ival
			inheritcompval = ival
			inheritparsval = rep(ival,ninherit)
		}

		newpars = list(
			tdzh=zhval,
			tdnzbins = nzbinsval,
			tdnzbinsmax = nzbinsmaxval,
			tdindex = indexval
		)
		if(ninherit > 0)
		{
			newpars$tdinheritcomp = inheritcompval
			newpars$tdinheritpars = inheritparsval
		}

		diskpars = disklists[[sublist]][[profile]]
		pars = names(diskpars)

		lists[[sublist]][[thickdiskname]] = c(
			disklists[[sublist]][[profile]], newpars)
		if(is.null(lists[[sublist]][[profile]])) lists[[sublist]][[profile]] = list()
		else
		{
			stopifnot(all(unlist(lapply(lists[[sublist]][[profile]][pars],length)) ==
				length(lists[[sublist]][[profile]][[pars[1]]])))
		}
		for(n in pars)
		{
			if(sublist == "tofit" || sublist == "intervals") subvar = ival
			else subvar = diskpars[[n]]
			lists[[sublist]][[profile]][[n]] = c(lists[[sublist]][[profile]][[n]],
				rep(subvar,2*nzbinsmax))
		}
	}

	return(lists)
}

makeProfitThickDisk <- function(modellist, tofit=NULL, tolog=NULL, intervals=NULL)
{
	if(!any(c(is.null(tofit),is.null(tolog),is.null(intervals))))
	{
		modelvec = unlist(modellist)
		tofitvec = unlist(tofit)
		tologvec = unlist(tolog)
		stopifnot(all(length(modelvec) == c(length(tofitvec),length(tologvec))))
		whichfit = which(tofitvec)
		parm = modelvec[whichfit]
		whichlog = which(tologvec[whichfit])
		parm[whichlog] = log10(parm[whichlog])
		modellist = profitRemakeModellist(parm,modellist = modellist, tofit = tofit, tolog = tolog, intervals = intervals)$modellist
	}
	whichthick = names(modellist)[which(endsWith(names(modellist),"thick"))]
	for(thickdiskname in whichthick)
	{
		thickdisk = modellist[[thickdiskname]]
		profile = strsplit(thickdiskname,"\\.")
		stopifnot(length(profile) == 1)
		stopifnot(length(profile[[1]]) == 2)
		profile = profile[[1]][1]
		stopifnot(all((unlist(lapply(thickdisk,length)) == 1) | (names(thickdisk) == "tdinheritpars")))

		nzbinstot = 2*thickdisk$tdnzbinsmax
		if(!is.null(thickdisk$tdinheritpars))
		{
			stopifnot(thickdisk$tdinheritcomp >= 1)
			comp = thickdisk$tdinheritcomp
			parnames = names(modellist[[profile]])
			stopifnot(all(thickdisk$tdinheritpars <= length(parnames)))
			for(name in parnames[thickdisk$tdinheritpars])
			{
				stopifnot(length(modellist[[profile]][[name]]) >= comp)
				thickdisk[[name]] = modellist[[profile]][[name]][comp]
			}
		}

		# y vector of minor axis
		axminor = thickdisk$ang*pi/180
		axminor = c(cos(axminor),sin(axminor))
		# assume ellipticity -> sin(inclination); now axminor is dx,dy per unit z
		axminor = axminor*(1-thickdisk$axrat)

		profilenames = names(modellist[[profile]])

		values = list()
		if(thickdisk$axrat < 1)
		{
			zbins = makeThickDiskZbins(nzbins = thickdisk$tdnzbins,
				zdisk = thickdisk$tdzh, nzbinsmax = thickdisk$tdnzbinsmax)
			nzbins = dim(zbins)[1]
			values$mag = thickdisk$mag -2.5*log10(zbins[,"mass"]/2)
			values$mag = c(values$mag,values$mag)
			values$xcen = thickdisk$xcen + axminor[1]*c(zbins[,"mid"],-zbins[,"mid"])
			values$ycen = thickdisk$ycen + axminor[2]*c(zbins[,"mid"],-zbins[,"mid"])
			extrabins = thickdisk$tdnzbinsmax - nzbins
			if(extrabins > 0)
			{
				indices = nzbinstot:(nzbinstot-2*extrabins+1)
				values$mag[indices] = Inf
				# These values don't really matter
				values$xcen[indices] = thickdisk$xcen
				values$ycen[indices] = thickdisk$ycen
			}
		} else {
			values$mag = c(thickdisk$mag,rep(Inf, nzbinstot-1))
		}
		indices = thickdisk$tdindex + (1:nzbinstot) - 1
		valuenames = names(values)

		for(name in names(thickdisk))
		{
			isprofile = name %in% profilenames
			isthick = name %in% c("tdzh","tdnzbins","tdnzbinsmax","tdindex","tdinheritcomp","tdinheritpars")
			stopifnot(isprofile || isthick)
			stopifnot(!(isprofile && isthick))
			# stopifnot ... length(thickprofile[[profile]] >= indices
			if(isprofile)
			{
				if(name %in% valuenames) modellist[[profile]][[name]][indices] = values[[name]]
				else modellist[[profile]][[name]][indices] = rep(thickdisk[[name]],nzbinstot)
			}
		}
	}
	return(modellist)
}

setupThickLists <- function(gid, nzbins = 10,finesample=1L, disktype="exp", fitalg = "cha",band="r")
{
	fits = file.path(samimodeldir,gid,"fits","pro",band,fitalg,paste0("serq",disktype,".R.dat"))
	stopifnot(file.exists(fits))
	fits = load(fits)
	fits = eval(parse(text=fits[1]))
	fits$Data$convopt$fft$fftwplan = NULL
	fits$Data$Nnotfit = prod(dim(fits$Data$region)) - sum(fits$Data$region)
	lists = list(modellist = NULL, intervals = NULL, tolog=NULL, tofit=NULL)
	disklists = lists
	for(name in names(lists))
	{
		lists[[name]] = fits$Data[[name]]
		disklists[[name]] = list(sersic=lapply(lists[[name]]$sersic, function (x) { x = x[2]}))
		# Keep only the bulge
		lists[[name]]$sersic = lapply(lists[[name]]$sersic, function (x) { x = x[1]})
	}

	thickmodel = makeProfitThickDiskLists(lists, disklists,
		zh=0.1*disklists$modellist$sersic$re,nzbins = nzbins, inheritcomp = 1)
	thickmodelmod = makeProfitThickDisk(thickmodel$modellist)

	rv = list(fits=fits, lists = thickmodel, modellist = thickmodelmod)

	img = profitMakeModel(disklists$modellist,dim=c(200,200),finesample = finesample)$z
	thickmodelmod$sersic$mag[1] = Inf
	imgt = profitMakeModel(thickmodelmod,dim=c(200,200),finesample = finesample)$z-thickmodelmod$sky$bg
	rv$imgs = list(thin=img, thick=imgt)

	return(rv)
}

setupThickData <- function(input = setupThickLists())
{
	init = input$fits
	best = init$Fit$Posterior1[which.max(init$Fit$Monitor[,"LP"]),]

	Data = profitSetupData(image = init$Data$image, region = init$Data$region, sigma = init$Data$sigma,
		segim = init$Data$segim, mask = init$Data$mask, modellist = input$lists$modellist, tofit = input$lists$tofit,
		tolog = input$lists$tolog, intervals = input$lists$intervals,
		constraints = makeProfitThickDisk, psf = init$Data$psf, finesample = init$Data$finesample, psffinesampled = init$Data$finesample>1,
		algo.func = init$Data$algo.func, like.func = init$Data$like.func)
	# This is crucial to apply intervals on the thick disk parameters first,
	# before profitLikeModel's call to remakemodellist applies constraints
	formals(Data$constraints)$tofit = Data$tofit
	formals(Data$constraints)$tolog = Data$tolog
	formals(Data$constraints)$intervals = Data$intervals
	matched = which(names(best) %in% names(Data$init))
	Data$init[names(best)[matched]] = best[matched]
	Data$init['sersic.thick.re'] = best['sersic.re2']
	Data$init['sersic.thick.mag'] = best['sersic.mag2']
	if(!is.na(Data$init["sersic.thick.nser"]) && !is.na(best["sersic.nser2"]))
		Data$init["sersic.thick.nser"] = best["sersic.nser2"]
	#profitLikeModel(Data$init,Data=Data,makeplots = TRUE, plotchisq = TRUE)

	return(Data)
}

fitthick <- function(thickdata = setupThickData(), init=NULL, nthreads=1)
{
	tofitidxs = which(unlist(thickdata$tofit))
	intervals = unlist(thickdata$intervals)
	lims = list(
		lower = intervals[endsWith(names(intervals),"lim1")][tofitidxs],
		upper = intervals[endsWith(names(intervals),"lim2")][tofitidxs]
	)
	tolog = unlist(thickdata$tolog)[tofitidxs]
	for(limname in names(lims))
	{
		lims[[limname]][tolog] = log10(lims[[limname]][tolog])
		stopifnot(all(!is.nan(lims[[limname]])))
	}

	cmasigmas = numeric(length(lims$lower)) + 0.1
	thickdata$algo.func = "CMA"

	cmafit = NULL
	#cmafit = cmaeshpc(thickdata$init, profitLikeModel, Data=thickdata, control=list(maxit=1e3, fnscale=-1.0, sigma=cmasigmas,
	#	diag.sigma=TRUE, diag.eigen=TRUE, diag.pop=TRUE, diag.value=TRUE, lambda=6, maxwalltime=Inf, trace=TRUE, stopfitness = 0,
	#	stop.tolx=1e-2*cmasigmas, nthreads=2), lower = lims$lower, upper=lims$upper)

	thickdata$verbose = TRUE
	thickdata$algo.func = "LD"
	thickdata$omp_threads = nthreads

	if(!is.null(init) && (length(init) == length(thickdata$init))) thickdata$init = init

	bestLP = -Inf
  converged = FALSE
  iter = 1
  maxiter = 50
  niters = 0
  baseiters = 5e2
  while(!converged)
  {
	  LDFit=LaplacesDemon(profitLikeModel, Initial.Values=thickdata$init, Data=thickdata,
			Iterations=baseiters,	Algorithm='CHARM', Thinning=1,	Specs=list(alpha.star=0.44))
	  newLP=which.max(LDFit$Monitor[,"LP"])
	  bestval=LDFit$Posterior1[newLP,]
	  newLP=LDFit$Monitor[newLP,"LP"]
	  thickdata$init = bestval
	  converged = ((newLP - bestLP) < exp(1)) || (iter >= maxiter)
	  bestLP = max(newLP,bestLP)
	  niters = niters + baseiters
	  iter = iter + 1
  }

	LDFit = LaplacesDemon(profitLikeModel, Data=thickdata, Initial.Values = thickdata$init, Iterations=10*baseiters, Thinning = 1,
		Algorithm = "CHARM", Specs = list(alpha.star = 0.44))
	fits = list(cma=cmafit, LD = LDFit)
	return(fits)
}

compareThickProfitMagriteFits <- function(gid, platescale=0.2,plotfits=TRUE)
{
	rv = list()
	if(!plotfits)
	{
		parmfr = par("mfrow")
		parmar = par("mar")
		par(mfrow=c(5,5),mar=c(0,0,0,0))
		bestpars = list()
	}
	for(disktype in c("exp","ser"))
	{
		thicklists = setupThickLists(gid=gid,disktype=disktype,nzbins = 25)
		thicklists$fits$Data$priors <- function(x,y) { return(0)}
		bestthin = thicklists$fits$Fit$Posterior1[which.max(thicklists$fits$Fit$Monitor[,"LP"]),]
		model = paste0(disktype,".thin")
		Data = thicklists$fits$Data

		if(plotfits) rv[[model]] = profitLikeModel(bestthin,Data,makeplots = T, plotchisq = T,maxsigma = 10)
		else
		{
			bestmodel = profitRemakeModellist(bestthin,modellist = Data$modellist,	tofit=Data$tofit,
				tolog = Data$tolog, intervals = Data$intervals, constraints = Data$constraints)
			bestpars[[model]] = bestmodel$parm
			tounlog = which(unlist(Data$tolog)[which(unlist(Data$tofit))])
			bestpars[[model]][tounlog] = 10^bestpars[[model]][tounlog]
			rv[[model]] = profitMakeModel(bestmodel$modellist,
				psf = Data$psf, finesample = Data$finesample, dim = dim(Data$image))$z
			if(!is.null(Data$modellist$sky$bg)) rv[[model]] = rv[[model]] - Data$modellist$sky$bg
		}
		thickdata = setupThickData(thicklists)
		if(disktype == "exp")
			best = as.numeric(strsplit("9.9924e+01         9.9772e+01         1.9323e+01         5.2915e-01        -1.6731e-01         1.3647e+02        -1.4957e-01         1.4919e+01         1.7332e+00        -3.4818e-01         1.4850e+00",split="     ")[[1]])
		else
			best = as.numeric(strsplit("9.9890e+01         9.9731e+01         1.9173e+01         5.5573e-01        -1.0492e-01         1.3677e+02        -1.1156e-01         1.4988e+01         1.6981e+00        -1.4988e-01        -2.3015e-01         1.1714e+00",split="     ")[[1]])
		model = paste0(disktype,".thick")
		Data = thickdata
		if(plotfits) rv[[model]] = profitLikeModel(best,thickdata,makeplots = T, plotchisq = T, maxsigma=10)
		else
		{
			bestmodel = profitRemakeModellist(best,modellist = Data$modellist,	tofit=Data$tofit,
				tolog = Data$tolog, intervals = Data$intervals, constraints = Data$constraints)
			bestpars[[model]] = bestmodel$parm
			tounlog = which(unlist(Data$tolog)[which(unlist(Data$tofit))])
			bestpars[[model]][tounlog] = 10^bestpars[[model]][tounlog]
			rv[[model]] = profitMakeModel(bestmodel$modellist, psf = Data$psf, finesample = Data$finesample, dim = dim(Data$image))$z
			if(!is.null(Data$modellist$sky$bg)) rv[[model]] = rv[[model]] - Data$modellist$sky$bg
		}
	}
	bestmr = readFITS(file.path(samimodeldir,gid,"best","analysis","maps","kidsdf","image_smoothed_0_r_star.fits"))$imDat
	magritedisk = strsplit(read.table(file.path(samimodeldir,gid,"best","disk.ini"),stringsAsFactors = FALSE)[,1],split = "=")
	bestpars[["magrite"]] = c(as.numeric(unlist(lapply(magritedisk, function(x) { return(x[2])}))),
		cos(read.table(file.path(samimodeldir,gid,"best","analysis","projections.dat"))[[2]]))
	names(bestpars[["magrite"]]) = paste0("magrite.",c(unlist(lapply(magritedisk, function(x) { return(x[1])})),"axrat"))
	if(plotfits) profitLikeModel(best,thickdata,makeplots = T, plotchisq = T, maxsigma=10) #, cutmod = bestmr)
	else
	{
		rv[["magrite"]] = bestmr
		cmap = rev(colorRampPalette(brewer.pal(9,'Greys'))(200))

		img = Data$image
		if(!is.null(Data$modellist$sky$bg)) img = img - Data$modellist$sky$bg
		dimimg = dim(img)

		sig = Data$sigma
		reg = Data$region
		if(!is.null(reg))
		{
    	cutimg=img[reg]
    	cutsig=sig[reg]
  	}
		resids = list()
		LLs = list()
		LL = list()
		chisq = list()
		chisqLL = list()
		dofs = sum(reg)
		mnames = names(rv)
		n = length(mnames)
		for(model in mnames)
		{
			resids[[model]] = (img-rv[[model]])/sig
			LLs[[model]] = resids[[model]]
			LLs[[model]][!reg] = 0
			LLs[[model]] = dnorm(LLs[[model]],log = TRUE)
			LL[[model]] = sum(LLs[[model]][reg])
			chisq[[model]] = sum((resids[[model]][reg])^2)
			chisqLL[[model]] = dchisq(chisq[[model]],df=dofs,log=TRUE)
		}
		residrange = range(unlist(lapply(resids,function(x) { range(x[reg])})))
		LLsrange = range(unlist(lapply(LLs,function(x) { range(x[reg])})))
		zlim = max(abs(residrange))*c(-1,1)
		zlimLL = max(abs(LLsrange))*c(-1,1)
		textcol="darkred"
		for(j in 1:n)
		{
			for(i in 1:n)
			{
				if(i < j)
				{
					diffLL = LLs[[i]]-LLs[[j]]
					zlimLLi = max(abs(range(diffLL)))*c(-1,1)
					magimage(diffLL,stretch="linear",col=cmap,zlim=zlimLLi,magmap=FALSE,axes=FALSE)
					text(x = dimimg[1]-1,y = 1,labels = bquote(paste(Delta,LL[norm],.(sprintf("=%.3e",LL[[j]]-LL[[i]])))),adj = c(1,0),col = textcol, cex=1.5)
					text(x = 1,y = dimimg[2]-1,labels = bquote(paste(Delta,LL,.(sprintf(":[%.1f,%.1f]",zlimLLi[1],zlimLLi[2])))),adj = c(0,1),col = textcol, cex=1.5)
				}
				else if(i == j)
				{
					magimage(resids[[i]],stretch="linear",col=cmap,zlim=zlim,magmap=FALSE,axes=FALSE)
					text(x = dimimg[1]-1,y = 1,labels = bquote(paste(chi[red]^2,.("="),.(sprintf("%.2f",chisq[[j]]/dofs)))),adj = c(1,0),col = textcol, cex=1.5)
					text(x = 1,y = dimimg[2]-1,labels = bquote(paste(chi^2,.("LL="),.(sprintf("%.4e",chisqLL[[j]])))),adj = c(0,1),col = textcol, cex=1.5)
					if(endsWith(mnames[[i]],"magrite")) inclzd = bestpars[[mnames[i]]][paste0("magrite.",c("axrat","zdisk"))]
					else
					{
						if(endsWith(mnames[[i]],"thin")) inclzd = c(bestpars[[mnames[i]]]["sersic.axrat2"],0)
						else if(endsWith(mnames[[i]],"thick")) inclzd = bestpars[[mnames[i]]][paste0("sersic.thick.",c("axrat","tdzh"))]
						inclzd[2] = inclzd[2]*platescale
					}
					inclzd[1] = acos(inclzd[1])*180/pi
					text(x = 1,y = 1,labels = bquote(paste(.(sprintf("[i,z]=[%.1f,%.2f]",inclzd[1],inclzd[2])))),adj = c(0,0),col = textcol, cex=1.5)
				}
				else
				{
					magimage(rv[[i]]-rv[[j]],stretch="linear",col=cmap,zlim=2.4e-11*c(-1,1),magmap=FALSE,axes=FALSE)
				}
				if(i == 1) mtext(mnames[[j]],side=2,padj=2,col=textcol)
				if(j == 1 && i > 1) mtext(mnames[[i]],side=3,padj=2,col=textcol)
			}
		}

		par(mfrow=parmfr,mar=parmar)
	}
	return(rv)
}
