makeProfoundSegPlot <- function(image, seg, mask=NULL,
	maxlevel=NULL, minsig=1, tintmap=FALSE, err=NULL)
{
	# https://www.r-bloggers.com/the-paul-tol-21-color-salute/
	cmap = rev(brewer.pal(9,'RdYlBu'))
	#cmap = c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
	uids = sort(unique(as.numeric(seg$segim)))
 	cmap = c("#D0D0D0","#000000",colorRampPalette(cmap)(sum(uids>0)))
  segim = seg$segim
  segim[segim<=0] = -1
  for(i in 2:length(uids)) segim[segim==uids[i]] = i-1
  segim[segim<=0] = 0
  if(!is.null(mask))
  {
  	segim[mask!=0] = -1
  	if(tintmap) image[mask != 0] = 0
  }
  # Doesn't work well yet. TBD.
  if(tintmap)
  {
  	if(is.null(maxlevel)) maxlevel = max(image,na.rm=TRUE)
  	hsv = rgb2hsv(col2rgb(cmap))
  	segimh = matrix(0,dim(segim)[1],dim(segim)[2])
  	segims = segimh
  	for(id in 2:length(uids))
  	{
  		segeq = segim == id
  		segimh[segeq] = hsv[1,id+1]
  		segims[segeq] = hsv[2,id+1]
  	}
  	segimv = image
  	segimv[segimv > maxlevel] = maxlevel
  	segimv = segimv/err
  	segimv[segimv < 0] = 0
  	segimv = log10(segimv+1)
  	segimv = segimv/max(segimv)
  	segimrgb = col2rgb(hsv(segimh,segims,segimv))/255
  	dim(segimrgb) = c(3,dim(segim))
  	segimrgb = aperm(segimrgb,c(2,3,1))
  	magimageRGB(segimrgb,stretch="lin",magmap=FALSE)
  } else {
  	magimage(segim+1,col=cmap,stretch="lin",magmap=FALSE)
  }
}

proFoundCutouts <- function(maps, plot=FALSE, bands=names(maps))
{
	stopifnot(is.list(maps) && (length(maps) >= 1))
	proc = list()
	bands = names(maps)
	for(band in bands)
	{
		proc[[band]] = profoundProFound(maps[[band]]$img, mask=maps[[band]]$mask,
			skyRMS=maps[[band]]$err, platescale = kidsinfo$platescale,
			tolerance = 5, smooth=0, plot=plot)
	}
	return(proc)
}

proFoundMultiband <- function(maps, bands=names(maps),
	...)
{
	stopifnot(all(bands %in% names(maps)))
	nbands = length(bands)
	for(bandi in 1:nbands)
	{
		band = bands[bandi]
		img = maps[[band]]$img
		err = maps[[band]]$err
		mask = maps[[band]]$mask
		if(bandi == 1)
		{
			dimimg = dim(img)
			allimg = matrix(0,dimimg[1],dimimg[2])
			allvar = allimg
			allmask = allimg
		}
		stopifnot(identical(dimimg,dim(img)) && identical(dimimg,dim(err))
			&& identical(dimimg,dim(mask)))
		allimg = allimg + img*maps[[band]]$gain_eff
		allvar = allvar + (err*maps[[band]]$gain_eff)^2
		# Variance weighted, not error weighted
		allmask = allmask + mask
	}
	allimg = allimg / sqrt(allvar)
	allmask = 0+(allmask>0)
	proc = profoundProFound(image = allimg, mask=allmask, axes=FALSE, mar=rep(0,4), ...)
	return(list(proc=proc,img=allimg, mask=allmask))
}

# Select point sources that have been fit
selectFittedPointsources <- function(fits, brightmag, cutoutheader, sigmacut=5, minmag=Inf,
	model="moffat1", select="auto")
{
	sds = c()
	bests=c()
	for(obj in names(fits))
	{
		fit = fits[[obj]]$data
		stopifnot(model %in% names(fit))
		fit = fit[[model]]
		Data = fit$Data
		Fit = fit$LDFit
		best = Fit$Posterior1[which.max(Fit$Monitor[,"LP"]),]
		sds = rbind(sds,Fit$Summary1[1:Fit$Parameters,"SD"])
		modcols = c(paste0("moffat.",c(paste0(c("x","y"),"cen"),"mag")),"sky.bg")
		magmod = -2.5*log10(Data$gain)
		xy1 = best[paste0("moffat.",c(paste0(c("x","y"),"cen")))]
		xyo = xy2radec(xy1[1],xy1[2],Data$header)
		xyo = radec2xy(xyo[1],xyo[2],header = cutoutheader)
		best[modcols] = best[modcols] - c(xy1-xyo,rep(magmod,2))
		bests = rbind(bests,best)
	}
	rownames(bests) = NULL
	bright = (bests[,"moffat.mag"] < brightmag)
	# why did I bother with this? & (bests[,"moffat.axrat"] < 0)
	stopifnot(sum(bright) > 0)
	smallestid = which.min(bests[bright,"moffat.fwhm"])
	smallestid = which(cumsum(bright) == smallestid)[1]
	smallest = bests[smallestid,]
	smallsd = sds[smallestid,]
	smalllo = smallest-sigmacut*smallsd
	smallhi = smallest+sigmacut*smallsd
	bestlo = bests - sigmacut*sds
	besthi = bests + sigmacut*sds
	probps = which(bests[,"moffat.mag"] < minmag)
	if(select == "auto")
	{
		for(var in paste0("moffat.",c("fwhm","con","axrat")))
		{
			probps = probps & (bestlo[,var] < smallhi[var]) &
				(besthi[,var] > smalllo[var])
		}
	}

	psids = names(fits)[probps]
	rv = list(ids=psids, means = bests[probps,], sds = sds[probps,],
		best=smallestid)
	return(rv)
}

# Select sources brighter than a threshold magnitude and
# estimate whether they are extended or not
selectSources <- function(stats, minedgedist, dimimg, minmag=20,
	brightpsmag=NULL, refps=NULL, reffmargindex=0.1)
{
	stopifnot(((!is.null(refps)) + (!is.null(brightpsmag))) == 1)
	reff = stats[,"R50"]*sqrt(stats[,"axrat"])
	close = stats[,"xcen"] < (dimimg[1]-minedgedist) & stats[,"xcen"] > minedgedist &
		stats[,"ycen"] < (dimimg[2]-minedgedist) & stats[,"ycen"] > minedgedist
	estreffps = min(reff[stats[,"mag"] < brightpsmag & close],na.rm=T)
	good = (stats[,"mag"] < minmag) & close
	psrc = reff < (estreffps * 10^reffmargindex)
	esrc = !psrc & good
	psrc = psrc & good
	rval = list(compact = stats$segID[psrc],
		extended = stats$segID[esrc])
	return(rval)
}

processImagesMultiband <- function(maps, bands=names(maps), ...)
{
	multiproc = proFoundMultiband(maps, bands=bands, ...)
	return(list(multi=multiproc))
}

processImageSources <- function(maps, minpsmags, brightpsmags,
	mask, segim, bands=names(maps), fitbands=bands, pixscale=1,
	minedgedist=15, addobjerr = TRUE)
{
	errcond = (mask==0) & (segim!=0)
	rv = list()
	for(band in bands)
	{
		bmaps = maps[[band]]
		dimimg = dim(bmaps$img)
		brightpsmag = brightpsmags[[band]]
		gain = bmaps$gain_eff
		stats = profoundSegimStats(image=bmaps$img, segim=segim, skyRMS = bmaps$err,
			mask = mask, pixscale = pixscale, gain = gain, header = bmaps$hdr)
		sources = selectSources(stats, minedgedist, dimimg, minmag=minpsmags[[band]],
			brightpsmag = brightpsmag)
		rv[[band]]=list(stats=stats, sources=sources, brightpsmag = brightpsmag)
		if(addobjerr)
		{
			err = bmaps$img
			err[err<0] = 0
			err[mask!=0] = 0
			err[errcond] = sqrt((bmaps$err[errcond]*gain)^2+err[errcond]*gain)/gain
			err[!errcond] = bmaps$err[!errcond]
			rv[[band]]$err = err
		}
	}
	return(list(single=rv))
}
