# Shamelessly stolen from stackoverflow
.colSdColMeans <- function(x, na.rm=TRUE) {
  if (na.rm) {
    n <- colSums(!is.na(x)) # thanks @flodel
  } else {
    n <- nrow(x)
  }
	cm = colMeans(x, na.rm=na.rm)
  colVar <- colMeans(x*x, na.rm=na.rm) - cm^2
  return(list(mean=cm,sd=sqrt(colVar * n/(n-1))))
}

.colindex <- function(data, colname)
{
	return(which(colnames(data) == colname))
}

magriteMagtriPointcols <- function(statistic, basecol, minsat = 0, maxval = 1)
{
	basergb = col2rgb(basecol)
	basehsv = rgb2hsv(basergb)
	stopifnot(basehsv["s",] > minsat)
	stopifnot(basehsv["v",] < maxval)
	cols = range(statistic)
	value = log10(1+9*(statistic-cols[1])/diff(cols))
	ncols = length(value)

	cols = hsv(basehsv["h",]+numeric(ncols),
		basehsv["s",] - value*(basehsv["s",]-minsat),
		basehsv["v",] - value*(basehsv["v",]-maxval)
		)

	return(cols)
}

# Chains is a list with the following elements:
# 	data - the data in a standard matrix or table format (rows, cols)
# 	samples - (optional) The number of data points to plot (all others ignored).
# 		Default: dim(data)[1], i.e. all points
# 	samptype - (optional) A sampling method.
# 		Default: "end"; starting from the end.
# 		Options: "ran" for random, "begin" from the start
# 	pointcols - (optional) colours for the points
# 	pointcol - (optional) colour for all points, if pointcols not specified
# 		Defaults to the value of defaultcols at the chain's index
# 	linecol - (optional) colour for lines. Defaults to pointcol if not specified
# 	concol - (optional) colour for contours. Defaults to linecol if not specified
# 	plotsd - (optional) plot lines at +/- 1 sigma
#
# 	column names will be read from the list of data.
# 	The default is to plot all columns in at least one data set
magtriplus=function(chains, plotlegend=FALSE, samples=NULL, samptype='end',
	sampallcon = FALSE, grid=FALSE, tick=FALSE, plotcols = NULL,  pdflog=FALSE,
	defaultcols = c("black","red","blue","green")
)
{
	# Current cex value
	ccex = par("cex")
	# 1D 1-sigma, 2D 1,2,3 sigma
	ltys = c(2,5,4,3)
	mltys = c(6,1)
	ptys = c(4,3)

  if(!is.list(chains))
  {
  	stopifnot(is.numeric(chains))
  	chains = list(list(data=chains))
  }
  nchains = length(chains)
  getplotcols = is.null(plotcols)
  nvecchains = 0
  allsamples = samples
  allsamptype = samptype
  for(i in 1:nchains)
  {
  	if(is.null(dim(chains[[i]]$data))) chains[[i]]$data = t(chains[[i]]$data)
  	chains[[i]]$data=as.data.frame(chains[[i]]$data)
  	chains[[i]]$dim = dim(chains[[i]]$data)
	  Nsamp=chains[[i]]$dim[1]
	  if(i == 1)
	  {
	  	Npar=chains[[i]]$dim[2]
	  	stopifnot(Npar >=1)
	  }
	  else stopifnot(chains[[1]]$dim[2] == Npar)
	  samples = chains[[i]]$samples
	  if(is.null(samples)) samples = allsamples
	  if(is.null(samples) || (samples > Nsamp) || (samples < 1)) samples = Nsamp
	  samptype = chains[[i]]$samptype
	  if(is.null(samptype)) samptype = allsamptype
	  if(is.null(samptype)) samptype = "end"
 	  if(samptype=='begin')
	  {
	    samples = 1:Nsamp
	  }
		else if(samptype=='ran') samples=sample(Nsamp,samples)
	  else if(samptype=='nth') samples=seq(1,Nsamp,by = samples)
		else samples=(chains[[i]]$dim[1] - Nsamp + 1):chains[[i]]$dim[1]

		chains[[i]]$samples = samples
		meansds = .colSdColMeans(chains[[i]]$data)
		chains[[i]]$means = meansds$mean
		chains[[i]]$sds = meansds$sd
		chains[[i]]$isvec = Nsamp == 1
		nvecchains = nvecchains + chains[[i]]$isvec

	  if(getplotcols)
  	{
  		stopifnot(!is.null(colnames(chains[[i]]$data)))
  		plotcols = union(plotcols,colnames(chains[[i]]$data))
	  }
  	if(is.null(chains[[i]]$pointcol))
  	{
  		if(i > length(defaultcols)) chains[[i]]$pointcol =
  			colors(distinct=TRUE)[floor(1+runif(1)*length(colors(distinct=TRUE)))]
  		else chains[[i]]$pointcol = defaultcols[i]
  	}
		if(is.null(chains[[i]]$pointcols)) chains[[i]]$pointcols = chains[[i]]$pointcol
		if(is.null(chains[[i]]$pointtypes)) chains[[i]]$pointtypes = 16
	  if(is.null(chains[[i]]$linecol)) chains[[i]]$linecol = chains[[i]]$pointcol
	  if(is.null(chains[[i]]$concol)) chains[[i]]$concol = chains[[i]]$linecol
		if(is.null(chains[[i]]$consd)) chains[[i]]$consd = NULL
		if(is.null(chains[[i]]$ltys)) chains[[i]]$ltys = ltys
		else stopifnot(length(chains[[i]]$ltys) == 4)
		if(is.null(chains[[i]][["lty"]])) chains[[i]][["lty"]] = mltys[1+chains[[i]]$isvec]
		else stopifnot(is.numeric(chains[[i]][["lty"]]) && length(chains[[i]][["lty"]]) == 1)
		if(is.null(chains[[i]]$pty)) chains[[i]]$pty = ptys[1+chains[[i]]$isvec]
		else stopifnot(length(chains[[i]]$pty) == 1)
  }
  stopifnot((nchains-nvecchains) > 0)
  ncols = length(plotcols)
  legheight = 0.4+0.2*(nchains-nvecchains+ceiling(nvecchains/ncols))
  layout(t(matrix(c(numeric(ncols*plotlegend)+(ncols^2+1),1:ncols^2),ncols,ncols+plotlegend))[(ncols+plotlegend):1,],
  	heights=c(rep.int(1,ncols),rep(legheight,1*plotlegend)))
  par(oma=c(4*(1-plotlegend),4,4,4))
  ranges = list()
  pdftext="PDF"
  if(pdflog) pdftext = "log10(PDF)"

  for(col in plotcols)
  {
  	rangei=NULL
  	for(chain in 1:nchains)
  	{
  		i = .colindex(chains[[chain]]$data, col)
  		rangei = range(c(rangei,chains[[chain]]$data[chains[[chain]]$samples,i]))
  	}
  	if(rangei[1] == rangei[2]) rangei = c(-0.01,0.01) + rangei[1]
    ranges[[col]] = rangei
  }
  plotcolsi = 1:ncols
  for(i in plotcolsi)
  {
  	row = plotcols[i]
    yrange = ranges[[row]]
    diffr = yrange[2] - yrange[1]
    yrange = yrange + 0.01*diffr*c(-1,1)

    for(j in plotcolsi)
    {
      par(mar=c(0,0,0,0))
    	col = plotcols[j]
      xrange = ranges[[col]]
      diffr = xrange[2] - xrange[1]
	    xrange = xrange + 0.01*diffr*c(-1,1)

      if(i == j)
      {
      	ylim = NULL
      	ytemp = list()
      	for(chaini in 1:nchains)
      	{
      		chain = chains[[chaini]]
      		if(!chain$isvec && row %in% colnames(chain$data))
      		{
		        xtemp = chain$data[,row]
		        if(!sampallcon) xtemp = xtemp[chain$samples]
		        ytemp[[chaini]] = density(xtemp)
		        ylim = range(c(ylim,ytemp[[chaini]]$y))
      		}
      	}
      	for(chaini in 1:nchains)
      	{
      		chain = chains[[chaini]]
	        if(chain$isvec) abline(v=chain$data[,row],lty=1,col=chain$linecol)
	        else if(!is.null(ytemp[[chaini]]))
	        {
	        	if(chaini == 1)
	        	{
	        		if(pdflog)
	        			plot(ytemp[[chaini]]$x,log10(ytemp[[chaini]]$y),axes=FALSE,main='',xlim=xrange,
									ylim=log10(c(max(ylim[1],ylim[2]/1e6),ylim[2])),yaxs='i',col=chain$linecol,type="l")
	        		else plot(ytemp[[chaini]],axes=FALSE,main='',xlim=xrange,ylim=ylim,yaxs='i',col=chain$linecol)
	        	}
	        	else
	        	{
	        		if(pdflog) lines(ytemp[[chaini]]$x,log10(ytemp[[chaini]]$y),col=chain$linecol)
	        		else lines(ytemp[[chaini]],col=chain$linecol)
	        	}
	        	meani = chain$means[row]
	        	sdi= chain$sds[row]
            abline(v=meani,lty=chain$lty,col=chain$linecol)
            if(is.null(chain$plotsd) || isTRUE(chain$plotsd))
            {
							abline(v=meani-sdi,lty=chain$ltys[1],col=chain$linecol)
							abline(v=meani+sdi,lty=chain$ltys[1],col=chain$linecol)
            }
	        }
      	}
      	magaxis(1,grid=grid, grid.col = 'lightgrey',labels=FALSE,tick=tick)
        box()
      }else{
        if(i<j){
          plot.new()
          plot.window(xlim=xrange,ylim=yrange)
          magaxis(1:2,grid=grid, grid.col = 'lightgrey',labels=FALSE,tick=tick)
          for(chaini in 1:nchains)
          {
          	chain = chains[[chaini]]
          	xtemp = chain$data[,col]
          	ytemp = chain$data[,row]
          	if(!sampallcon)
          	{
          		xtemp = xtemp[chain$samples]
          		ytemp = ytemp[chain$samples]
          	}
          	if(!chain$isvec && !(all(xtemp == xtemp[1]) || all(ytemp == ytemp[1])))
          	{
          		if(is.null(chain$consd)) magcon(xtemp,ytemp, dobar=FALSE, doim=FALSE,add=TRUE,
								lty=chain$ltys[2:4], xlim=xrange,ylim=yrange, col=chain$concol, ngrid=1e3)
          		else magcon(xtemp,ytemp, chain$consd*c(diff(xrange),diff(yrange)), dobar=FALSE, doim=FALSE,add=TRUE,
								lty=chain$ltys[2:4], xlim=xrange,ylim=yrange, col=chain$concol, ngrid=1e3)
          	}
	          points(chain$means[col], chain$means[row],col=chain$pointcol,pch=chain$pty,cex=2*ccex)
	          if(chaini == 1)
	          {
	          	box()
	          	if(nchains > 1) par(new=TRUE)
	          }
	          meani = chain$means[col]
	          #abline(v=meani,lty=chain$lty,col=chain$linecol)

	          if(FALSE && !chain$isvec && !(!is.null(chain$plotsd) && isTRUE(chain$plotsd)))
	          {
	          	sdi = chain$sds[col]
	          	if(sdi > 0)
	          	{
								abline(v=meani-sdi,lty=chain$ltys[1],col=chain$linecol)
								abline(v=meani+sdi,lty=chain$ltys[1],col=chain$linecol)
	          	}
	          }
          }
        }
      	else
        {
          plot.new()
          plot.window(xlim=xrange,ylim=yrange)
          magaxis(1:2,grid=grid, grid.col = 'lightgrey',labels=FALSE,tick=tick)
          for(chain in chains)
          {
          	if(is.null(chain$pointcols)) pch='.'
          	else pch=chain$pointtypes
          	pcex = ccex/(1+0.5*!is.null(chain$pointcols))
            pointcol = chain$pointcols
          	if(length(pointcol) >= max(chain$samples)) pointcol = pointcol[chain$samples]
	          points(chain$data[chain$samples,c(col,row)],pch=pch,col=pointcol, cex=pcex)
	          points(chain$means[col], chain$means[row],col=chain$pointcol,pch=chain$pty,cex=2*ccex)
          }
          box()
        }
      }
      right = j == Npar
     	if(j == 1 || right) {
     		side= 2+2*right
     		magaxis(side,cex.axis=1.25,mgp=c(0,(i%%2)+0.5*right,0),cex.axis=1)
     		axtext = plotcols[i]
     		if((j == 1 && i == 1) || (right && i == Npar)) axtext = "PDF"
     		mtext(axtext,side=side, padj=(-2.5-0.5*right)*(-1)^(right))
    	}
			top = i == Npar
			if(i == 1 || top) {
     		side= 1+2*top
     		magaxis(side,cex.axis=1.25,mgp=c(0,0.9*(j%%2)+0.5*!top,0),cex.axis=1)
     		axtext = plotcols[j]
     		mtext(axtext,side=side, padj=(2.3+0.5*!top)*(-1)^(top))
    	}
    }
  }
  if(plotlegend)
  {
  	#basenames = c("Samples","Mean","Mean±1σ")
  	namespre = c("","<","<")
  	namespost = c("",">",">±1σ")
  	allnames = c()
  	ltys = c()
  	pchs = c()
  	ncols = 6
  	pcols = c()
  	for(chaini in 1:nchains)
  	{
  		chain = chains[[chaini]]
  		if(!chain$isvec)
  		{
  			allnames = c(allnames,
					paste0(namespre,chain$name,namespost),
					"50%","68%","95%")
  			ltys = c(ltys,NA,chain$lty,chain$ltys)
  			pchs = c(pchs,chain$pointtype,chain$pty,rep(NA,4))
  			pcols = c(pcols, chain$pointcol, rep(chain$linecol,5))
  		}
  	}
  	nlines = 0
  	for(chaini in 1:nchains)
  	{
  		chain = chains[[chaini]]
  		if(chain$isvec)
  		{
  			allnames = c(allnames,chain$name)
  			ltys = c(ltys,chain$lty)
  			pchs = c(pchs,chain$pty)
  			pcols = c(pcols, chain$linecol)
  			nlines = nlines + 1
  		}
  	}
  	if(nlines > 0)
  	{
  		minord = length(allnames)
  		while((minord %% ncols) != 0) minord = minord + 1
  		topad = minord - length(allnames)
  		if(topad > 0)
  		{
  			ltys = c(ltys,rep(NA,topad))
  			pchs = c(pchs,rep(NA,topad))
  			pcols = c(pcols,rep(NA,topad))
  			allnames = c(allnames,rep(NA,topad))
  		}
  	}
  	ord = t(matrix(1:length(allnames),ncols,length(allnames)/ncols))
  	nrows = dim(ord)[1]
  	# What a pain in the @#%$
  	# See http://stackoverflow.com/questions/21262472/adjust-spacing-between-text-in-horizontal-legend
  	namemat =  matrix(allnames[ord],nrow=nrows,ncol=ncols)
  	namemat = apply(X = namemat, MARGIN = c(1,2), nchar)
  	namemat[!is.finite(namemat)] = 0
  	textwidths = (apply(namemat,MARGIN=c(2),max)+2)
  	# Don't ask me why 1.5. The legend has some margins that I have had no luck removing, possibly set in layout.
  	textwidths = matrix(cumsum(textwidths)-textwidths,nrows,ncols,byrow=TRUE)/(1.5*sum(textwidths))
  	print(textwidths[1,])
  	plot(c(),c(),xlim=c(0,1),ylim=c(0,1),axes = F,xlab = '', ylab='',mar=c(0,0,0,0))
  	textwidths = textwidths/matrix((1:ncols)-1,nrows,ncols,byrow=TRUE)
		textwidths[,1] = 0
		textwidths = as.vector(textwidths)
		#  0.01465461 0.01465461 0.01465461 0.21412500 0.21412500 0.21412500 0.41359539 0.41359539 0.41359539 0.61306579 0.61306579 0.61306579 0.81253618 0.81253618 0.81253618
# [16] 1.01200658 1.01200658 1.01200658
		#legend2(x=0, y=0.5, legend=allnames[ord],col=pcols[ord],lty=ltys[ord],pch=pchs[ord],ncol=ncols, cex=1.5,bty="n")
		print(pchs)
		legend(x=0, y=(legheight-0.3)/legheight, legend=allnames[ord],
			col=pcols[ord],lty=ltys[ord],pch=pchs[ord],
			text.width=textwidths, ncol=ncols, cex=1.5,bty="n")
  }
  # No output for now
}
