###########################################################
# FitsToInputData.R
# Plot observed inputs and predicted fits for select data

###########################################################


##### TotC
setwd(fig.dir)
if(fig.type == 'gui')    {dev.new(height = fig.ht.3, width = fig.wdth.3)}
if(fig.type == 'wmf') {win.metafile("./TotC_AllSp.wmf" , height = fig.ht.3, width = fig.wdth.3)}
par(mfrow=c(3,3))
par(mar=c(2, 2.6, 0.3, 1) +0.1);  par(oma=c(2.5,2.5,2.3,0)) # Horizontal y-axis
for (i in 1:nsp)  {
  par(las = 1) # horizontal text
  ytmp<-c(obs.totc[[i]],pred.totc[[i]])
    plot(yrs[i,],obs.totc[[i]],type="p",xlab="",ylab="",xlim=yr.lim,ylim=c(0,max(ytmp)),axes=F,cex=p.size,lwd=p.width)
    lines(yrs[i,],pred.totc[[i]],col=msp.col,lwd=l.width, lty = lty.key[msp.lty])
    axis(side=2, at=axTicks(2), labels=TRUE, cex.axis=a.size,lwd=a.width)
    axis(side=1, at=yr.ticks, labels=FALSE, cex.axis=a.size,lwd=a.width)
    if(i >= 7)  {axis(side=1, at=yr.ticks, labels=TRUE, cex.axis=a.size,lwd=a.width)}
  mtext(fig.names[i], side=3, cex=fig.name.cex, line=0.3, font=2)
  box(lwd=a.width)
  } # end of species loop
  par(las = 0)
  mtext("Year",side=1,line=1,cex=label.size,outer=T)
  mtext("Thousands of metric tons", side=2, line=0.7, cex=label.size,outer=T) 
  if (fig.type == 'wmf') {dev.off()}



###########################################################
##### TotFIC

fig.cols <- 3

setwd(fig.dir)
if(fig.type == 'gui')    {dev.new(height = fig.ht.4, width = fig.wdth.4)}
if(fig.type == 'wmf') {win.metafile("./TotFIC.wmf" , height = fig.ht.4, width = fig.wdth.4)}
par(mfrow=c(4,fig.cols))
par(mar=c(3, 2.3, 1.0, 1) +0.1);  par(oma=c(1.5,3.3,1.5,0)) # Horizontal y-axis
for (j in 1:sp9.rep$nFIC)  {
  # Determine species to plot for specific season
    # Select age-structured species with TSwt > 0 in survey j
  TSseas <- TSwt[length.struct==0,j]
  TSpos <- TSseas[TSseas!=0]
  sp.seas <- sp.order[names(TSpos)]
  nplots <- ceiling(length(sp.seas)/fig.cols)*fig.cols

  for (spct in 1:length(sp.seas))
    {
    i <- as.integer(sp.seas[spct])
    sp <- names(sp.seas)[spct]
    par(las = 1)
    obs.tmp <- obs.totfic[[i]][j,]
    pred.tmp <- pred.totfic[[i]][j,]
    ytmp<-c(obs.tmp,pred.tmp)
    plot  (yrs[i,],obs.tmp,type="p",xlab="",ylab="",xlim=yr.lim,ylim=c(0,max(ytmp)),axes=F,cex=p.size,lwd=p.width)
    lines (yrs[i,],pred.tmp,col=msp.col,lwd=l.width, lty = lty.key[msp.lty])
    axis(side=2, at=axTicks(2), labels=TRUE, cex.axis=a.size,lwd=a.width)
    axis(side=1, at=yr.ticks, labels=TRUE, cex.axis=a.size,lwd=a.width)
    box(lwd=a.width)
    mtext(fig.names[i], side=3, cex=fig.name.cex, line=0.3, font=2)
    if (spct == 1)  {
      par(xpd=NA) #now can write text in outer margin areas
			par(las = 0)
			key.coords <-cnvrt.coords(key.x.coord,key.y.coord,'plt')$usr  #converts from plt to usr coordinate system
			text(key.coords,labels=type.key[j], font=2, cex=1.2)
			}
    } # end of sp.seas loop

  for (i in 1:(nplots-length(sp.seas)))  {
    plot(yrs[i,],rep(1,length(yrs[i,])),type="l",col=par("bg"),axes=F,xlab="",ylab="")  # leave empty slot
    }  # end of empty loop
  } # end of survey loop
par(las = 0)
mtext("Year",side=1,line=0,cex=label.size,outer=T)
mtext("Number/tow", side=2, line=0.5, cex=label.size,outer=T) 

if (fig.type == 'wmf') {dev.off()}
      
  

