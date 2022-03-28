

Greys <- colors()[round(seq(from=153,to=227,length=10),0)]

###########################################################
##### M2
M2.fig.ht <- fig.ht.2 + 0.3
setwd(fig.dir)
if(fig.type == 'gui')    {dev.new(height = M2.fig.ht, width = fig.wdth.2)}
if(fig.type == 'wmf') {win.metafile("./M2_AllPrey.wmf" , height = M2.fig.ht, width = fig.wdth.2)}
par(mfrow=c(2,3))
par(mar=c(2, 2.6, 0.3, 1) +0.1);  par(oma=c(5.0,2.5,2.3,0)) # Horizontal y-axis
for (sp in 1:length(prey.sp))  {
  py <- prey.sp[sp]
  y.range = c(0,round(max.M2[[py]],3))
  par(las = 1) # horizontal text
  plot( yrs[py,],M2[[py]][,1], type="l",col=Greys[1],lty=lty.key[1],  xlab="",ylab="",xlim=yr.lim,ylim=y.range,lwd=l.width,axes=F)
  for (a in 2: nage[py])  {
    lines(yrs[py,],M2[[py]][,a],col='blue',type="l",lwd=l.width, lty=lty.key[a])
    lines(yrs[py,],M2[[py]][,a],col=Greys[a],type="l",lwd=l.width, lty=lty.key[a])
    }  # end of age loop
  axis(side=1, at=yr.ticks, labels=FALSE, cex.axis=a.size,lwd=a.width)
  axis(side=2, at=axTicks(2), labels=TRUE, cex.axis=a.size,lwd=a.width)
  if(sp >= 4)  {axis(side=1, at=yr.ticks, labels=TRUE, cex.axis=a.size,lwd=a.width)}
  mtext(fig.names[py], side=3, cex=fig.name.cex, line=0.3, font=2)
  box(lwd=a.width)
  if(sp==5 ) {
    par(xpd=NA) #now can write text in outer margin areas
    legend(x =(yrs[1]-45), y=-0.73, 1:max(nage[prey.sp]), bty="n", col=Greys[1:max(nage[prey.sp])],lty=lty.key[1:max(nage[prey.sp])],cex=1.0,lwd=1,ncol=max(nage[prey.sp])) 
    }
  } # end of species loop
  par(las = 0)
  mtext("Year",side=1,line=1,cex=label.size,outer=T)
  mtext("Predation mortality", side=2, line=1, cex=label.size,outer=T) 
  if (fig.type == 'wmf') {dev.off()}



###########################################################
##### Compare Conusmed biomass and fisheries catch
# sumCon.sp = mil of kg = tmt = TotC
sumCon.sp <- lapply(sum.Con,rowSums)
setwd(fig.dir)
if(fig.type == 'gui')    {dev.new(height = fig.ht.2, width = fig.wdth.2)}
if(fig.type == 'wmf') {win.metafile("./Consum.Vs.Catch.wmf" , height = fig.ht.2, width = fig.wdth.2)}
par(mfrow=c(2,3))
par(mar=c(2, 2.6, 0.3, 1) +0.1);  par(oma=c(2.3,2.5,2.3,0)) # Horizontal y-axis
for (sp in 1:length(prey.sp))  {
  py <- prey.sp[sp]
  ytmp <- c(obs.totc[[py]] , sumCon.sp[[py]])
  par(las = 1) # horizontal text
  plot(names(obs.totc[[py]]), obs.totc[[py]], type="l", lty=1, xlab="", ylab="", xlim=yr.lim, ylim=c(0,max(ytmp)), lwd=l.width, axes=F)
  lines(names( sumCon.sp[[py]]),  sumCon.sp[[py]], type="l", lty=2, lwd=l.width)
  axis(side=1, at=yr.ticks, labels=FALSE, cex.axis=a.size,lwd=a.width)
  axis(side=2, at=axTicks(2), labels=TRUE, cex.axis=a.size,lwd=a.width)
  if(sp >= 4)  {axis(side=1, at=yr.ticks, labels=TRUE, cex.axis=a.size,lwd=a.width)}
  mtext(fig.names[py], side=3, cex=fig.name.cex, line=0.3, font=2)
  box(lwd=a.width)
  } # end of species loop
  par(las = 0)
  mtext("Year",side=1,line=1,cex=label.size,outer=T)
  mtext("Thousands of metric tons", side=2, line=1, cex=label.size,outer=T) 
  if (fig.type == 'wmf') {dev.off()}



###################################################################
## Plot biomass of modeled prey species consumed by predator species in each year
# Plot based on consum.pysp

#display.brewer.pal(8,"RdYlBu")
bar.col <- brewer.pal(8,"RdYlBu")[-3]
setwd(fig.dir)
if(fig.type == 'gui')    {dev.new(height = fig.ht.4, width = fig.wdth.4)}
if(fig.type == 'wmf') {win.metafile("./ConsumedB_ByPred.wmf" , height = fig.ht.4, width = fig.wdth.4)}
par(mfrow=c(3,3))
par(mar=c(4.5, 2.0, 0.0, 1) +0.1);  par(oma=c(0.5,2.5,2.3,0)) # Horizontal y-axis
for (i in 1:length(pred.sp))  {
  par(las = 1) # horizontal text
  pd <- pred.sp[i]
  tmp <- consum.pysp[[pd]]
    rownames(tmp) <- as.character(fyr[pd]:lyr[pd])
    colnames(tmp) <- c(fig.names,"Other")
  barplot2(t(tmp)[fig.names[prey.sp],], legend.text=F, col=bar.col, axis.lty=1) #, density=seq(from=20,by=15,length=length(prey.sp)))
  box(lwd=a.width)
  mtext(fig.names[pd], side=3, cex=fig.name.cex, line=0.3, font=2)
  } # end of species loop

  # Add legend:  Make fake barplot and then add legend so that it can be positioned correctly
  tmp[]<-NA
  barplot2(t(tmp)[fig.names[prey.sp],], axes=F, axisnames=F, col=bar.col, ylim=c(5,15))#, density=seq(from=20,by=15,length=length(prey.sp)))
  par(xpd=NA) #now can write text in outer margin areas  
  legend("center",fig.names[prey.sp], col=bar.col, pch=15, bty="n", border="black") # density=seq(from=20,by=15,length=length(prey.sp)) , bty="n")

  par(las = 0)
  mtext("Year",side=1,line=-1.5,cex=label.size,outer=T)
  mtext("Thousands of metric tons", side=2, line=0.9, cex=label.size,outer=T) 
  if (fig.type == 'wmf') {dev.off()}



###########################################################
# 2-15-12: Modified code to plot the biomass of each species and consumption of just predator sp
# Total consumption by predator species (including other food) and predator B
# sum.consum.yr[pdsp,yr]
# TotB[[sp]]

setwd(fig.dir)
if(fig.type == 'gui')    {dev.new(height = fig.ht.3, width = fig.wdth.3)}
if(fig.type == 'wmf') {win.metafile("./PdConsum.And.TotB.wmf" , height = fig.ht.3, width = fig.wdth.3)}
par(mfrow=c(3,3))
par(mar=c(2, 2.6, 0.3, 1) +0.1);  par(oma=c(2.5,2.5,2.3,0)) # Horizontal y-axis

for (i in 1:nsp)  {
  pd <- sp.names[i]
  ytmp <- c( sum.consum.yr[pd,] , TotB[[pd]] )
  par(las = 1) # horizontal text

  plot(names(TotB[[pd]]), TotB[[pd]], type="l", lty=1, xlab="", ylab="", xlim=yr.lim, ylim=c(0,max(ytmp)), lwd=l.width, axes=F)
  if(length(pred.sp[pred.sp==pd])>0)  {  lines( names(TotB[[pd]]),  sum.consum.yr[pd,], type="l", lty=2, lwd=l.width) }
  axis(side=1, at=yr.ticks, labels=FALSE, cex.axis=a.size,lwd=a.width)
  axis(side=2, at=axTicks(2), labels=TRUE, cex.axis=a.size,lwd=a.width)
  if(i > 6)  {axis(side=1, at=yr.ticks, labels=TRUE, cex.axis=a.size,lwd=a.width)}
  mtext(fig.names[pd], side=3, cex=fig.name.cex, line=0.3, font=2)
  box(lwd=a.width)
  } # end of species loop
par(las = 0)
mtext("Year",side=1,line=1,cex=label.size,outer=T)
mtext("Thousands of metric tons", side=2, line=1, cex=label.size,outer=T) 
if (fig.type == 'wmf') {dev.off()}


# **Check: read cb (CB Ratios) for a particular predator, pd, from clipboard
# pd <- 'Cod'
# max(abs(rowSums(cb*B[[pd]]) - sum.consum.yr[pd,]))



###########################################################
# Species preference coefficients for each predator species
# 02-15-12: Per Chris's request, 1) Modified species order within the legend box
#                                       and 2) made all y-axes have the same range

setwd(fig.dir)
if(fig.type == 'gui')    {dev.new(height = (fig.ht.4-0.7), width = fig.wdth.4)}
if(fig.type == 'wmf') {win.metafile("./Species.Pref.Coeffs.wmf" , height = (fig.ht.4-0.7), width = fig.wdth.4)}
par(mfrow=c(3,3))
par(mar=c(4.0, 2.6, 0.3, 0) +0.1);  par(oma=c(1.0,1.7,2.3,0.3)) # Horizontal y-axis
par(las=1)
yaxis.max <- do.call(rbind,Rho.pdpy)$value+do.call(rbind,Rho.pdpy)$std
for (i in 1:length(Rho.pdpy))  {
  # i <- 1
  pd <- names(Rho.pdpy)[i]
  rho.tmp <- rbind(Rho.pdpy[[i]], t(as.matrix(rho.other)))
  npoints <- nrow(rho.tmp)
  x.axis  <- seq(from = 1, to = npoints, by = 1)
  x.labels <- as.vector(rho.tmp$prey.names)
  xmin <- min(x.axis) - 0.25
  xmax <- max(x.axis) + 0.25
  y.axis <- as.numeric(rho.tmp$value)
  y.max <- y.axis + as.numeric(rho.tmp$std)
  y.min <- y.axis - as.numeric(rho.tmp$std)

  plot(x.axis, y.axis, pch = "-", xlim = c(xmin,xmax), ylim=c(0,max(yaxis.max,na.rm=T)), lwd=2, cex=2, axes=F, xlab="", ylab="")
  errbar(x.axis, y.axis, y.max, y.min, add=TRUE, pch="", cap=0.02)
  axis(side=2, at=axTicks(2), labels=TRUE, cex.axis=a.size,lwd=a.width)
  axis(side=1, at=x.axis, labels = x.labels, cex.axis=(a.size),lwd=a.width, padj = -0.2)
  mtext(fig.names[pd], side=3, cex=fig.name.cex, line=0.3, font=2)
  box()
  } # end of predator species loop
par(las=0)
# Add legend:  Make fake barplot and then add legend so that it can be positioned correctly
plot(1:10,1:10, type="n", axes=F,xlab="",ylab="")
par(xpd=NA) #now can write text in outer margin areas  
legend("center", c(
                              "G = Goosefish",
                              "C = Cod",
                              "M = Mackerel",
                              "W = White hake",
                              "S = Silver hake",
                              "H = Herring",
                              "O = Other food"
                              ),
    col=bar.col, bty="n", border="black", title="Prey species") # density=seq(from=20,by=15,length=length(prey.sp)) , bty="n")
mtext("Prey species", side=1, line= -0.7, cex=label.size,outer=T)
mtext("Log (species-preference coefficient)", side=2, line=0.1, cex=label.size,outer=T)
if(fig.type == 'wmf' | fig.type == 'pdf') {dev.off()}


  
  
  
  
   