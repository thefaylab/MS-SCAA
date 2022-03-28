consum.pysp

# load Brewer color and Gplots
library(RColorBrewer)  # Color schemes
library(gplots)  # Barplot2





display.brewer.pal(8,"RdYlBu")
#bar.col <- brewer.pal(7,"RdYlBu")
bar.col <- brewer.pal(8,"RdYlBu")[-3]
setwd(fig.dir)
if(fig.type == 'gui')    {dev.new(height = fig.ht.4, width = fig.wdth.4)}
if(fig.type == 'wmf') {win.metafile("./ConsumedB_ByPred.wmf" , height = fig.ht.4, width = fig.wdth.4)}
par(mfrow=c(3,3))
par(mar=c(4.5, 2.0, 0.3, 1) +0.1);  par(oma=c(0.5,2.5,2.3,0)) # Horizontal y-axis
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


# ----------------------------------- #

bar.col <- 'black'
bar.col <- brewer.pal(7,"RdYlBu")
setwd(fig.dir)
if(fig.type == 'gui')    {dev.new(height = fig.ht.4, width = fig.wdth.4)}
if(fig.type == 'wmf') {win.metafile("./ConsumedB_ByPred.wmf" , height = fig.ht.4, width = fig.wdth.4)}
par(mfrow=c(3,3))
par(mar=c(4.5, 2.0, 0.3, 1) +0.1);  par(oma=c(0.5,2.5,2.3,0)) # Horizontal y-axis
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

  # Add legend
  tmp[]<-NA
  barplot2(t(tmp)[fig.names[prey.sp],], legend=fig.names[prey.sp], axes=F, axisnames=F, col=bar.col, ylim=c(5,15),) # density=seq(from=20,by=15,length=length(prey.sp)))

  par(las = 0)
  mtext("Year",side=1,line=-1.5,cex=label.size,outer=T)
  mtext("Thousands of metric tons", side=2, line=0.9, cex=label.size,outer=T) 
  if (fig.type == 'wmf') {dev.off()}


