###################################################################
### Extra data manipulation ###

## Calculate average diet
    #Average over years to obtain an average obs diet (by age)
    # avg.pred.fh.age was previously calculated within Organizing.MSP.Outputs
  
colMeans.na <- function(x){colMeans(x, na.rm=T)}
avg.obs.fh.age <- list.bysp
 for (i in 1:nsp)  {
  avg.obs.fh.age[[i]] <- t(apply(obs.fh.age[[i]],3,colMeans.na))
  rownames(avg.obs.fh.age[[i]]) <- as.character(1:nage[i])
  colnames(avg.obs.fh.age[[i]]) <- FH.sp.names
  }

# All prey names, including other food
FH.fig.names <- c(fig.names,"Other")
  names(FH.fig.names) <- FH.sp.names
prey.sp.other <- c(prey.sp,"Other")



###################################################################
## Plot %Wt of modeled prey species consumed by predator species in each year.....
  # i.e. as proportion of total annual consumption
# Plot based on sc.consum.pysp

#display.brewer.pal(8,"RdYlBu")
bar.col <- brewer.pal(8,"RdYlBu")[-3]
setwd(fig.dir)

if(fig.type == 'gui')    {dev.new(height = fig.ht.4, width = fig.wdth.4)}
if(fig.type == 'wmf') {win.metafile("./Scaled_ConsumedB_ByPred.wmf" , height = fig.ht.4, width = fig.wdth.4)}
par(mfrow=c(3,3))
par(mar=c(4.5, 2.0, 0.0, 1) +0.1);  par(oma=c(0.5,2.5,2.3,0)) # Horizontal y-axis
for (i in 1:length(pred.sp))  {
  par(las = 1) # horizontal text
  pd <- pred.sp[i]
  tmp <- sc.consum.pysp[[pd]]
    rownames(tmp) <- as.character(fyr[pd]:lyr[pd])
    colnames(tmp) <- fig.names
  barplot2(t(tmp)[fig.names[prey.sp],], legend.text=F, col=bar.col, border=colors()[176], axis.lty=1) #, density=seq(from=20,by=15,length=length(prey.sp)))
  box(lwd=a.width)
  mtext(fig.names[pd], side=3, cex=fig.name.cex, line=0.3, font=2)
  } # end of species loop

# Add legend:  Make fake barplot and then add legend so that it can be positioned correctly
tmp[]<-NA
barplot2(t(tmp)[fig.names[prey.sp],], axes=F, axisnames=F, col=bar.col, ylim=c(5,15))#, density=seq(from=20,by=15,length=length(prey.sp)))
par(xpd=NA) #now can write text in outer margin areas  
legend("center",fig.names[prey.sp], col=bar.col, pch=15, bty="n", border="black", title="Prey species") # density=seq(from=20,by=15,length=length(prey.sp)) , bty="n")

par(las = 0)
mtext("Year",side=1,line=-1.5,cex=label.size,outer=T)
mtext("Proportion by weight", side=2, line=1.3, cex=label.size,outer=T) 
if (fig.type == 'wmf') {dev.off()}



###################################################################
### From the perspective of the prey species
# Proportion of total biomass consumed by each predator

#display.brewer.pal(8,"RdYlBu")
bar.col <- brewer.pal(8,"RdYlBu")[-3]
setwd(fig.dir)

if(fig.type == 'wmf') {win.metafile("./BiomConsumed_ByPrey.wmf" , height = fig.ht.4, width = fig.wdth.4)}
if(fig.type == 'gui')    {dev.new(height = fig.ht.4, width = fig.wdth.4)}
par(mfrow=c(3,3))
par(mar=c(4.5, 2.0, 0.0, 1) +0.1);  par(oma=c(0.5,2.5,2.3,0)) # Horizontal y-axis
for (i in 1:length(prey.sp))  {
  par(las = 1) # horizontal text
  py <- prey.sp[i]
  tmp <- consum.pdsp[,,py]
    rownames(tmp) <- as.character(fyr[py]:lyr[py])
    colnames(tmp) <- fig.names
  barplot2(t(tmp)[fig.names[pred.sp],], legend.text=F, col=bar.col, axis.lty=1) #, density=seq(from=20,by=15,length=length(prey.sp)))
  box(lwd=a.width)
  mtext(fig.names[py], side=3, cex=fig.name.cex, line=0.3, font=2)
  } # end of species loop

# Add legend:  Make fake barplot and then add legend so that it can be positioned correctly
tmp[]<-NA
barplot2(t(tmp)[fig.names[pred.sp],], axes=F, axisnames=F, col=bar.col, ylim=c(5,15))#, density=seq(from=20,by=15,length=length(prey.sp)))
par(xpd=NA) #now can write text in outer margin areas  
legend("center",fig.names[pred.sp], col=bar.col, pch=15, bty="n", border="black",title="Predator species") # density=seq(from=20,by=15,length=length(prey.sp)) , bty="n")

par(las = 0)
mtext("Year",side=1,line=-1.5,cex=label.size,outer=T)
mtext("Thousands of metric tons", side=2, line=1.3, cex=label.size,outer=T) 
if (fig.type == 'wmf') {dev.off()}

  
  
###################################################################
## Plot avg observed and predicted diet
# Plots based on avg.obs.fh.age and avg.pred.fh.age

# Arrange observed and predicted datasets for loop
fh.fig.names <- c("./obs.AvgDiet.at.age.wmf","./pred.AvgDiet.at.age.wmf")
  names(fh.fig.names) <- c("obs","pred")
fh.avg.datasets <- c('avg.obs.fh.age','avg.pred.fh.age')
  names(fh.avg.datasets) <- names(fh.fig.names)

# Loop over obs/pred
for (ds in 1:length(fh.fig.names)) {
  diet.age <-  get(fh.avg.datasets[ds])  # avg.obs.fh.age or avg.pred.fh.age
  
  #display.brewer.pal(8,"RdYlBu")
  bar.col <- brewer.pal(8,"RdYlBu")[-3]
  setwd(fig.dir)
  if(fig.type == 'gui')    {dev.new(height = fig.ht.4, width = fig.wdth.4)}
  if(fig.type == 'wmf') {win.metafile(fh.fig.names[ds] , height = fig.ht.4, width = fig.wdth.4)}
  par(mfrow=c(3,3))
  par(mar=c(4.5, 2.0, 0.0, 1) +0.1);  par(oma=c(0.5,2.5,2.3,0)) # Horizontal y-axis

  # Loop over predator sp
  for (i in 1:length(pred.sp))  {
    par(las = 1) # horizontal text
    pd <- pred.sp[i]

    tmp <- diet.age[[pd]]
      colnames(tmp) <- c(fig.names,"Other")
    barplot2(t(tmp)[FH.fig.names[prey.sp.other],], legend.text=F, col=bar.col, axis.lty=1)  #, density=seq(from=20,by=15,length=length(prey.sp.other)))
    box(lwd=a.width)
    mtext(fig.names[pd], side=3, cex=fig.name.cex, line=0.3, font=2)
    } # end of species loop

  # Add legend:  Make fake barplot and then add legend so that it can be positioned correctly
  tmp[]<-NA
  barplot2(t(tmp)[FH.fig.names[prey.sp.other],], axes=F, axisnames=F, col=bar.col, ylim=c(5,15))#, density=seq(from=20,by=15,length=length(prey.sp)))
  par(xpd=NA) #now can write text in outer margin areas  
  legend("center",FH.fig.names[prey.sp.other], col=bar.col, pch=15, bty="n", border="black", title="Prey species") # density=seq(from=20,by=15,length=length(prey.sp)) , bty="n")

  par(las = 0)
  mtext("Age",side=1,line=-1.5,cex=label.size,outer=T)
  mtext("Proportion by weight", side=2, line=0.9, cex=label.size,outer=T) 
  if (fig.type == 'wmf') {dev.off()}

  } # End of datset (ds) loop
  




 