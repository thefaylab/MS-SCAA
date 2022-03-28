# PlotSizePreference.R

# Based off of:  C:\\Users\\Kiersten L Curti\\Documents\\Dissertation\\Data\\NEFSC Food Habits Data\\Indep_of_TimeStep\\PreyLengths-9Species\\SizePreferenceResults.R

rm(list=ls())
ls()

fig.type <- c('wmf')
fig.dir <- c("C:\\My_Programs\\Dissertation_Outputs\\Chapter2_Full_GB_Model\\Figures")
table.dir <- c("C:\\My_Programs\\Dissertation_Outputs\\Chapter2_Full_GB_Model\\Tables")


my.programs <- c("C:/My_Programs/")
path=paste(my.programs,"R/Functions/reptoRlist.R",sep="")
reptoRlist<-dget(path);

spref.dir <- c("C:/Users/Kiersten L Curti/Documents/Dissertation/Data/NEFSC Food Habits Data/Indep_of_TimeStep/PreyLengths-9Species/")
spref.names <- c("SDog" , "WSk" , "Goose" , "Cod" , "Pol" , "WHake" , "SHake")

fig.names <- c("Spiny dogfish", "Winter skate", "Goosefish", "Cod", "Pollock", "White hake", "Silver hake")
  names(fig.names) <- spref.names

#Figure height and widths for figs of par(mfcol=c(3,3))
spref.fig.ht.3 <- 6.8 # Middle x-axes removed
fig.wdth.3 <- 6.2

#line/point details
l.width <- 1.5
p.size <- 1.5
p.width <- 0.8
a.width <- 1
a.size <- 1.0

label.size <- 0.8
fig.name.cex <- 0.8


####### ------------------------------------------------------------------------------------------------------------- #######

bins <- vector('list',length(spref.names))
  names(bins) <- spref.names
prop <- bins
prophat <- bins
eta <- bins
sig1 <- bins
sig2 <- bins

for (i in 1:length(spref.names))  {
  sp <- spref.names[i]
  dataset <- paste(sp,"AllPy",sep="_") 
  sp.dir <- paste(spref.dir,dataset,sep="")
  setwd(sp.dir)
  pd.rep<-reptoRlist(paste(sp,"SizePref.rep",sep=""))
  
  bins[[sp]] <- pd.rep$bin
  prop[[sp]] <- pd.rep$Prop  
  prophat[[sp]] <- pd.rep$scProphat

  eta[[sp]] <- pd.rep$eta
  sig1[[sp]] <- pd.rep$sig1
  sig2[[sp]] <- pd.rep$sig2
  } # end of species loop


# Plot observed and predicted proportions, scaled
# 02-15-12: Modified figure so that each subplot has the same axis scale

setwd(fig.dir)

if(fig.type=='gui') {dev.new(height=spref.fig.ht.3,width=fig.wdth.3)}
if(fig.type=='wmf') {win.metafile("./SizePref_AllPredSp.wmf", height=spref.fig.ht.3, width=fig.wdth.3) }
par(mfrow=c(3,3))
par(mar=c(4, 2.6, 0.3, 1) +0.1);  par(oma=c(0.6,2.5,2.2,0)) # Horizontal y-axis
# par(mar=c(3, 2.6, 0.3, 1) +0.1);  par(oma=c(1.3,2.5,2.0,0)) # Horizontal y-axis

xmax <- max(unlist(bins))
ymax <- max(unlist(c(prop,prophat)))
for (i in 1:length(spref.names))  {
  sp <- spref.names[i]
  par(las = 1) # horizontal text
  # plot(bins[[sp]],prop[[sp]], type="p",xlab="",ylab="",axes=F, xlim=c(0,max(bins[[sp]])), ylim=c(0,max(c(prop[[sp]],prophat[[sp]]))), cex=p.size, lwd=p.width)
  plot(bins[[sp]],prop[[sp]], type="p",xlab="",ylab="",axes=F, xlim=c(0,xmax), ylim=c(0,ymax), cex=(p.size-0.5), lwd=p.width)
  lines(bins[[sp]],prophat[[sp]],lwd=l.width)
  axis(side=2, at=axTicks(2), labels=TRUE, cex.axis=a.size,lwd=a.width, hadj = 0.87)
  axis(side=1, at=axTicks(1), labels=TRUE, cex.axis=a.size,lwd=a.width, padj = -0.5)
  mtext(fig.names[i], side=3, line=0.3, cex=fig.name.cex,outer=F,font=2)	
  box(lwd=a.width)
  } # end of species loop

par(las = 0)
mtext("Log (Predator Wt / PreyWt)",side=1,line=-1.4,cex=label.size,outer=T)
mtext("Probability", side=2, line=1.1, cex=label.size,outer=T) 
if (fig.type == 'wmf') {dev.off()}



# Create table of size preference values

setwd(table.dir)
spref <- cbind(unlist(eta), unlist(sig1), unlist(sig2))
  colnames(spref) <- c("Eta","Sigma1","Sigma2")
  rownames(spref) <- fig.names[rownames(spref)]



write.csv(spref,"Size_Preference_Parameters.csv", row.names=T)


