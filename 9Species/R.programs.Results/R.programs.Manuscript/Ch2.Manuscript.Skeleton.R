##Analyzing the results from the MSP Model

rm(list=ls())
ls()

Run <- 'Run83_Mod_With'

my.programs <- c("C:\\Users\\kiersten.curti\\My_Programs\\")
run.dir <- c(paste(my.programs,"\\ADMB\\9Species\\",Run,sep=""))
setwd(run.dir)
load("OrganizedResults.RData")

fig.type <- c('wmf')

#fig.dir <- c("C:\\My_Programs\\Dissertation_Outputs\\Chapter2_Full_GB_Model\\Figures")
#table.dir <- c("C:\\My_Programs\\Dissertation_Outputs\\Chapter2_Full_GB_Model\\Tables")
fig.dir <- paste(run.dir,"\\Summary_Figures",sep="")
table.dir <- paste(run.dir,"\\Summary_Tables",sep="")



#Figure height and widths for figs of par(mfcol=c(3,3))
#fig.ht.3 <- 6.8 # Original - with middle x-axes
fig.ht.3 <- 6.2 # Middle x-axes removed
fig.wdth.3 <- 6.2
#Figure height and widths for figs of par(mfcol=c(2,3))
fig.ht.2 <- 4.35 # Middle x-axes removed
fig.wdth.2 <- 6.2
# Figure height and widths for figs of par(mfcol=c(4,3))
fig.ht.4 <- 8.1 # Middle x-axes removed
fig.wdth.4 <- 6.2


#line/point details
l.width <- 1.5
p.size <- 1.5
p.width <- 0.8
a.width <- 1
a.size <- 1.0

label.size <- 0.8
fig.name.cex <- 0.8

lty.key <- 1:10

# ID'g which species had length structured ssp runs
length.struct <- rep(0,length(sp.names))
  names(length.struct) <- sp.names
length.struct[c('SDog','WSk')] <- 1

# Species order
sp.order <- 1:9
  names(sp.order) <- sp.names

# Number years
nyr <- length(unique(yrs))

library(TeachingDemos)  #For figure legend
type.key <- c('a)' , 'b)' , 'c)')
key.x.coord <- -0.45
key.y.coord <- 1.11


###########################################################

prey.sp <- names(max.M2)[max.M2!=0]
pred.sp <- names(sp.order)[unique(pred)[order(unique(pred))]]
nint <- length(pred)

# load Brewer color and Gplots
library(RColorBrewer)  # Color schemes
library(gplots)  # Barplot2
library(Hmisc) # Errbar

# Consumption
source(paste(my.programs,"\\R\\Code_for_Chapters\\Chapter2_Results\\OrganizeConsumedB.R",sep="") )
# Organize consum 4darray;

# Fits to input data
source(paste(my.programs,"\\R\\Code_for_Chapters\\Chapter2_Results\\Ch2Figs.FitsToInputData.R",sep="") )
  # Plot observed inputs and predicted fits for select data

# Species preference coefficients
source(paste(my.programs,"\\R\\Code_for_Chapters\\Chapter2_Results\\Organize.SpeciesPrefCoeffs.R",sep="") )
# Create list of parameter estimates and corresponding standard deviations

# Predicted indices
source(paste(my.programs,"\\R\\Code_for_Chapters\\Chapter2_Results\\Ch2Figs.Predictions.R",sep="") )
# Plot predicted indices: M2, Consumption figures

# Additional figures: Food for thought
source(paste(my.programs,"\\R\\Code_for_Chapters\\Chapter2_Results\\Ch2Figs.AdditionalFigures.R",sep="") )
# Plot additional FH figures; Also includes additional data manipulations


###########################################################

# Size preference curves
  # This program stands alone
# source("C:\\My_Programs\\R\\Code_for_Chapters\\Chapter2_Results\\PlotSizePreference.R")





