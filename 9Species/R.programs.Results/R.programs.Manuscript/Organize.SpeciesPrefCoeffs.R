###########################################################
# Organize.rho.R
# Species preference coefficients for each predator species


setwd(run.dir)

#Reading in the standard deviation file
std.file <- read.table("9species.std", header=TRUE, fill=TRUE)
std.file$"dev" <- NULL

#Extract variable names
param.names <- as.vector(std.file$name)
param.names <- do.call(rbind, strsplit(param.names,"\\["))
param.names[,2] <- do.call(rbind, strsplit(param.names[,2],"\\]"))
  colnames(param.names) <- c("parameter","param.ct")
# Combine variable names with stdevs
param.stdevs <- cbind(param.names, std.file)

prey.names <- sp.names
  names(prey.names) <- sp.names
prey.names["Goose"] <- "G"
prey.names["Cod"] <- "C"
prey.names["Mack"] <- "M"
prey.names["Pol"] <- "P"
prey.names["WH"] <- "W"
prey.names["SH"] <- "S"
prey.names["Her"] <- "H"

# Extract rho estimates
Rho.stdevs <- param.stdevs[param.stdevs$"parameter"=="iRho",]
  col.tmp <- colnames(Rho.stdevs)
Rho.stdevs <- cbind(Rho.stdevs, pred,sp.names[pred], prey,prey.names[prey])
  colnames(Rho.stdevs) <- c(col.tmp, "pred","pred.names","prey","prey.names")
  

# Order by prey; split into lists by predator for plotting
ordered.tmp <- Rho.stdevs[order(Rho.stdevs$prey),]
Rho.pdpy <- split(ordered.tmp,ordered.tmp$pred)
  pd.tmp <- as.numeric(names(Rho.pdpy))
  names(Rho.pdpy)<- sp.names[pd.tmp]


# Create line for other food reference prey species coefficient
rho.other <- rep(NA,ncol(Rho.pdpy[[1]]))
  names(rho.other) <- colnames(Rho.pdpy[[1]])
  rho.other["value"] <- 0
  rho.other["std"] <- NA
  rho.other["prey.names"] <- "O"
  rho.other["prey"] <- "Other"


