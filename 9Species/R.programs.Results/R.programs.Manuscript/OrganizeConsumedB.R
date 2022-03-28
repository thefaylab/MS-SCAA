#### Organize Consumed Biomass (consum 4darray in ADMB)


# Concatenate species and ages for rownames
sp.ages <- list.bysp
for (i in 1:nsp)  {  sp.ages[[i]] <- paste(sp.names[i],rep(1:nage[i]),sep=".")  }
sp.ages <- as.vector(c(unlist(sp.ages),"Other"))

consum <- list.bysp
for (i in 1:nsp)
  {
  consum[[i]] <- sp9.rep[[paste('consum.',i,sep='')]]
    yr.tmp <- rep(1:nyr,each=nage[i])
    age.tmp <- rep(1:nage[i],nyr)
    row.tmp <- paste("yr",yr.tmp,".",age.tmp,sep="")
    rownames(consum[[i]]) <- row.tmp
    colnames(consum[[i]]) <- sp.ages
  }
  
  
# dim(consum)
  # rows = nage[pd]*nyear
  # cols =  sum(nage)+1
    # SDog:  930 111
    # WSk:   630 111
    # Goose:  300 111
    # Cod:  300 111
    # Mack:  300 111
    # Pol:  270 111
    # WH:  210 111
    # SH:  180 111
    # Her:  180 111

# rem: 4darray Consum(1,nsp,FHfyr2,FHlyr2,1,nage,1,sage+1);



#########################################################
####### Consumption by predator species, summed over all prey  ########

# Sum over prey species/age to obtain annual consumption by each predator species/age
sum.consum.pdage <- lapply(consum,rowSums)  

# Sum over predator age to obtain annual consumption by each predator species
consumpt.sum.over.age <- function(con.tmp) {
  # con.tmp <- sum.consum.pdage[['WH']]
  yr.key <- do.call(rbind, strsplit(names(con.tmp), "\\."))[,1]
  con.tmp.yr <- lapply(split(con.tmp,yr.key),sum)
  con.tmp.yr <- unlist(con.tmp.yr)
  con.tmp.yr[paste("yr",1:30,sep="")] # Order rows
  }
sum.consum.yr <- lapply(sum.consum.pdage, consumpt.sum.over.age)
sum.consum.yr <- do.call(rbind, sum.consum.yr)



#########################################################
####### Consumption by predator species of each prey species ########

# Sum over predator age
consum.sum.pdage <- list.bysp
for (i in 1:nsp)  {
  pd <- sp.names[i]  #pd <- 'WH'
  con.tmp <- consum[[pd]]
  yr.key <- do.call(rbind, strsplit(rownames(con.tmp), "\\."))[,1]
  con.tmp.pd <- split(con.tmp, yr.key)
  con.tmp.pd.yr <- lapply(con.tmp.pd, function(tmp) {
                                                        tmp1 <- matrix(tmp,nrow=nage[pd],ncol=(sum(nage)+1))
                                                        colnames(tmp1) <- sp.ages
                                                        tmp2 <- colSums(tmp1) # Sum over pd age
                                                        tmp2
                                                        } )
  con.tmp.pd.yr<- do.call(rbind, con.tmp.pd.yr) # unlist
  con.tmp.pd.yr
  consum.sum.pdage[[i]] <- con.tmp.pd.yr                                                      
  } # end of species loop

# Split and sum by prey species
py.tmp <- do.call(rbind, strsplit(sp.ages,"\\."))[,1]
consum.pysp <- list.bysp
for (i in 1:nsp)  {
   pd <- sp.names[i]  
    con.tmp <- consum.sum.pdage[[pd]]
    yr.order <- rownames(con.tmp)
    con.tmp.pyage <- split(t(con.tmp),py.tmp)  # Split by prey sp
    con.tmp.py <- con.tmp.pyage  # object for sum over prey age    
    for (j in 1:nsp)  { # loop over prey
      py <- sp.names[j]
      tmp <- con.tmp.pyage[[py]]
      tmp1 <- t(matrix(tmp, nrow=nage[py], ncol=nyr))
      rownames(tmp1) <- yr.order
      tmp2 <- rowSums(tmp1)
      con.tmp.py[[py]] <- tmp2
      }  # end of prey loop
  # Annual consumption of each prey species by each predator species
  consum.pysp[[pd]] <- t(do.call(rbind, con.tmp.py))
  consum.pysp[[pd]] <- consum.pysp[[pd]][paste("yr",1:30,sep=""), FH.sp.names] # Order rows and columns
  } # end of species loop



#########################################################
####### Calculation check ########

# Calculate total by predator species for comparison to sum.consum.yr
  # they best be the same
con.bypd.tmp <- do.call(rbind, lapply(consum.pysp, rowSums))
sum.consum.yr - con.bypd.tmp
max(abs(sum.consum.yr - con.bypd.tmp))



#########################################################
####### Scaled consumption by predator species of each prey species ########
# Scale consum.pysp by year so that 
   # predator consumption **of modeled fish speices** within each year sums to 1
scale.byrow <- function(y){
  result <- t(apply(y,1,function(x){x[sp.names]/sum(x[sp.names])}))
  result}
sc.consum.pysp <- lapply(consum.pysp[pred.sp],scale.byrow)


#########################################################
####### From the prey's perspective.... 
  # Reshape consum.pysp: Aggregate predators for each prey species
  # Who are the primary predators of each prey species?

# Convert list to 3darray
consum.pysp.array <- array(unlist(consum.pysp),c(nrow(consum.pysp[[1]]),ncol(consum.pysp[[1]]),length(consum.pysp)))
  dimnames(consum.pysp.array)[[1]] <- rownames(consum.pysp[[1]])
  dimnames(consum.pysp.array)[[2]] <- colnames(consum.pysp[[1]])
  dimnames(consum.pysp.array)[[3]] <- names(consum.pysp)
# Reshape array so by prey species, not predator
consum.pdsp <- aperm(consum.pysp.array,c(1,3,2))
# Keep only prey species
consum.pdsp <- consum.pdsp[,,prey.sp]

