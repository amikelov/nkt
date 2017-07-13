library(tcR)
library(dplyr)
library(openxlsx)
source("nkt_functions.R")
load("all_clonesets.RData")

## parsing
YB.CMV <- parse.folder("/home/artem/Desktop/IBH/NKT/new_nkt/fractions/TRB/YB_CMV/",.format = "mixcr")
YB.CMV <- set.rank(YB_CMV)

samples <- c('DL','IrZv', 'TA','YB','TV','Kar','Kov', 'LY')

for(samp in samples) { 
  assign(paste0(samp,".nkt"),
         parse.folder(paste0('/home/artem/Desktop/IBH/NKT/new_nkt/fractions/TRB/',samp), 
                      .format = 'mixcr'))
  assign(paste0(samp,".nkt"),set.rank(get(paste0(samp,".nkt"))))
}

for(samp in samples) { 
  assign(samp,
         parse.folder(paste0('/home/artem/Desktop/IBH/NKT/new_nkt/fulls/export/',
                            samp), 
                      .format = 'mixcr'))
  assign(samp,set.rank(get(samp)))
}




## YB - processing
YB.intersect <- nkt.intersect(YB.nkt,YB,"YB")
YB.intersect[,1:6] <- NULL
YB.intersect[,2] <- NULL
YB.intersect <- YB.intersect[,c(2,1,3:ncol(YB.intersect))]

## IFN - crossess
IFNhi2neg.NKT <- merge(YB.CMV[[1]],YB.intersect, by = "CDR3.nucleotide.sequence",all.x = T)
IFNmed2neg.NKT <- merge(YB.CMV[[2]],YB.intersect, by = "CDR3.nucleotide.sequence",all.x =T)
IFNneg2neg.NKT <- merge(YB.CMV[[3]],YB.intersect, by = "CDR3.nucleotide.sequence",all.x =T)
IFNneg2pos.NKT <- merge(YB.CMV[[4]],YB.intersect, by = "CDR3.nucleotide.sequence",all.x =T)

hs1 <- createStyle(fgFill = "#DCE6F1", halign = "CENTER", textDecoration = "bold", border = "TopBottom")
write.xlsx(IFNhi2neg.NKT,"/home/artem/Desktop/IBH/NKT/new_nkt/analysis/IFN/IFNhi2neg.NKT.xlsx",headerStyle=hs1)
write.xlsx(IFNmed2neg.NKT,"/home/artem/Desktop/IBH/NKT/new_nkt/analysis/IFN/IFNmed2neg.NKT.xlsx",headerStyle=hs1)
write.xlsx(IFNneg2neg.NKT,"/home/artem/Desktop/IBH/NKT/new_nkt/analysis/IFN/IFNneg2neg.NKT.xlsx",headerStyle=hs1)
write.xlsx(IFNneg2pos.NKT,"/home/artem/Desktop/IBH/NKT/new_nkt/analysis/IFN/IFNneg2pos.NKT.xlsx",headerStyle=hs1)

IFNhi2neg.NKT<-IFNhi2neg.NKT[,c(1,6:7,4:5,17,25,23,24,22,26,21,20,19,33,34,27)]
IFNmed2neg.NKT <- IFNmed2neg.NKT[,c(1,4:7,17,25,23,24,22,26,21,20,19)]
IFNneg2neg.NKT <- IFNneg2neg.NKT[,c(1,4:7,17,25,23,24,22,26,21,20,19)]
IFNneg2pos.NKT <- IFNneg2pos.NKT[,c(1,4:7,17,25,23,24,22,26,21,20,19)]

## all processing
for (samp in samples){
  assign(paste0(samp,".intersect"),
         nkt.intersect(get(paste0(samp,".nkt")),
                       get(samp),
                       samp))
}

## add p.values
  
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
clusterExport(cl,"Kar")
clusterExport(cl,"Kar.nkt")
Kar.nkt.p <- lapply(Kar.nkt,function(x){ clone.count.prob(full = Kar$Kar_F, fraction = x, adj.method = "fdr",cl)} )
stopCluster(cl)


for (samp in samples[samples!="LY" & samples!="Kov"]){
  
  assign(paste0(samp,".nkt"),
         lapply(get(paste0(samp,".nkt")),
                function(x){
                  print()
                  return(clone.count.prob(get(samp)[3],x, "fdr"))
                  }
         ))
}