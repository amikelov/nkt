library(tcR)
library(dplyr)
library(openxlsx)
library(parallel)
options("openxlsx.numFmt" = NULL)

# without V.gene 
nkt.intersect <- function(nkt, fulls,sample_name) {
  nkt <- set.rank(nkt)
  fulls <- set.rank(fulls)
  nkt.shared <- shared.repertoire(nkt[2:4], .min.ppl = 1, .type='n0rr')  # беру все последовательности с 3p56p2Cp, 3p56p, 3p2Cp
  
  nkt.in.fulls  <- find.clonotypes(.data = c(nkt,fulls),
                                   .targets =  nkt.shared$CDR3.nucleotide.sequence,
                                   .method = "exact" ,
                                   .col.name = c("CDR3.amino.acid.sequence","Read.count","Read.proportion","Rank", "V.gene"),
                                   .target.col = "CDR3.nucleotide.sequence" ) # делаю кросс
  
  nkt.in.fulls<- nkt.in.fulls[ , order(names(nkt.in.fulls))]  # упорядочиваю колонки чтобы одна и та же инфа была рядом
  
 
  nkt.in.fulls <- nkt.in.fulls[c(grep("CDR3.nucleotide.sequence", names(nkt.in.fulls)),
                                 grep(paste0("CDR3.amino.acid.sequence.",sample_name,"_3p56p2Cp"), names(nkt.in.fulls)),
                                 grep("Rank",names(nkt.in.fulls)),
                                 grep("Read.proportion",names(nkt.in.fulls)),
                                 grep("Read.count",names(nkt.in.fulls)),
                                 grep("V.gene",names(nkt.in.fulls))
                 )]
  nkt.in.fulls <- arrange(nkt.in.fulls, nkt.in.fulls[,grepl(paste0("Rank.",sample_name,"_3p56p2Cp"),
                                                            names(nkt.in.fulls))]) # сортирую дата фрейм по рэнку в 3p56p2Cp
  nkt.in.fulls[is.na(nkt.in.fulls)] <- 0
 
  ## расщепление фенотипа
  
  ncol.3p56p2Cp <- grep("Read.count.*3p56p2Cp",names(nkt.in.fulls))
  ncol.3p56p <- grep("Read.count.*3p56p$",names(nkt.in.fulls))
  ncol.3p2cp <- grep("Read.count.*3p2сp$",names(nkt.in.fulls))
  ncol.3p <- grep("Read.count.*3p$",names(nkt.in.fulls))
  
  
  fit <- lm_fc(nkt.in.fulls)
  nkt.in.fulls$lm_full <- fit$coefficients[1]*nkt.in.fulls[,ncol.3p56p2Cp]+fit$coefficients[2]*nkt.in.fulls[,ncol.3p56p]+fit$coefficients[3]*nkt.in.fulls[,ncol.3p2cp]+fit$coefficients[4]*nkt.in.fulls[,ncol.3p]
  nkt.in.fulls$split3p56p2Cp <- fit$coefficients[1]*nkt.in.fulls[,ncol.3p56p2Cp]/nkt.in.fulls$lm_full
  nkt.in.fulls$split3p56p <- fit$coefficients[2]*nkt.in.fulls[,ncol.3p56p]/nkt.in.fulls$lm_full
  nkt.in.fulls$split3p2Cp <-fit$coefficients[3]*nkt.in.fulls[,ncol.3p2cp]/nkt.in.fulls$lm_full
  nkt.in.fulls$split3p <- fit$coefficients[4]*nkt.in.fulls[,ncol.3p]/nkt.in.fulls$lm_full
  ## сохраняем в красивом виде 
  
  wb <- createWorkbook()
  addWorksheet(wb,"Intersections")
  hs1 <- createStyle(fgFill = "#DCE6F1", halign = "CENTER", textDecoration = "bold", border = "TopBottom")
  
 
 
  s <- createStyle(numFmt = "# ###",halign = "CENTER")
  addStyle(wb,1,style =s, cols=c(grep("Read.count",names(nkt.in.fulls)),grep("Rank",names(nkt.in.fulls))), rows=2: nrow(nkt.in.fulls),gridExpand = T)
  
  s <- createStyle(border = "Left",borderStyle = "thick")
  addStyle(wb,1,style =s,cols=c(min(grep("Read.count",names(nkt.in.fulls))),
                                min(grep("Rank",names(nkt.in.fulls))),
                                min(grep("Read.proportion",names(nkt.in.fulls))),
                                min(grep("V.gene",names(nkt.in.fulls)))),
           rows=2: nrow(nkt.in.fulls),
           gridExpand = T,
           stack = T)
  
  s <- createStyle(border = "Right",borderStyle = "thick")
  addStyle(wb,1,style =s,cols=c(max(grep("Read.count",names(nkt.in.fulls))),
                                max(grep("Rank",names(nkt.in.fulls))),
                                max(grep("Read.proportion",names(nkt.in.fulls))),
                                max(grep("V.gene",names(nkt.in.fulls)))),
           rows=2: nrow(nkt.in.fulls),
           gridExpand = T,
           stack = T)
                             
  
  s <- createStyle(numFmt = "#0.###%",halign = "CENTER")
  addStyle(wb,1,style =s, cols=c(grep("Read.proportion",names(nkt.in.fulls)),grep("split",names(nkt.in.fulls))), rows=2: nrow(nkt.in.fulls),gridExpand = T, stack =  T)
  
  writeData(wb,1,nkt.in.fulls, headerStyle = hs1)
  setColWidths(wb,1,cols=1,widths = 57)
  setColWidths(wb,1,cols=2,widths = 25)
  setColWidths(wb,1,cols=c(3:14),widths = 19)
  
  
  
  
  
  
  saveWorkbook(wb, paste0("/home/artem/Desktop/IBH/NKT/new_nkt/analysis/",sample_name,".xlsx"), overwrite = TRUE)
  
  
  
  return(nkt.in.fulls) 
}

# with V.gene
# nkt.intersect.V <- function(nkt, fulls,sample_name) {
#   nkt <- set.rank(nkt)
#   fulls <- set.rank(fulls)
#   nkt.shared <- shared.repertoire(nkt[2:4], .min.ppl = 1, .type='nvrr')  # беру все последовательности с 3p56p2Cp, 3p56p, 3p2Cp
#   
#   nkt.in.fulls  <- find.clonotypes(.data = c(nkt,fulls), 
#                                    .targets =  data.frame(cbind(CDR3.nucleotide.sequence = nkt.shared$CDR3.nucleotide.sequence, 
#                                                                 V.gene = nkt.shared$V.gene)), 
#                                    .method = "exact",
#                                    .col.name = c("CDR3.amino.acid.sequence","Read.count","Read.proportion","Rank"),
#                                    .target.col = c("CDR3.nucleotide.sequence", "V.gene")) # делаю кросс
#   
#   nkt.in.fulls<- nkt.in.fulls[ , order(names(nkt.in.fulls))]  # упорядочиваю колонки чтобы одна и та же инфа была рядом
#   
#   
#   nkt.in.fulls <- nkt.in.fulls[,c(2,1,ncol(nkt.in.fulls),3:(ncol(nkt.in.fulls)-1))]
#   nkt.in.fulls <- arrange(nkt.in.fulls, nkt.in.fulls[,grepl(paste0("Rank.",sample_name,"_3p56p2Cp"),
#                                                            names(nkt.in.fulls))]) # сортирую дата фрейм по рэнку в 3p56p2Cp
#   nkt.in.fulls[is.na(nkt.in.fulls)] <- 0
#   return(nkt.in.fulls) 
#   
# }
  

## Построение линейной модели для расчета расщепления фенотипа 

  lm_fc <- function(cloneset){
    reads.3p56p2Cp <- cloneset[,grep("Read.count.*3p56p2Cp",names(cloneset))]
    reads.3p56p <- cloneset[,grep("Read.count.*3p56p$",names(cloneset))]
    reads.3p2cp <- cloneset[,grep("Read.count.*3p2сp$",names(cloneset))]
    reads.3p <- cloneset[,grep("Read.count.*3p$",names(cloneset))]
    reads.f <- cloneset[,grep("Read.count.*_F$",names(cloneset))[1]]
    model <- lm(reads.f~reads.3p56p2Cp+reads.3p56p+reads.3p2cp+reads.3p+0)
    return(summary(model))
  }


## Функция случайной выборки фракции (например, из фулла)
      
  vyb <- function(from, len,repl){
  sample.frac <- as.data.frame(table(sample(from$CDR3.nucleotide.sequence,
                               len,
                               prob = from$Umi.proportion,
                               replace = repl)))
  sample.frac <- sample.frac[order(sample.frac$Freq,decreasing = T),]
  return(sample.frac)
  }
  
  
## Функция расчета вероятности вытянуть столько или больше UMI каждого клона при случайной выборке
  
  clone.count.prob <- function(full, fraction, adj.method,cl) {

    # no_cores <- detectCores() - 1
    # cl <- makeCluster(no_cores)
    # clusterExport(cl,strsplit(deparse(substitute(fraction)),'\\$')[[1]][1]  )
    # clusterExport(cl,strsplit(deparse(substitute(full)),'\\$')[[1]][1] )
    
    plist <- parApply(fraction,1, cl=cl,function(clone){
      m = full[full$CDR3.nucleotide.sequence == as.character(clone[5]),]$Umi.count # UMI.Count этого клона в фулле (белые шары в коробке)
      if (length(m)==0) m=1
      n = sum(full$Umi.count)-m # суммарный UMI.Count всех остальных клонов в фулле (черные шары)
      k = sum(fraction$Umi.count) # общее количество UMI в этой фракции (кол-во вытягиваемых шаров)
      x = as.integer(clone[1])  # для скольки вытянутых белых шаров считаем вероятность
      p <-  phyper(x,m,n,k,lower.tail = FALSE) # вероятность вытянуть x или более белых шаров
    
      return(p)

    })
    # stopCluster(cl)
    
    fraction$p <- plist
    fraction$p.adj <- p.adjust(p=plist, method = adj.method,n =  length(plist)) # поправка на множественные сравнения
    
    
    colnames(fraction)[ncol(fraction)] <- paste0("p.adj (", adj.method, ")")
    fraction <- fraction[,c(1:4,ncol(fraction),5:(ncol(fraction)-1) )]
    
    return(fraction)
  }


  
 
  
 # YB.CMV.p <- lapply(YB.CMV.nkt,function(x){ clone.count.prob(full = YB$YB_NKT_F, fraction = x, adj.method = "fdr",cl)} )
  