##### determine differentially expressed genes using MAST #####
# author: Maria Tosches
# edited by: Juliska E Boer
# date: 03 Nov 2020

DE_Gene_Union = function(SeuratObject,clusters,min.pct=0.4,thresh.use=log(2), data.info.col){
  DEgenes = list(character(),character(),character(),character())
  names(DEgenes) = c('union_fdr<0.05','union_fdr<0.01','union_fdr<0.005','union_fdr<1e-5')
  combinations = t(combn(as.character(clusters),2)) 
  colnames(combinations) = c("ident.1","ident.2")
  nDEgenes = matrix(nrow=nrow(combinations),ncol=4)
  colnames(nDEgenes) = c('union_fdr<0.05','union_fdr<0.01','union_fdr<0.005','union_fdr<1e-5')
  rownames(nDEgenes) = paste('cluster',combinations[,1],'cluster',combinations[,2],sep='_')
  pb   <- txtProgressBar(1, 100, style=3)
  for (i in 1:nrow(combinations)){
    a = combinations[i,1]
    b = combinations[i,2]
    suppressMessages(diff.genes <-  FindMarkers.MAST(SeuratObject,a,b,min.pct,thresh.use,data.info.col))
    #matrix from FindMarkers.MAST only contains genes that have a fdr<0.05
    DEgenes[[paste("cluster",a,"cluster",b, sep="_")]] = diff.genes
    DEgenes[['union_fdr<0.05']] = union(DEgenes[['union_fdr<0.05']],rownames(diff.genes))
    DEgenes[['union_fdr<0.01']] = union(DEgenes[['union_fdr<0.01']],rownames(diff.genes[diff.genes$fdr<0.01,]))
    DEgenes[['union_fdr<0.005']] = union(DEgenes[['union_fdr<0.005']],rownames(diff.genes[diff.genes$fdr<0.005,]))
    DEgenes[['union_fdr<1e-5']] = union(DEgenes[['union_fdr<1e-5']],rownames(diff.genes[diff.genes$fdr<1e-5,]))
    nDEgenes[i,1] = nrow(diff.genes)
    nDEgenes[i,2] = nrow(diff.genes[diff.genes$fdr<0.01,])
    nDEgenes[i,3] = nrow(diff.genes[diff.genes$fdr<0.005,])
    nDEgenes[i,4] = nrow(diff.genes[diff.genes$fdr<1e-5,])
    rm(list = c('diff.genes','a','b'))
    setTxtProgressBar(pb, (i*100)/nrow(combinations))
  }
  DEgenes[['nDEgenes']] = nDEgenes
  DEgenes[['nDEgene_union']] = sapply(DEgenes[1:4],length)
  return(DEgenes)
}

FindMarkers.MAST <- function(object, id1, id2, min.pct, thresh.use,data.info.col){
  cells.1 <- names(object@active.ident[object@active.ident == id1])
  cells.2 <- names(object@active.ident[object@active.ident == id2])
  cells.to.compare <- c(cells.1,cells.2)
  
  raw.neur.counts <- as.matrix(object@assays$RNA@counts[,rownames(object[[]])])
  neur.alldata <- log(raw.neur.counts + 1)/log(2) # i.e. log2 base
  
  genes.use = rownames(object@assays$RNA@data)
  thresh.min = 0
  data.temp1=round(apply(neur.alldata[genes.use, cells.1, drop = F],1,function(x)return(length(x[x>thresh.min])/length(x))),3)
  data.temp2=round(apply(neur.alldata[genes.use, cells.2, drop = F],1,function(x)return(length(x[x>thresh.min])/length(x))),3)
  data.alpha=cbind(data.temp1,data.temp2); colnames(data.alpha)=c("pct.1","pct.2")
  alpha.min=apply(data.alpha,1,max)
  names(alpha.min)=rownames(data.alpha)
  genes.use=names(which(alpha.min>min.pct))
  
  neur.data <- neur.alldata[genes.use, cells.to.compare]
  
  symbolid <- rownames(neur.data)
  primerid <- symbolid
  rownames(neur.data) <- primerid
  neur.cData <- as.data.frame(cbind(as.character(object@active.ident), object@meta.data[,data.info.col]))
  colnames(neur.cData) <- c("cluster",colnames(object@meta.data)[data.info.col])
  rownames(neur.cData) <- names(object@active.ident)
  neur.cData$cluster <- as.character(neur.cData$cluster) 
  neur.cData <- neur.cData[cells.to.compare,]
  neur.cData$wellKey <- as.character(cells.to.compare)
  
  neur.fData <- cbind(primerid,symbolid)
  neur.fData <- as.data.frame(neur.fData)
  neur.fData[,1] <- as.character(neur.fData[,1]) 
  
  neur.to.compare <- FromMatrix(neur.data, neur.cData, neur.fData)
  
  colData(neur.to.compare)$cngeneson <- scale(colSums(assay(neur.to.compare)>0))
  colData(neur.to.compare)$cluster <- as.factor(colData(neur.to.compare)$cluster)
  colData(neur.to.compare)$cluster <- relevel(colData(neur.to.compare)$cluster, as.character(id1))
  
  zlmCond <- zlm(~ cluster + cngeneson, neur.to.compare)
  summaryCond <- summary(zlmCond, doLRT=paste0("cluster", id2)) 
  
  summaryDt <- summaryCond$datatable  
  summaryDt2 <- summaryDt	
  set1 <- summaryDt2[(contrast==paste0("cluster", id2) & component == "H"),c(1,4)]
  colnames(set1) <- c("primerid","Pr")
  set2 <- summaryDt2[(contrast==paste0("cluster", id2) & component == "logFC"), .(primerid, coef, ci.hi, ci.lo)]
  
  fcHurdle <- merge(set1,set2, by="primerid")
  fcHurdle[,fdr:=p.adjust(Pr, "fdr")]
  
  expMean=function(x) {
    return(log(mean(exp(x)-1)+1))
  }  
  
  data.avg.diff <- as.matrix(object@assays$RNA@data[rownames(object@assays$RNA@data) %in% genes.use, colnames(object@assays$RNA@data) %in% cells.to.compare])
  data1 <- data.avg.diff[,cells.1]
  data2 <- data.avg.diff[,cells.2]
  genes.use2<-rownames(data.avg.diff)
  avg_diff=unlist(lapply(genes.use2,function(x)(expMean(as.numeric(data1[x,]))-expMean(as.numeric(data2[x,])))))
  names(avg_diff)<-rownames(data.avg.diff)
  
  fcHurdle2 <- as.data.frame(fcHurdle)
  rownames(fcHurdle2)<-fcHurdle2$primerid
  avg_diff <- avg_diff[rownames(fcHurdle2)]
  
  fcHurdle3 <- cbind(fcHurdle2, avg_diff)
  
  fcHurdleSig <- fcHurdle3[fcHurdle3$fdr<0.05 & abs(fcHurdle3$avg_diff)>thresh.use,]
  fcHurdleSig <- fcHurdleSig[complete.cases(fcHurdleSig),]
  fcHurdleSig <- cbind(fcHurdleSig[,c(2:3,6:7)],data.alpha[row.names(fcHurdleSig),])
  fcHurdleSig <- fcHurdleSig[order(abs(fcHurdleSig$avg_diff),decreasing=T),]
  
  return(fcHurdleSig)
  
}