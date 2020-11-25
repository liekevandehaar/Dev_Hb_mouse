##### determine correlations and signficance #####
# author: Maria Tosches
# edited by: Juliska E Boer
# date: 03 Nov 2020

SpPermute = function(ExpressionTable1, DEgenes1, ExpressionTable2, DEgenes2, nPermutations = 1000, genes.use='intersect', corr.method = 'spearman'){
  
  #Step1: Take intersect of DEgenes for analysis
  DEgenes = intersect(DEgenesSpecies1,DEgenesSpecies2)
  
  #Step2: Prune Expression Tables & Remove rows with no expression
  Sp1 = ExpressionTable1[rownames(ExpressionTable1) %in% DEgenes1,]
  Sp1 = Sp1[rowSums (Sp1)!=0,]
  Sp2 = ExpressionTable2[rownames(ExpressionTable2) %in% DEgenes2,]
  Sp2 = Sp2[rowSums (Sp2)!=0,]
  
  #Step3: Scale Expression Tables by gene expression sum
  Sp1 <- sweep(Sp1,MARGIN=1,FUN="/",STATS=rowSums(Sp1))
  Sp2 <- sweep(Sp2,MARGIN=1,FUN="/",STATS=rowSums(Sp2))
  
  #Step4: Merge Expression Tables
  geTable = merge(Sp1,Sp2, by='row.names', all=F)
  rownames(geTable) = geTable$Row.names
  geTable = geTable[,2:ncol(geTable)]
  
  #Step5:  Correlation
  #5a:  Correlation
  Corr.Coeff.Table = cor(geTable,method=corr.method)
  
  #5b:  Shuffle data
  shuffled.cor.list = list()
  pb   <- txtProgressBar(1, 100, style=3)
  
  for (i in 1:nPermutations){
    shuffled = apply(geTable[,1:ncol(Sp1)],1,sample)
    shuffled2 = apply(geTable[,(ncol(Sp1)+1):ncol(geTable)],1,sample)
    shuffled = cbind(t(shuffled),t(shuffled2))
    shuffled.cor = cor(shuffled,method=corr.method)
    shuffled.cor.list[[i]] = shuffled.cor
    rm(list=c('shuffled','shuffled2','shuffled.cor'))
    if ((i %% 100) ==0){
      setTxtProgressBar(pb, (i*100)/nPermutations)
    }
  }
  
  p.value.table = matrix(ncol=ncol(geTable), nrow = ncol(geTable))
  rownames(p.value.table) = colnames(geTable)
  colnames(p.value.table) = colnames(geTable)
  
  shuffled.mean.table = matrix(ncol=ncol(geTable), nrow = ncol(geTable))
  rownames(shuffled.mean.table) = colnames(geTable)
  colnames(shuffled.mean.table) = colnames(geTable)
  
  a = combn(1:ncol(geTable),2)
  for (i in 1:ncol(a)){
    cor.scores = sapply(shuffled.cor.list,"[",a[1,i],a[2,i])
    shuffled.mean.table[a[1,i],a[2,i]] = mean(cor.scores)
    shuffled.mean.table[a[2,i],a[1,i]] = mean(cor.scores)
    p.value = mean(abs(cor.scores)>=abs(Corr.Coeff.Table[a[1,i],a[2,i]]))
    p.value.table[a[1,i],a[2,i]] = p.value
    p.value.table[a[2,i],a[1,i]] = p.value
    rm(list=c('cor.scores','p.value'))
    setTxtProgressBar(pb, (i*100)/ncol(a))
  }
  
  #Step6: Return variables
  
  list.to.return = list(Corr.Coeff.Table,shuffled.mean.table,p.value.table,DEgenes,rownames(geTable),length(DEgenes),nDESp1,nDESp2)
  names(list.to.return) = c('corr.coeff','shuffled_correlation_score_means','p.value','DEgenes_intersect','DEgenes_in_analysis','nDEgenes_in_analysis','nDEgenes_Sp1','nDEgenes_Sp2')
  return(list.to.return)
  
}