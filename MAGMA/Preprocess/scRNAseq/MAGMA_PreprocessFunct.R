##### ANOVA and gene conversion functions for preprocessing scRNAseq data for MAGMA #####
# author: Juliska E Boer
# date: 04 Nov 2020

ANOVA_infogenes <- function(object){
  exp = object@assays$RNA@data
  if(class(exp[1,1])=="character"){
    exp = as.matrix(exp)
    storage.mode(exp) <- "numeric"
  }
  clusters = as.character(object@active.ident)
  summed = apply(exp,1,sum)
  exp = exp[summed!=0,]
  mod  = model.matrix(~clusters)
  fit = lmFit(exp,mod)
  eb = eBayes(fit)
  pF = p.adjust(eb$F.p.value,method="BH")
  exp = exp[pF<0.00001,]
  DEGs = rownames(exp)
  return(DEGs)
}

MMUgenes_toHSgenes <- function(mmu){
  ncbi_mouse$Synonyms = paste0("|", ncbi_mouse$Synonyms, "|")
  mmu$mm.ensg = ncbi_mouse$ensg[match(mmu$mm.symbol, ncbi_mouse$Symbol)]
  mmu$mm.ensg[is.na(mmu$mm.ensg)] = sapply(mmu$mm.symbol[is.na(mmu$mm.ensg)], function(x){
    tmp = unique(ncbi_mouse$ensg[grepl(paste0("\\|",x,"\\|"), ncbi_mouse$Synonyms)]);
    if(length(tmp)==1){tmp}else{NA}
  })
  mmu$hs.ensg = mm2hs$hs.ensg[match(mmu$mm.ensg, mm2hs$mm.ensg)]
  dup = unique(mmu$hs.ensg[duplicated(mmu$hs.ensg)])
  mmu$hs.ensg[mmu$hs.ensg %in% dup] = NA
  return(mmu)
}