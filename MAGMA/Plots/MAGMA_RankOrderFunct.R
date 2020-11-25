##### function for calculating the rank-order correlations between two performed MAGMA protocols #####
# author: Juliska E Boer
# date: 04 Nov 2020

calc_rank_cor <- function(kyoko, skene, all){
  kyoko$R = rank(kyoko$P)
  skene$R = rank(skene$P)
  if (all==TRUE){
    colnames(kyoko) = c("VARIABLE", "P-kyoko", "Dataset", "Rkyoko")
    colnames(skene) = c("VARIABLE", "p-skene", "Rskene")
  }
  else {
    colnames(kyoko) = c("VARIABLE", "P-kyoko", "Rkyoko")
    colnames(skene) = c("VARIABLE", "p-skene", "Rskene")
  }
  both = full_join(kyoko, skene, by="VARIABLE")
  both$d = abs(both$`Rkyoko`-both$`Rskene`)
  both$d2 = (both$d^2)
  sum_d2 = sum(both$d2)
  n = dim(both)[1]
  rs = (1 - ((6*sum_d2)/(n*((n^2) - 1))))
  to.return = list(rs, as.data.frame(both))
  return(to.return)
}