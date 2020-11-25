##### plot confusion matrix holding classficiations performed by RandomForest classifier #####
# author: Juliska E Boer
# date: 03 Nov 2020

plotConfusionMatrix = function(X,row.scale=TRUE, col.scale=FALSE, max.size=5, ylab.use="Known", xlab.use="Predicted"){
  if (!col.scale & row.scale){ X = t(scale(t(X), center=FALSE, scale=rowSums(X)));  X=X*100 }
  if (col.scale & !row.scale){ X = scale(X, center=FALSE, scale=colSums(X)); X = X*100 }
  if(col.scale & row.scale){
    print("Only one of row.scale or col.scale should be true. performing row scaling by default")
    X = t(scale(t(X), center=FALSE, scale=rowSums(X)))
    X=X*100
  }
  X = X[rev(1:dim(X)[1]),]
  X = melt(X, id.vars="Known")
  colnames(X) = c("Known", "Predicted", "Percentage")
  X$Known = factor(X$Known, levels=rev(unique(X$Known)));
  X$Predicted = as.factor(X$Predicted)
  p = ggplot(X, aes(y = Known, x = Predicted)) + geom_point(aes(colour = Percentage,  size = Percentage)) + 
    scale_color_gradient(low ="lightgrey", high = "red", limits=c(0, 100)) + scale_size(range = c(1, max.size), labels = c(0,25,50,75,100)) + theme_bw() #+nogrid
  p = p + xlab(xlab.use) + ylab(ylab.use) + theme(axis.text.x=element_text(size=12, face="italic", hjust=1,angle=45)) + 
    theme(axis.text.y=element_text(size=12, face="italic"))
  print(p)
}