library(ROCR)

###Auxiliary function to color the dots


string.to.colors = function (string, colors = NULL) 
{
  if (is.factor(string)) {
    string = as.character(string)
  }
  if (!is.null(colors)) {
    if (length(colors) != length(unique(string))) {
      (break)("The number of colors must be equal to the number of unique elements.")
    }
    else {
      conv = cbind(unique(string), colors)
    }
  }
  else {
    conv = cbind(unique(string), rainbow(length(unique(string))))
  }
  unlist(lapply(string, FUN = function(x) {
    conv[which(conv[, 1] == x), 2]
  }))
}



color_convertion=function(x,max_scale=NULL,min_scale=NULL) {
  f <- colorRamp(c("white","red"))
  x=as.numeric(x)
  if (is.null(max_scale)) {
    max_scale=quantile(x,0.99,na.rm = T)
  }
  if (is.null(min_scale)) {
    min_scale=quantile(x,0.01,na.rm = T)
  }
  
  x_prime=ifelse(x>max_scale,max_scale,x)
  x_prime=ifelse(x<min_scale,min_scale,x)
  
  x_prime=(x_prime-min_scale)/(max_scale-min_scale)
  x_color=f(x_prime)/255
  x_color[!complete.cases(x_color),]=c(0,0,0)
  x_color=rgb(x_color)
  return(x_color)
}


##Base function to create a logistic classifier for each gene of interest

logistic_classifier=function(x,labels) {
  x=as.numeric(x)
  m=glm(labels~x,family = "binomial")
  pred=prediction(m$fitted.values, labels, label.ordering = NULL)
  auc=performance(pred, measure = "auc")
  auc=auc@y.values[[1]]
  reg_coef=coef(m)[2]
  results=c(auc,reg_coef)
  names(results)=c("AUC","Coef")
  return(results)
}

##Function to compute the AUC for all genes

Compute_AUC_marker = function(Expression_matrix,Clustering) {
  cat(paste("Starting identification of the markers using initially",ncol(Expression_matrix),"genes \n"))
  if (length(unique(as.character(Clustering)))>2) {
    stop("Please provide a clustering vector with only two groups")
  }
  if (!is.logical(Clustering)) {
    warning("The provided vector is not boolean. \n The less common group will be used as the target one.")
    Counting_cluster = table(Clustering)
    Clustering = Clustering==(names(which.max(Counting_cluster)))
  }
  
  AUC_gene=apply(Expression_matrix,MARGIN = 2,FUN = function(x) {logistic_classifier(x,Clustering)})
  AUC_gene = as.data.frame(t(AUC_gene))
  rownames(AUC_gene)=colnames(Expression_matrix)
  return(AUC_gene)
}

##Functions to visualize the AUC curve 
Plot_ROC = function(Expression_matrix,gene,Clustering) {
  if (length(unique(as.character(Clustering)))>2) {
    stop("Please provide a clustering vector with only two groups")
  }
  if (!is.logical(Clustering)) {
    warning("The provided vector is not boolean. \n The less common group will be used as the target one.")
    Counting_cluster = table(Clustering)
    Clustering = Clustering==(names(which.max(Counting_cluster)))
  }
  
  x=as.numeric(Expression_matrix[,gene])
  m=glm(labels~x,family = "binomial")
  pred=prediction(m$fitted.values, Clustering, label.ordering = NULL)
  
  tpr=performance(pred, measure = c('tpr'))
  tpr = tpr@y.values[[1]]
  
  fpr=performance(pred, measure = c('fpr'))
  fpr = fpr@y.values[[1]]
  par(las=1)
  plot(1-fpr,1-tpr,xaxs="i",yaxs="i",xlim=c(0,1),ylim=c(0,1),type="l",lwd=2,xlab="False positive rate",
       ylab="True positive rate",main=gene,cex.main = 1.3,cex=1.5)
  abline(0,1,lwd=2,lty=2,col="red4")
}

##Functions to visualize which genes are good candidates

Diagnostic_plot = function(AUC_vector,Mean_expression_vector,AUC_treshold = 0.75,Mean_expression_threshold = 20) {
  
  Mean_expression_vector = Mean_expression_vector[AUC_vector>0]
  AUC_vector = AUC_vector[AUC_vector>0]
  
  High_quality_genes = names(which(Mean_expression_vector < Mean_expression_threshold & AUC_vector>AUC_treshold))
  par(las=1)
  if (length(High_quality_genes)<1) {
    cat("No marker genes were found to fulfill the criterion ")
    plot(AUC_vector,Mean_expression_vector,xlim=c(0.5,1),xaxs="i",yaxs="i",
         xlab="AUC score",ylab="Mean expression",cex.lab=1.3,cex=1.3,
         pch=21,bg="grey")
    
  }
  if (length(High_quality_genes) >= 1) {
    plot(AUC_vector,Mean_expression_vector,xlim=c(0.5,1),xaxs="i",yaxs="i",
         xlab="AUC score",ylab="Mean expression",cex.lab=1.3,cex=1.3,ylim=c(0,1.2*max(Mean_expression_vector)),
         pch=21,bg=string.to.colors(names(AUC_vector)%in%High_quality_genes,c("orange","grey")))
    text(AUC_vector[High_quality_genes],Mean_expression_vector[High_quality_genes],High_quality_genes)
    abline(v=AUC_treshold,lwd=2,lty=2,col="grey")
    abline(h=Mean_expression_threshold,lwd=2,lty=2,col="grey")
    
  }
  
}


