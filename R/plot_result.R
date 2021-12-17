#' Plots the cutsite probabilities (lambda parameter) for EM results
#'
#'@param lambda a vector of lambda probabilities derived from the multinomial distribution
#'@param name the name of the run for the title of the plot
#'@param out the output file to save plot to
#'@examples
#'plot_cutsites(result$lambda,"my experiment","my_experiment_cutsites.png")
#'@export

plot_cutsites <- function(lambda,name,out){
  idx <- length(lambda)/2
  f_lambda <- lambda[1:idx]
  r_lambda <- lambda[idx+1:length(lambda)]
  r_lambda <- r_lambda[1:idx]
  motif_start <- 101
  motif_end <- length(lambda)/2 - 200 + motif_start

  png(file=out, width = 600, height = 400)
  plot(r_lambda,type="l",col="red", xaxt="n", xlab = "Dist. to motif(bp)", ylab="Cut-site probability",main=paste("Cutsite Probabilities: ",name, sep=" "))
  lines(f_lambda, col="blue")
  abline(v=c(motif_start,motif_end), col=c("black", "black"), lty=c(2,2), lwd=c(1, 1))
  axis(1,at=c(0,50,motif_start,motif_end,motif_end+50,motif_end+100),labels = c(-100,-50,0,0,50,100))
  legend(1,95,legend=c("Reverse strand","Forward strand"),col=c("red","blue"))
  dev.off()


}

#' Calculate statistcs between EM results and ChIP-seq data
#'
#' @param results a binary vectory of bound status predictions
#' @param region_labels dataframe of coordinates corresponding to results with columns 'chr','start','end' and rows as the coordinates in "chr:start-end" format
#' @param chip_file a bed file of ChIP-seq peaks
#' @return A list of statistics including precision, FDR, FPR, sensitivity,and a dataframe labeled_df which is results with added columns for predicted labels and ChIP-seq labels
#'
#' @import GenomicRanges
#'@export
compare_results<-function(results,region_labels, chip_file){
  #load in our results and label them with what regions they came from
  results = cbind(region_labels,results)
  colnames(results) = c("chr","start","end","bound_status")
  bound_results = results[results$bound_status == 1,]
  unbound_results = results[results$bound_status == 0,]
  chip_df = read.table(chip_file)[,c(1:3)]
  colnames(chip_df) = c("chr","start","end")
  #drop duplicate rows
  chip_df = chip_df[!duplicated(chip_df),]
  rownames(chip_df) = paste(chip_df$chr,":",chip_df$start,"-",chip_df$end, sep="")

  #only consider chip-seq peaks that overlap our regions as true positives and visa versa
  motif_gr = as(region_labels,"GRanges")
  chip_gr = as(chip_df, "GRanges")
  chip_pos_idx = unique(queryHits(findOverlaps(motif_gr, chip_gr))) # indeces in our regions that should be positive
  chip_neg = region_labels[-chip_pos_idx,] #should not be bound according to chip-seq
  chip_pos = region_labels[chip_pos_idx,] #should be bound according to chip-seq

  #construct a df with a predictions and labels columns for roc plot
  label = rep(0,dim(region_labels)[1])
  label[chip_pos_idx] = 1
  results = cbind(results,label)

  #calculate metrics based on just the above definition
  n_true_pos = sum(rownames(chip_pos) %in% rownames(bound_results))
  n_false_pos = dim(bound_results)[1] - n_true_pos

  n_true_neg = sum(rownames(chip_neg) %in% rownames(unbound_results))
  n_false_neg = dim(unbound_results)[1] - n_true_neg

  #calculate sensitivity
  sensitivity = n_true_pos/(n_true_pos+n_false_neg)

  #calculate the false positive rate
  FPR = n_false_pos/(n_false_pos + n_true_neg)

  #calculate the FDR
  FDR = n_false_pos/(n_true_pos+n_false_pos)

  #get precision
  precision = n_true_pos/(n_true_pos+n_false_pos)
  return(list(precision=precision,FDR=FDR,FPR=FPR,sensitivity=sensitivity,labeled_df=results))
}
#' Plots ROC curve in precrec's style
#'
#' @param labelled_df a dataframe outputted by compare_results or with columns for bound_status (predicted labels) and label
#' @import precrec
#' @import ggplot2
#'
plot_roc<- function(labelled_df){
  precrec_obj <- evalmod(scores=labelled_df$bound_status,labels=labelled_df$label)
  autoplot(precrec_obj)
}
#' Plots ROC curve in PRROC's style
#'
#'@param labelled_df a dataframe outputted by compare_results or with columns for bound_status (predicted labels) and label
#'@import PRROC
#'@import ggplot2
#'
plot_roc2 <- function(labelled_df){
  PRROC_obj <- roc.curve(scores.class0 = labelled_df$bound_status, weights.class0=labelled_df$label,
                         curve=TRUE)
  plot(PRROC_obj)

}
