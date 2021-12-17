#' Run CENTIPEDE model
#'
#' @param input_list List containing PWMScan, TSS data, and X
#' @param cons_file Optional matrix including conservation scores
#' @return The results from running the CENTIPEDE model (bound status)
#' @import edgeR
#' @import AnnotationDbi
#' @import org.Hs.eg.db
#' @import GenomicRanges
#' @export
run_centipede<-function(input_list,cons_file=NULL){ #Input pwm scan file, DNase data, bam... Conservation data should be optional
  get_TSSs <- function(pwm_scan_file){ 
    motif_loci <- pwm_scan_file
    motif_loci = cbind(rep("hi",nrow(motif_loci)),motif_loci)
    colnames(motif_loci) = c("name","chr","start","end","score","strand","pval","match")
    nearest = c() 
    midpoints = floor((motif_loci$end - motif_loci$start)/2) + motif_loci$start
    nearest_outs = nearestTSS(chr = motif_loci$chr, locus = midpoints) #This is optimized to search many loci at once so now its much much faster
    nearest = 1/(1+(abs(nearest_outs$distance)/1000)) 
    return(nearest)
  }
  
  TSS<-get_TSSs(input_list$regions[,-8])
  centi_X <- as.matrix(input_list$mat)
  Xlist = list(DNase = input_list$mat)
  init_pl <- as.numeric(rowSums(Xlist[[1]])> quantile(rowSums(Xlist[[1]]),0.9))

  pwm_scan = input_list$regions
  pwm_scan <- pwm_scan[,-8] #get rid of q-value column
  colnames(pwm_scan) = c("motif","chr","start","end","strand","score","pval","seq")
  if(missing(cons_file)){
    comb_priors<- data.frame(cbind(pwm_scan$score,TSS))#Only account for PWM and TSS data, in same order
    comb_priors <- as.matrix(comb_priors)
  } else {
    cons = read.table(cons_file)
    comb_priors<- data.frame(cbind(pwm_scan$score,TSS,cons$V4))#PWM, TSS, and cons data
    comb_priors <- comb_priors[-45,]
    comb_priors <- apply(comb_priors, 2, as.numeric)
    comb_priors <- comb_priors[which(comb_priors[,1]>=13.288),]
    comb_priors <- as.matrix(comb_priors)
    init_pl <- init_pl[-45]
    init_pl <- init_pl[which(comb_priors[,1]>=13.288)]
    centi_X <- centi_X[-45,]
    centi_X <- centi_X[which(comb_priors[,1]>=13.288),]
  }

  #run Dnase dataset
  #dnase_input$mat[,220] = dnase_input$mat[,221]
  dnase_EM <- CentipedeMixEm$new(centi_X, num_components = as.integer(2), prior_info = comb_priors, as.vector(init_pl))
  dnase_res <- dnase_EM$run.EM(loglik_tol = 1e-5)

  return(dnase_res)
}








