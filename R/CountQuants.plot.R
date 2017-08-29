#' Plot the counterfactual quantiles
#' @param CountQuants The counterfactual quantiles matrix
#' @param SDMat The standard errors matrix

#' @export

CountQuants.plot=function(CountQuants,SDMat){
  library(ggplot2)
  library(reshape2)
  cqs=melt(CountQuants)
  sd=melt(SDMat)
  sqs=merge(cqs,sd,by=c("Var1","Var2"))
  colnames(sqs) <- c("quants","ctvals","cqs","sd")
  CountQuants.plot=ggplot2::ggplot(data=sqs, aes(ctvals, cqs, ymax=cqs+1.96*sd,
                                                 ymin=cqs-1.96*sd)) +
    ggplot2::geom_line(aes(ctvals, cqs, group=quants, color=quants)) +
    # ggplot2::geom_errorbar(size=.3, width=.02) +
    ggplot2::geom_line(aes(ctvals, cqs+1.96*sd,group=quants, color=quants), linetype="dashed") +
    ggplot2::geom_line(aes(ctvals, cqs-1.96*sd,group=quants, color=quants), linetype="dashed") +
    ggplot2::geom_point(aes(ctvals, cqs, group=quants,color=quants))  +
    ggplot2::scale_y_continuous("Unconditional Counterfactual Quantiles") + ##, limits=c(-.4, .4)) +
    ggplot2::scale_x_continuous("Treatments") +   ##limits=c(0,1), breaks=c(.1,.3,.5,.7,.9)) +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.border = element_rect(colour = 'black', size=1,
                                               fill=NA,
                                               linetype='solid'))
  CountQuants.plot
}
