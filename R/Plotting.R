# Purpose : Functions for plotting estimated densities and the local fdr

#' Plot Estimated Densities
#' 
#' @param z Observed z scores
#' @param f.mix Mixture density, optional
#' @param f.en Empirical null density, optional
#' @param k Bins for histogram
#' @import ggplot2 
#' @export

plotDens = function(z,f.mix,f.en,k=120){
  # Range
  r = Range(z);
  # Evaluation grid
  x = seq(from=r[1],to=r[2],length.out=1001);
  # Plot
  y = NULL;
  q = ggplot() + geom_histogram(data=data.frame(z),aes(x=z,y=..density..),alpha=0.8,fill="gray",bins=k);
  q = q + theme_bw() + xlab("Z Scores") + ylab("Density");
  if(!missing(f.mix) & !missing(f.en)){
    D.mix = data.frame("x"=x,"y"=f.mix(x));
    D.en = data.frame("x"=x,"y"=f.en(x));
    q = q + geom_line(data=D.mix,aes(x=x,y=y,color="Mixture",linetype="Mixture"));
    q = q + geom_line(data=D.en,aes(x=x,y=y,color="eNull",linetype="eNull"));
    q = suppressMessages(q + scale_color_manual(name="Curve",values=c("coral","royalblue")));
    q = suppressMessages(q + scale_linetype_manual(name="Curve",values=c("dotted","solid")));
  } else if (!missing(f.mix)){
    D.mix = data.frame("x"=x,"y"=f.mix(x));
    q = q + geom_line(data=D.mix,aes(x=x,y=y,color="Mixture",linetype="Mixture"));
    q = suppressMessages(q + scale_color_manual(name="Curve",values=c("royalblue")));
    q = suppressMessages(q + scale_linetype_manual(name="Curve",values=c("solid")));
  } else if (!missing(f.en)){
    D.en = data.frame("x"=x,"y"=f.en(x));
    q = q + geom_line(data=D.en,aes(x=x,y=y,color="eNull",linetype="eNull"));
    q = suppressMessages(q + scale_color_manual(name="Curve",values=c("coral")));
    q = suppressMessages(q + scale_linetype_manual(name="Curve",values=c("dotted")));
  }
  return(q);
}

#' Plot False Discovery Rates
#' 
#' @param z Observed z scores. 
#' @param lfdr Local fdr function
#' @param tFDR Tail fdr function
#' @param k Bins for histogram
#' 
#' @import ggplot2
#' @importFrom graphics hist
#' @export

plotFDR = function(z,lfdr,tFDR,k=120){
  # Estimate histogram
  D = hist(x=z,breaks=k,plot=F);
  # Rescaled histogram
  D = data.frame("x"=D$mids,"y"=D$density/max(D$density));
  # Evaluation grid
  x = seq(from=min(D$x),to=max(D$x),length.out=201);
  # Base Plot
  y=NULL;
  q = ggplot() + geom_bar(data=D,aes(x=x,y=y),stat="identity",fill="gray",color="gray",alpha=0.8);
  q = q + theme_bw() + xlab("Z Scores") + ylab("Scaled Density");
  q = q + theme(axis.title=element_text(size=12),title=element_text(size=14));
  # False discovery rate
  if(!missing(lfdr) & !missing(tFDR)){
    D.lfdr = data.frame("x"=x,"y"=lfdr(x));
    D.tFDR = data.frame("x"=x,"y"=tFDR(x));
    q = q + geom_line(data=D.lfdr,aes(x=x,y=y,color="lfdr",linetype="lfdr"));
    q = q + geom_line(data=D.tFDR,aes(x=x,y=y,color="tFDR",linetype="tFDR"));
    q = suppressMessages(q + scale_color_manual(name="Curve",values=c("darkseagreen","goldenrod")));
    q = suppressMessages(q + scale_linetype_manual(name="Curve",values=c("dashed","dashed")));
  } else if (!missing(lfdr)){
    D.lfdr = data.frame("x"=x,"y"=lfdr(x));
    q = q + geom_line(data=D.lfdr,aes(x=x,y=y,color="lfdr",linetype="lfdr"));
    q = suppressMessages(q + scale_color_manual(name="Curve",values=c("darkseagreen")));
    q = suppressMessages(q + scale_linetype_manual(name="Curve",values=c("dashed")));
  } else if (!missing(tFDR)){
    D.tFDR = data.frame("x"=x,"y"=tFDR(x));
    q = q + geom_line(data=D.tFDR,aes(x=x,y=y,color="tFDR",linetype="tFDR"));
    q = suppressMessages(q + scale_color_manual(name="Curve",values=c("goldenrod")));
    q = suppressMessages(q + scale_linetype_manual(name="Curve",values=c("dashed")));
  }
  # Output
  return(q);
}