\name{DTR}
\alias{DTR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of drug target response for primary or secondary targets
}
\description{
Generates a plot of drug target response for primary or secondary targets based on a drug of interest and a gene of interest.
}
\usage{
DTR(DN,GN,Pred,Exp,DRS,GES,CutOff= 3,plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DN}{
Drug of interest
}
  \item{GN}{
Gene of interest
}
  \item{Pred}{
Prediction object, an output result from prediction whether it is a primary target or secondary target
}
  \item{Exp}{
Expression matrix
}
  \item{DRS}{
Drug response matrix
}
  \item{GES}{
Gene Effect Scores
}
  \item{plot}{
whether users want to plot, default is true
}
  \item{CutOff}{
cutoff value for gene expression of gene of interest high or low
}
}

\value{
vialbility after KO vs drug response of gene of interest low vs high
}

\author{
 sanjusinha7
}

\examples{
library(parallel)
data (OntargetM)
library(parallel)
set.seed (12345)
S.Drugs <- c('K70301465','K09951645')
KO.GES <- OntargetM$avana_CRISPR
sec.prism <- OntargetM$secondary_prism
d.expr <- OntargetM$expression_20Q4
sim <- mclapply(S.Drugs,function(x) GetSim(x,sec.prism,KO.GES),mc.cores = 2)
names(sim) <- S.Drugs
Meta.data <- OntargetM$DrugMetadata
DrugTargetSim <- PredTarget(sim,Meta.data)
Drug.Gene.max.sim <- PredMaxSim(sim,Meta.data)
identical ( DrugTargetSim[,1],Drug.Gene.max.sim[,1] )
Pred.d <-cbind (DrugTargetSim,Drug.Gene.max.sim )
DOI = 'ibrutinib'
GOI ='BTK'
DTR(DOI,GOI,Pred.d,d.expr,sec.prism,KO.GES,CutOff= 2)
}

