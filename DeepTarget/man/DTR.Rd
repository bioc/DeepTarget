\name{DTR}
\alias{DTR}
\title{
Predicting Drug Target Response (DTR) for primary or secondary targets
}
\description{
Predicting whether the drug is likely response to primary or secondary targets and also generates the plot for visualization.
}
\usage{
DTR(DN,GN,Pred,Exp,DRS,GES,CutOff= 3,plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{DN}{Drug of interest}
    \item{GN}{Gene of interest}
    \item{Pred}{Prediction object, an output result from prediction whether it is a primary target or secondary target}
    \item{Exp}{Expression matrix}
    \item{DRS}{Drug response matrix}
    \item{GES}{Gene Effect Scores}
    \item{plot}{whether users want to plot, default is true}
    \item{CutOff}{cutoff value for gene expression of gene of interest high or low}
}
\value{vialbility after KO vs drug response of gene of interest low vs high cut-off values set by users}
\author{sanjusinha7, Trinh Nguyen}
\examples{
library(BiocParallel)
data (OntargetM)
set.seed (12345)
S.Drugs <- c('K70301465','K09951645')
KO.GES <- OntargetM$avana_CRISPR
sec.prism <- OntargetM$secondary_prism
d.expr <- OntargetM$expression_20Q4
sim.out <- bplapply(S.Drugs,function(x) computeCor(x,sec.prism,KO.GES))
names(sim.out) <- S.Drugs
Meta.data <- OntargetM$DrugMetadata
DrugTargetSim <- PredTarget(sim.out,Meta.data)
Drug.Gene.max.sim <- PredMaxSim(sim.out,Meta.data)
identical ( DrugTargetSim[,1],Drug.Gene.max.sim[,1] )
Pred.d <-cbind (DrugTargetSim,Drug.Gene.max.sim )
DOI = 'ibrutinib'
GOI ='BTK'
DTR(DOI,GOI,Pred.d,d.expr,sec.prism,KO.GES,CutOff= 2)
}

