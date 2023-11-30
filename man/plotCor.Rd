\name{plotCor}
\alias{plotCor}
\title{
Plot the correlation
}
\description{
Plot the correlation of a predicted target
}
\usage{
plotCor(DN,GN,Pred,DRS,GES,plot=TRUE)
}
\arguments{
    \item{DN}{Drug Name}
    \item{GN}{Gene Name}
    \item{Pred}{Output from prediction object}
    \item{DRS}{Drug response score}
    \item{GES}{Gene Effect scores}
    \item{plot}{ default is plot=TRUE}
}
\value{Correlation plot}
\author{
sanjusinha7, Trinh Nguyen
}
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
plotCor (DOI,GOI,Pred.d,sec.prism,KO.GES)
}

