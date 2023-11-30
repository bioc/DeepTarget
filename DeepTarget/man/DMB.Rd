\name{DMB}
\alias{DMB}
\title{ Predicting Drug Mutant Binding for mutant or non-mutant form}
\description{
Predicting whether the drug is likely bind to mutant or non-mutant form and also generates the plot for visualization.
}
\usage{
DMB(DrugName,GOI,Pred,Mutant,DRS,GES,plot=TRUE)
}
\arguments{
    \item{DrugName}{Drug of interest}
    \item{GOI}{Gene of interest}
    \item{Pred}{Prediction object resulting from both PredTarget and PredMaxSim functions to predict whether it is a primary target or secondary target}
    \item{Mutant}{Mutant matrix}
    \item{DRS}{Drug response matrix}
    \item{GES}{Gene Effect Scores}
    \item{plot}{Default is TRUE for plotting}
}
\value{
The plot of viability after KO as the X-axis vs drug response in a mutant target as the Y-axis.
}
\author{
sanjusinha7, Trinh Nguyen
}
\examples{
library(BiocParallel)
data (OntargetM)
S.Drugs <- c('K70301465','K09951645')
KO.GES <- OntargetM$avana_CRISPR
sec.prism <- OntargetM$secondary_prism
d.mt <- OntargetM$mutations_mat
sim.out <- bplapply(S.Drugs,function(x) computeCor(x,sec.prism,KO.GES))
names(sim.out) <- S.Drugs
Meta.data <- OntargetM$DrugMetadata
DrugTargetSim <- PredTarget(sim.out,Meta.data)
Drug.Gene.max.sim <- PredMaxSim(sim.out,Meta.data)
identical ( DrugTargetSim[,1],Drug.Gene.max.sim[,1])
Pred.d <-cbind (DrugTargetSim,Drug.Gene.max.sim)
DOI = 'dabrafenib'
GOI = 'BRAF'
DMB (DOI,GOI,Pred.d,d.mt,sec.prism,KO.GES)
}

