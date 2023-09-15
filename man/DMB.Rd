\name{DMB}
\alias{DMB}
\title{Plot of drug target response for mutant and non-mutant groups}
\description{
Generates a plot of Drug target response for primary or secondary targets based on a drug of interest.
}
\usage{
DMB(DrugName,GOI,Pred,Mutant,DRS,GES,plot=TRUE)
}
\arguments{
    \item{DrugName}{Drug of interest}
    \item{GOI}{Gene of interest}
    \item{Pred}{Prediction object, an output result from prediction whether it is a primary target or secondary target}
    \item{Mutant}{Mutant matrix}
    \item{DRS}{Drug response matrix}
    \item{GES}{Gene Effect Scores}
    \item{plot}{Default is TRUE for plotting}
}
\value{
vialbility after KO vs drug response in mutant.
}
\author{
sanjusinha7, Trinh Nguyen
}
\examples{
library(parallel)
data (OntargetM)
S.Drugs <- c('K70301465','K09951645')
KO.GES <- OntargetM$avana_CRISPR
sec.prism <- OntargetM$secondary_prism
d.mt <- OntargetM$mutations_mat
sim <- mclapply(S.Drugs,function(x) GetSim(x,sec.prism,KO.GES),mc.cores = 2)
names(sim) <- S.Drugs
Meta.data <- OntargetM$DrugMetadata
DrugTargetSim <- PredTarget(sim,Meta.data)
Drug.Gene.max.sim <- PredMaxSim(sim,Meta.data)
identical ( DrugTargetSim[,1],Drug.Gene.max.sim[,1] )
Pred.d <-cbind (DrugTargetSim,Drug.Gene.max.sim )
DOI = 'dabrafenib'
GOI = 'BRAF'
DMB (DOI,GOI,Pred.d,d.mt,sec.prism,KO.GES)
}

