\name{DoInteractMutant}
\alias{DoInteractMutant}
\title{
Compute interaction between the drug and KO expression in term of mutant vs non-mutant
}
\description{
Compute interaction between the drug and KO expression in term of mutant vs non-mutant
}
\usage{
DoInteractMutant(Predtargets,Mutant,DRS,GES)
}

\arguments{
    \item{Predtargets}{a dataframe of drugs information and their most targeted gene with stats of correlation}
    \item{Mutant}{Mutant matrix}
    \item{DRS}{Drug scores matrix}
    \item{GES}{Gene effect scores matrix from KO method}
}
\value{
A list of drug names with their interaction values from two groups mutant and non-mutant
    \item{drug1 }{interaction with estimate and P vals from the linear model}
    \item{drug2 }{interaction with estimate and P vals from the linear model}
    \item{drugN }{interaction with estimate and P vals from the linear model}
}
\author{sanjusinha7, Trinh Nguyen}
\examples{
library(BiocParallel)
data (OntargetM)
set.seed (12345)
All.Drugs <- OntargetM$DrugMetadata[,"broad_id_trimmed"]
S.Drugs <- sample(All.Drugs, 5)
KO.GES <- OntargetM$avana_CRISPR
sec.prism <- OntargetM$secondary_prism
sim <- bplapply(S.Drugs,function(x) computeCor(x,sec.prism,KO.GES))
names(sim) <- S.Drugs
Meta.data <- OntargetM$DrugMetadata
DrugTargetSim <- PredTarget(sim,Meta.data)
d.mt <- OntargetM$mutations_mat
MutantInteract <- DoInteractMutant (DrugTargetSim,d.mt,sec.prism,KO.GES)
}

