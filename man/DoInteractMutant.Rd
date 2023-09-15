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
    \item{Predtargets}{a table of drug information and their most target gene with stats of correlation}
    \item{Mutant}{Mutant matrix}
    \item{DRS}{Drug scores}
    \item{GES}{Gene effect scores from KO method}
}
\value{
A list of drug names with their interaction values from two groups mutant and non-mutant
    \item{drug1 }{interaction with estimate and P vals}
    \item{drug2 }{interaction with estimate and P vals}
    \item{drugN }{interaction with estimate and P vals}
}
\author{sanjusinha7, Trinh Nguyen}
\examples{
library(parallel)
data (OntargetM)
set.seed (12345)
All.Drugs <- OntargetM$DrugMetadata[,"broad_id_trimmed"]
S.Drugs <- sample(All.Drugs, 5)
KO.GES <- OntargetM$avana_CRISPR
sec.prism <- OntargetM$secondary_prism
sim <- mclapply(S.Drugs,function(x) GetSim(x,sec.prism,KO.GES),mc.cores = 2)
names(sim) <- S.Drugs
Meta.data <- OntargetM$DrugMetadata
DrugTargetSim <- PredTarget(sim,Meta.data)
d.mt <- OntargetM$mutations_mat
MutantInteract <- DoInteractMutant (DrugTargetSim,d.mt,sec.prism,KO.GES)
}

