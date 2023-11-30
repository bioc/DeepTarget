\name{DoInteractExp}
\alias{DoInteractExp}
\title{
Compute the interaction between the drug and KO expression
}
\description{
Computes interaction between the drug and KO expression in term of lower vs higher expression using linear model.
}
\usage{
DoInteractExp(Predtargets,Exp,DRS, GES,CutOff=3)
}
\arguments{
    \item{Predtargets}{a dataframe of drugs information and their most targeted gene with stats of correlation}
    \item{Exp}{Expression matrix}
    \item{DRS}{Drug scores matrix}
    \item{GES}{Gene effect scores matrix from KO method}
    \item{CutOff}{desired cut-off for low expression}
}
\value{
A list of drug names with their interaction values from two groups low and high expression based on the desired cut-off.
    \item{drug1 }{interaction with estimate and P vals from the linear model}
    \item{drug2 }{interaction with estimate and P vals from the linear model}
    \item{drugN }{interaction with estimate and P vals from the linear model}
}

\author{
sanjusinha7, Trinh Nguyen
}

\examples{
library(BiocParallel)
data (OntargetM)
set.seed (12345)
All.Drugs <- OntargetM$DrugMetadata[,"broad_id_trimmed"]
S.Drugs <- sample(All.Drugs, 5)
KO.GES <- OntargetM$avana_CRISPR
sec.prism <- OntargetM$secondary_prism
sim.out <- bplapply(S.Drugs,function(x) computeCor(x,sec.prism,KO.GES))
names(sim.out) <- S.Drugs
Meta.data <- OntargetM$DrugMetadata
DrugTargetSim <- PredTarget(sim.out,D.M = Meta.data)
d.expr <- OntargetM$expression_20Q4
ExpInteract <- DoInteractExp (DrugTargetSim,d.expr,sec.prism,KO.GES,CutOff = 2)
}
