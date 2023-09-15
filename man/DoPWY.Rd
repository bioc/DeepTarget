\name{DoPWY}
\alias{DoPWY}
\title{Provide a probability score for each pathway to be the MOA of a drug}
\description{
Predicts a Primary Target at a pathwya Level. It next finds the pathways that are most enriched in the genes with high DKS scores. It does this by performing a pathway enrichment test on the ranked gene list by DKS score. The output is a matrix of pathway-level probabilities for each drug to be the primary MOA
}

\usage{
DoPWY(Sim.GES.DRS,D.M)
}

\arguments{
    \item{Sim.GES.DRS}{The list of result from "CorCRISPRDrugR" function.}
    \item{D.M}{meta data from drug}
}

\value{
a dataframe contain the pathway level probability to be a primary MOA
    \item{drug1 }{a dataframe contain the pathway level probability to be a primary MOA}
    \item{drug2 }{a dataframe contain the pathway level probability to be a primary MOA}
    \item{drugN }{a dataframe contain the pathway level probability to be a primary MOA}
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
Pwy.Enr <- DoPWY(sim,Meta.data)
}

