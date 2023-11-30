\name{PredMaxSim}
\alias{PredMaxSim}
\title{
Predict the most similar gene to the drug response
}
\description{
Predicts the gene that has the most similarity associated with drug’s response scores from the set of all genes.
}
\usage{
PredMaxSim (Sim.GES.DRS,D.M)
}
\arguments{
    \item{Sim.GES.DRS}{similarity between Drug’s response scores and Gene effect scores from Knock-out method such as CRISPR}
    \item{D.M}{Drug Metadata}
}
\value{
a dataframe of drug(s) information with the most predicted gene(s) with the max corelation value(s), P value(s), and FDR value(s).
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
sim.out <- bplapply(S.Drugs,function(x) computeCor(x,sec.prism,KO.GES))
names(sim.out) <- S.Drugs
Meta.data <- OntargetM$DrugMetadata
Drug.Gene.max.sim <- PredMaxSim(sim.out,Meta.data)
}

