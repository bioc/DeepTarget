\name{PredTarget}
\alias{PredTarget}
\title{
Prediction  of the most similar known targeted gene.
}
\description{
Predicts the gene that has the most similarity to a drug’s response scores. This is done based on selecting a gene that has the most correlation across the known targeted genes by their drug.
}
\usage{
PredTarget(Sim.GES.DRS,D.M)
}
\arguments{
    \item{Sim.GES.DRS}{similarity between Drug’s response scores and Gene effect scores from Knock-out method such as CRISPR.}
    \item{D.M}{Drug Metadata}
}
\value{
a dataframe of drug(s) information with the most known predicted gene(s) with the max corelation value(s), P value(s), and FDR value(s).
}
\author{
sanjusinha7, Trinh Nguyen
}
\examples{
library(BiocParallel)
data(OntargetM)
set.seed (12345)
All.Drugs <- OntargetM$DrugMetadata[,"broad_id_trimmed"]
S.Drugs <- sample(All.Drugs, 5)
KO.GES <- OntargetM$avana_CRISPR
sec.prism <- OntargetM$secondary_prism
sim.out <- bplapply(S.Drugs,function(x) computeCor(x,sec.prism,KO.GES))
names(sim.out) <- S.Drugs
Meta.data <- OntargetM$DrugMetadata
DrugTargetSim <- PredTarget(sim.out,Meta.data)
}
