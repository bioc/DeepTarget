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
  \item{Sim.GES.DRS}{
similarity between Drug’s response scores and Gene effect scores from Knock-out method such as CRISPR.
}
  \item{D.M}{
      Drug Metadata
}
}
\value{
a table of drug information and the most correlated gene across all genes.
}

\author{
 sanjusinha7, Trinh Nguyen
}
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
Drug.Gene.max.sim <- PredMaxSim(sim,Meta.data)
}

