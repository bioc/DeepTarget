\name{PredTarget}
\alias{PredTarget}
\title{
Prediction  of the most similar known targeted gene.
}
\description{
Predicts the gene that has the most similarity to a drug’s response scores across targeted genes
}

\usage{

PredTarget(Sim.GES.DRS,D.M)
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
a table of drug information and their most targeted gene with stats of correlation.
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
DrugTargetSim <- PredTarget(sim,Meta.data)
}
