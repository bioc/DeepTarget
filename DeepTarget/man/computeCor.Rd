\name{computeCor}
\alias{computeCor}
\title{Compute a correlation between the every gene vs each drug response}
\description{
Compute correlations between the viability of cell lines after CRISPR Knock Out of each gene and of the same cell lines after drug treatment.
}
\usage{computeCor(DrugName,DRS,GES)}
\arguments{
    \item{DrugName}{Drug Name}
    \item{DRS}{Drug’s response scores}
    \item{GES}{Gene effect scores from Knock-out method such as CRISPR.}
}
\value{
a list of matrices for the interesting drugs, where each matrix containing gene names with the correlation values and P values associated with response scores from a given drug ID.
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
names(sim.out ) <- S.Drugs
head(sim.out)
}

