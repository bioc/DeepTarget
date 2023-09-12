\name{GetSim}
\alias{GetSim}
\title{Compute a correlation between the every gene crispr KO vs each drug response}
\description{
Computes correlations between the viability of cell lines after CRISPR knockout of each gene and the viability of the same cell lines after drug treatment.
}
\usage{GetSim(DrugName,DRS,GES)}
\arguments{
  \item{DrugName}{
     Drug Name
}
  \item{DRS}{
      Drug’s response scores
}
  \item{GES}{
      Gene effect scores from Knock-out method such as CRISPR.
}
}
\value{
A dataframe containing genes with their correlation values and P values associated with response scores from a given drug ID.
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
head(sim)
}

