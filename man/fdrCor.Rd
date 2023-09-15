\name{fdrCor}
\alias{fdrCor}
\title{Calculate FDR}
\description{this function is used to calculate for FDR based on the list of P val}
\usage{fdrCor(test_list)}
\arguments{
    \item{test_list}{ a list of P val}
}
\value{
a table with FDR values
}
\author{
sanjusinha7, Trinh Nguyen
}
\examples{
library(parallel)
data (OntargetM)
set.seed (12345)
All.Drugs <- OntargetM$DrugMetadata[,"broad_id_trimmed"]
S.Drugs <- sample(All.Drugs, 1)
KO.GES <- OntargetM$avana_CRISPR
sec.prism <- OntargetM$secondary_prism[S.Drugs,]
out.Sim.GES.DRS <- sapply(1:nrow(KO.GES),function(x)
    unlist(cor.test_trimmed_v0.default(sec.prism, KO.GES[x,])))
out.Sim.GES.DRS <- t(out.Sim.GES.DRS)
row.names(out.Sim.GES.DRS) <- row.names(KO.GES)
FDRVal <- fdrCor(out.Sim.GES.DRS[,1])
}