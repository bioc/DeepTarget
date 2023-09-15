\name{cor.test_trimmed_v0.default}
\alias{cor.test_trimmed_v0.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Trimmed cor test that is faster than standard cor.test
}
\description{
Calculation the corelation. This will be called inside other function such as Getsim{}
}
\usage{
cor.test_trimmed_v0.default(x, y, alternative = c("two.sided", "less", "greater"), method = c("pearson"), exact = NULL, continuity = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{x}{Drug response matrix}
    \item{y}{avana_CRISPR Matrix}
    \item{alternative}{option can be chosen as "two.sided", "less", "greater"}
    \item{method}{default is "pearson"}
    \item{exact}{default is NULL}
    \item{continuity}{default is FALSE}
    \item{\dots}{ default is NULL}
}
\details{
this function is called by another main function from DeepTarget
}
\references{
The code is adopted from https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/cor.test.R
}
\value{
corellation values of the drug scores and Knock-out effect scores.
}
\author{
SurajGupta
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
colnames(out.Sim.GES.DRS) <- row.names(KO.GES)
}
