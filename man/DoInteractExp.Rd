\name{DoInteractExp}
\alias{DoInteractExp}
\title{
Compute the interaction between the drug and KO expression
}
\description{
Computes interaction between the drug and KO expression in term of lower vs higher expression
}
\usage{
DoInteractExp(Predtargets,Exp,DRS, GES,CutOff=3)
}
\arguments{
 \item{Predtargets}{a table of drug information and their most targeted gene with stats of correlation}

\item{Exp}{Expression matrix}

\item{DRS}{Drug scores}

\item{GES}{Gene effect scores from KO method}

\item{CutOff}{desired cut-off for low expression}
}

\value{
A list of drug names with their interaction values from two groups low and high expression.
    \item{drug1 }{interaction with estimate and P vals}
    \item{drug2 }{interaction with estimate and P vals}
     \item{drugN }{interaction with estimate and P vals}
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
DrugTargetSim <- PredTarget(sim, D.M = Meta.data)
d.expr <- OntargetM$expression_20Q4
ExpInteract <- DoInteractExp (DrugTargetSim,d.expr,sec.prism,KO.GES,CutOff = 2)
}
