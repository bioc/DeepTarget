\name{plotSim}
\alias{plotSim}
\title{
Plot the similarty between corelation values and P val;
}
\description{
Plot the similarty between corelation values and P val;
}
\usage{
plotSim(dx,dy,clr=NULL, plot=TRUE)
}
\arguments{
    \item{dx}{a matrix of p vals}
    \item{dy}{a matrix of correlation vals}
    \item{clr}{Desired color}
    \item{plot}{ default plot =TRUE}
}
\value{a plot of similarity}
\author{
Ying Hu,Trinh Nguyen
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
cor=sapply(sim, function(x) x[,2])
Pval=sapply(sim, function(x) x[,1])
par(mar=c(4,4,5,2), xpd=TRUE, mfrow=c(3,3));
color <- colorRampPalette(c("lightblue",'darkblue'));
plotSim(Pval,cor,clr=color);
}
