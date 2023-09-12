\name{line2user}
\alias{line2user}
\title{
A utility function to correctly place the colorBar.
}
\description{
This function is used in conjunction with the colorBar function to create color bar placed in top middle of a plot.
}
\usage{
line2user(line, side)
}
\arguments{
  \item{line}{
}
  \item{side}{
%%     ~~Describe \code{side} here~~
}
}
\author{
Ying Hu <yhu@mail.nih.gov>
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
color <- colorRampPalette(c("lightblue",'darkblue'));
x <- -log10 (Pval[,1])
y <- cor[,1]
y.clr <- color(10)[as.numeric(cut(y,breaks = 10))]
par(mar=c(4,4,5,2), xpd=TRUE, mfrow=c(1,1));
plot(x, y, pch = 19, col = y.clr , xlab="-log10(Pval)", ylab= "correlation strength");
x1   <- (min(x) + max(x))/2 - line2user(line=2, side=4) * 0.25;
x2   <- (min(x) + max(x))/2 + line2user(line=2, side=4) * 0.25;
y1   <- line2user(line=2, side=3);
y2   <- line2user(line=3, side=3);
colorBar(x1, y1, x2, y2, col.l=color(100), v.min=min(y), v.max=max(y))
}