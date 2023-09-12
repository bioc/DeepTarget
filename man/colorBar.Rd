\name{colorBar}
\alias{colorBar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
create a key bar for a plot
}
\description{
create a key bar for a plot
}
\usage{
colorBar(xl, yb, xr, yt, v.min, v.max, col.l, alpha=1, color.spec="hsv", title ="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xl}{
x axis, left
}
  \item{yb}{
y axis, bottom
}
  \item{xr}{
x, right
}
  \item{yt}{
y, top
}
  \item{v.min}{
mininum value
}
  \item{v.max}{
max value
}
  \item{col.l}{
%%     ~~Describe \code{col.l} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{color.spec}{
%%     ~~Describe \code{color.spec} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
}
\author{
Ying Hu, Trinh Nguyen
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
plot(x, y, pch = 19, col = y.clr , xlab="-log10(Pval)", ylab="correlation strength");
x1   <- (min(x) + max(x))/2 - 2;
x2   <- (min(x) + max(x))/2 + 2;
y1   <- max(y) + (min(abs(y)) + max(abs(y))) * 0.40;
y2   <- max(y) + (min(abs(y)) + max(abs(y))) * 0.55;
colorBar(x1, y1, x2, y2, col.l=color(100), v.min=min(y), v.max=max(y))
}
