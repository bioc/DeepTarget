\name{Depmap2DeepTarget}
\alias{Depmap2DeepTarget}
\title{
Retrieve input data required from Depmap to Deeptarget package.
}
\description{
Retrieve gene expression, Cripr, and mutation data and then preperation the matrix compatible as input for Deeptarget.
}
\usage{
Depmap2DeepTarget(DataType,version)
}
\arguments{
    \item{DataType}{Desired Datatype: rna (TPM), cripr, or mutation}
    \item{version}{Version of data}
}
\value{
a matrix of desired datatype.
}
\author{
Ying Hu, Trinh Nguyen
}

\examples{
library("depmap")
library("ExperimentHub")
dat.m <- Depmap2DeepTarget("mutationCalls","20Q4")
head(dat.m)
dat.TPM.rna <- Depmap2DeepTarget("TPM","20Q4")
head(dat.TPM.rna)
dat.crispr <- Depmap2DeepTarget("crispr","19Q1")
head(dat.crispr)
}

