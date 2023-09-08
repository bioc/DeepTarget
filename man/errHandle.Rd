\name{errHandle}
\alias{errHandle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Custom Error handling
}
\description{
The err_handle() function is a simple function that tries to evaluate an expression x, and if an error occurs it returns NA instead. This is known as error handling, and it allows the code to continue running even if an error is encountered.
}
\usage{
errHandle(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}


\author{
sanjusinha7
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
    tryCatch(x, error = function(e) {
        NA
    })
  }
}

