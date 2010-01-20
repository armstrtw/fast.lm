\name{fast.lm}
\alias{fast.lm}
\title{
  no frills linear regression
}
\description{
  a no frills version of lm.  designed to return only the coefficients
  and nothing else. (and potentially some basic stats in a future version)
}
\usage{
fast.lm(A, x)
}
\arguments{
  \item{A}{
    design matrix of dimension \code{n * p}.
  }
  \item{x}{
    vector of observations of length \code{n}, or a matrix with \code{n} rows.
  }
}
\value{
  the length \code{p} coefficient vector.
}
\author{
  Whit Armstrong
}
\seealso{
  \code{\link{lm}} the standard lm method.
}
\examples{
NC <- 5
NR <- 1e3
A <- matrix(rnorm(NR*NC),nrow=NR,ncol=NC)
b <- rnorm(NR)
flm.time <- system.time(ans <- fast.lm(A,b))
lm.fit.time <- system.time(lm.ans <- as.vector(lm.fit(A,b)$coefficients))
cat("all.equal:", all.equal(ans,lm.ans),"\n")

bmark.times <- rbind(flm.time,
                     lm.fit.time,
                     lm.fit.time/flm.time)
print(bmark.times)
ans.compare <- cbind(ans,lm.ans)
}
\keyword{regression}
\keyword{array}
