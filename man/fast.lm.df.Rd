\name{fast.lm.df}
\alias{fast.lm.df}
\title{
  no frills linear regression
}
\description{
  a no frills version of lm.  designed to return only the coefficients
  and nothing else. (and potentially some basic stats in a future version)
}
\usage{
fast.lm.df(panel, right.hand.side, left.hand.sides)
}
\arguments{
  \item{panel}{
    a dataframe containing the relevant vectors
  }
  \item{right.hand.side}{
    a string referencing the column containing the observation vector of
    length \code{n}.
  }
  \item{left.hand.sides}{
    a vector of strings referencing the vectors to use as the design
    matrix (\code{n * p}).
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
NC <- 10
NR <- 1e5
A <- as.data.frame(asofdate=1:NR,matrix(rnorm(NR*NC),nrow=NR,ncol=NC))
colnames(A) <- letters[1:NC]
flm.time <- system.time(ans <- fast.lm.df(panel=A, right.hand.side="a",left.hand.sides=c("b","c","d","e")))
lm.time <- system.time(lm.ans <- as.vector(lm("a ~ 0 + b + c + d + e",data=A)$coefficients))
cat("all.equal:", all.equal(ans,lm.ans),"\n")

bmark.times <- rbind(flm.time,
                     lm.time,
                     lm.time/flm.time)
print(bmark.times)
ans.compare <- cbind(ans,lm.ans)
}
\keyword{regression}
\keyword{array}
