\name{expanding.lm.df}
\alias{expanding.lm.df}
\title{
  no frills rolling linear regression
}
\description{
  a no frills rolling version of lm.  designed to return only the coefficients
  and nothing else. (and potentially some basic stats in a future version)

}
\usage{
expanding.lm.df(panel, right.hand.side, left.hand.sides, min.rows)
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
  \item{min.rows}{
    required number of rows of observations before rolling lm can start.
  }
}
\value{
  a matrix of the coefficients over the expanding window.
}
\author{
  Whit Armstrong
}
\seealso{
  \code{\link{lm}} the standard lm method.
}
\examples{

do.running.lm <- function(panel, right.hand.side, left.hand.sides, min.obs=30L) {
    NR <- nrow(panel)
    ans <- matrix(NA, NR - min.obs + 1, length(left.hand.sides))
    formula <- paste(right.hand.side,paste(c("0",left.hand.sides),collapse="+"),sep="~")
    for(i in min.obs:NR) {
        ans[i-min.obs+1,] <- as.vector(lm(formula=formula,data=panel[1:i,])$coefficients)
    }
    ans
}

NC <- 5
NR <- 1e3
A <- as.data.frame(asofdate=1:NR,matrix(rnorm(NR*NC),nrow=NR,ncol=NC))
colnames(A) <- letters[1:NC]
flm.time <- system.time(ans <- expanding.lm.df(panel=A, right.hand.side="a",left.hand.sides=c("b","c","d","e")))
lm.time <- system.time(lm.ans <- do.running.lm(panel=A, right.hand.side="a",left.hand.sides=c("b","c","d","e")))
cat("all.equal:", all.equal(ans,lm.ans),"\n")
bmark.times <- rbind(flm.time,
                     lm.time,
                     lm.time/flm.time)
print(bmark.times)

}
\keyword{regression}
\keyword{array}
