\name{expanding.panel.df}
\alias{expanding.panel.df}
\title{
  no frills expanding window linear regression for panel data.
}
\description{
  a no frills version of lm.  designed to return only the coefficients
  and nothing else. (and potentially some basic stats in a future
  version)

  a regression is performed for each unique time step in the data set.
  and the coefficients are returned as a time series.
}
\usage{
expanding.panel.df(panel, right.hand.side, left.hand.sides, asofdate.column = "asofdate", min.dates = 5L)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{asofdate.column}{
    string referencing the column of \code{panel} to use as the dates vector.
  }
  \item{min.dates}{
    minimum number of dates required to start the expanding window.
  }
}
\value{
  a time series of the coefficients over the expanding window.
}
\author{
  Whit Armstrong
}
\seealso{
  \code{\link{lm}} the standard lm method.
}
\examples{
library(fast.lm)

do.running.lm <- function(panel, right.hand.side, left.hand.sides, asofdate.colname="asofdate",min.dates=5L) {
    dts <- panel[,asofdate.colname]
    unique.dates <- unique(dts)
    NR <- nrow(panel)
    ans <- matrix(NA, length(unique.dates)  - min.dates + 1, length(left.hand.sides))
    formula <- paste(right.hand.side,paste(c("0",left.hand.sides),collapse="+"),sep="~")
    for(i in min.dates:length(unique.dates)) {
        ans[i-min.dates+1,] <- as.vector(lm(formula=formula,data=panel[dts <= unique.dates[i],])$coefficients)
        cat(".")
    }
    cat("\n")
    ans
}

NC <- 5
NR <- 1e3
panels <- 10
asofdate <- as.double(rep(1:(NR/panels),each=panels))
A <- data.frame(asofdate=asofdate,matrix(rnorm(NR*NC),nrow=NR,ncol=NC))
colnames(A) <- c("asofdate",letters[1:NC])
flm.time <- system.time(ans <- expanding.panel.df(panel=A, right.hand.side="a",left.hand.sides=c("b","c","d","e")))
lm.time <- system.time(lm.ans <- do.running.lm(panel=A, right.hand.side="a",left.hand.sides=c("b","c","d","e")))
cat("all.equal",all.equal(ans,lm.ans),"\n")
bmark.times <- rbind(flm.time,
                     lm.time,
                     lm.time/flm.time)
print(bmark.times)

}
\keyword{regression}
\keyword{array}
