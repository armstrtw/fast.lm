fast.lm <- function(A,x) {
    .Call("fast_lm",A,x,PACKAGE="fast.lm")
}

fast.lm.df <- function(panel, right.hand.side, left.hand.sides) {
    .Call("fast_lm_dataframe",panel, right.hand.side, left.hand.sides, PACKAGE="fast.lm")
}

expanding.lm.df <- function(panel, right.hand.side, left.hand.sides, min.rows=30L) {
    .Call("expanding_lm_dataframe",panel, right.hand.side, left.hand.sides, as.integer(min.rows), PACKAGE="fast.lm")
}

expanding.panel.df <- function(panel, right.hand.side, left.hand.sides, asofdate.column="asofdate", min.dates=5L) {
    .Call("expanding_panel_dataframe", panel, right.hand.side, left.hand.sides, asofdate.column, as.integer(min.dates), PACKAGE="fast.lm")
}

group.lm.df <- function(panel, right.hand.side, left.hand.sides, groups, blend) {
    .Call("group_lm_dataframe", panel, right.hand.side, left.hand.sides, groups, blend, PACKAGE="fast.lm")
}


summary.lm <- function(object)
{
    z <- object
    p <- z$rank
    Qr <- object$qr
    n <- NROW(Qr$qr)
    rdf <- n - p
    p1 <- 1L:p
    ## do not want missing values substituted here
    r <- z$residuals
    f <- z$fitted.values
    mss <- if (attr(z$terms, "intercept")) sum((f - mean(f))^2) else sum(f^2)
    rss <- sum(r^2)
    resvar <- rss/rdf
    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    se <- sqrt(diag(R) * resvar)
    est <- z$coefficients[Qr$pivot[p1]]
    tval <- est/se
    ans <- z[c("call", "terms")]
    ans$residuals <- r
    ans$coefficients <- cbind(est, se, tval, 2*pt(abs(tval), rdf, lower.tail = FALSE))
    ans
}
