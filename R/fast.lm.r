fast.lm <- function(A,x) {
    .Call("fast_lm",A,x,PACKAGE="fastlm")
}

panel.lm <- function(panel, right.hand.side, left.hand.sides, asofdate.column="asofdate") {
    .Call("panel_lm",panel, right.hand.side, left.hand.sides, asofdate.column)
}

expanding.panel <- function(panel, right.hand.side, left.hand.sides, asofdate.column="asofdate", min.dates=30L) {
    .Call("expanding_panel",panel, right.hand.side, left.hand.sides, asofdate.column, as.integer(min.dates))
}
