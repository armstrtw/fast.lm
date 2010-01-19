fast.lm <- function(A,x) {
    .Call("fast_lm",A,x,PACKAGE="fast.lm")
}

fast.lm.df <- function(panel, right.hand.side, left.hand.sides) {
    .Call("fast_lm_dataframe",panel, right.hand.side, left.hand.sides, PACKAGE="fast.lm")
}

expanding.lm.df <- function(panel, right.hand.side, left.hand.sides, min.rows) {
    .Call("expanding_lm_dataframe",panel, right.hand.side, left.hand.sides, min.rows, PACKAGE="fast.lm")
}

expanding.panel.df <- function(panel, right.hand.side, left.hand.sides, asofdate.column="asofdate", min.dates=5L) {
    .Call("expanding_panel_dataframe", panel, right.hand.side, left.hand.sides, asofdate.column, min.dates, PACKAGE="fast.lm")
}
