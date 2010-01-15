fastlm <- function(A,x) {
    .Call("fast_lm",A,x,PACKAGE="fastlm")
}
