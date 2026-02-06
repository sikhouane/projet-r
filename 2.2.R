f_val <- function(x, y) exp(x*y) - exp(x) + 2

f_grad <- function(u) {
  x <- u[1]; y <- u[2]
  c(y*exp(x*y) - exp(x),
    x*exp(x*y))
}

