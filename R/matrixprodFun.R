##
## C implementation of traces of products
##

trX <- function(x, PACKAGE="gRcox") {
  .Call("tr", x, PACKAGE="gRcox")
}

trXY <- function(x, y, PACKAGE="gRcox") {
  .Call("trProd", x, y, PACKAGE="gRcox")
}

trXYXY <- function(x, y, PACKAGE="gRcox") {
  .Call("trProd2", x, y, PACKAGE="gRcox")
}

trXYZY <- function(x, y, z, PACKAGE="gRcox") {
  .Call("trProd3", x, y, z, PACKAGE="gRcox")
}

trXYZ <- function(x, y, z, PACKAGE="gRcox") {
  .Call("trXYZ", x, y, z, PACKAGE="gRcox")
}


