#' Dichotomization
#' 
#' This function takes a real value and converts it to 1 if it is positive and
#' else to 0.
#' 
#' 
#' @param val A real value
#' @return 0 or 1.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@math.unistra.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link{ifelse}}
#' @keywords utilities
#' @examples
#' 
#' dimX <- 6
#' Astar <- 4
#' (dataAstar4 <- t(replicate(10,simul_data_YX(dimX,Astar))))
#' 
#' dicho(dataAstar4)
#' 
#' rm(list=c("dimX","Astar"))
#' 
#' @export dicho
dicho <- function(val) {ifelse(val>0,1,0)}
