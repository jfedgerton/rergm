is.inCH <- function(p, M) {
  p <- as.vector(p)
  if (!is.matrix(M)) 
    stop("Second argument must be a matrix.")
  if (length(p) != NCOL(M)) 
    stop("Number of columns in matrix (2nd argument) is not equal to dimension ",
         "of first argument.")
  R = NROW(M)
  C=length(p)+1
  if(require("Rglpk", quietly = TRUE)) {
    ans <- Rglpk_solve_LP(obj=c(p,-1), mat=cbind(rbind(p,M),-1),
                          dir=as.vector(rep("<=",R+1)), rhs=as.vector(c(1,rep(0,R))),
                          max=TRUE, bounds=list(lower=list(ind=1:C,val=rep(-Inf,C))))
  }else{
    stop("The 'Rglpk' package must be installed if the ergm control parameter 'style'='stepping' ")
  }
  
  if(ans$optimum==0){x<-TRUE}  #if the max is zero, the point q is in the CH of the points p
  else{x<-FALSE}              #if the max is strictly positive, q is not in the CH of p
  x
}
