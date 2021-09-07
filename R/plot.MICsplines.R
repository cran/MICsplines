plot.MICsplines <-
function(x, ...)
{
  obj=x
  mat=obj$mat
  x=obj$x
  matplot(x[order(x)],mat[order(x),],type="l",xlab="",ylab="",main=paste(obj$type,"Splines Basis"))
  
}
