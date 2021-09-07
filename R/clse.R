clse <-
function(dat.obj)   #constrained least-squares estimation
{
  lam=dat.obj$lam
  ncoef=length(lam)
  beta.vec=double(ncoef)
  xmats=dat.obj$mat
  y=dat.obj$y

  aa1=as.matrix(xmats[,lam==0])
  aa2=as.matrix(xmats[,lam==1])
  aa=cls(y=y,X=aa2-Px(aa1)%*%aa2)

  beta2=aa$betahat
  beta1=solve(t(aa1)%*%aa1)%*%t(aa1)%*%(y-aa2%*%beta2)

  beta.vec[lam==0]=beta1
  beta.vec[lam==1]=beta2

  yhat=xmats%*%beta.vec
  yhat=as.vector(yhat)

  res=list(dat.obj=dat.obj, beta.vec=beta.vec, yhat=yhat)
  
  return(res)
}
