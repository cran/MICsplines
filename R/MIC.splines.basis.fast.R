MIC.splines.basis.fast <-
function(x, df = NULL, knots = NULL, boundary.knots=NULL, type="Ms", degree = 3, delta=0.01, eq.alloc=FALSE)
{
  if(is.null(df)&is.null(knots))
  {
    stop("either df or knots needs to be supplied")
  }

  if(!is.null(df))
  { 
    if(df<=degree)
    {
      stop("df should be larger than degree")
    }
      
    cc=df-degree+1
    knots=quantile(x,1:(cc-1)/cc)
  }

  if(is.null(boundary.knots))
  {
   boundary.knots=range(x)
   boundary.knots[2]=boundary.knots[2]+0.0001
  }

  if(eq.alloc==T)
  {
    if(is.null(df))
    {
      stop("df needs to be supplied for equal allocation")
    }
    
    if(df<=degree)
    {
      stop("df should be larger than degree")
    }
    
    cc=df-degree+1
    knots=boundary.knots[1]+(1:(cc-1))*(boundary.knots[2]-boundary.knots[1])/cc
  }

  tt=c(rep(boundary.knots[1],degree),knots,rep(boundary.knots[2],degree))
  nn=length(x)
  mm=degree+length(knots)
  mat=matrix(0, nrow=nn,ncol=mm)

  c.delta=delta
  c.ii=integer(nn)
  c.min.tt=min(tt)
  c.max.tt=max(tt)
  xtt=seq(c.min.tt, c.max.tt, by=c.delta)
  nxtt=length(xtt)
  ytt=double(nxtt)
  c.res=double(nn)

  #M splines basis
  if(type=="Ms")
  {
    c.type=1
    c.Cs=0
    ctmp=.C("MIC_splines_basis_C", as.double(x), as.integer(c.ii), as.integer(nn), as.double(tt), as.double(c.min.tt), as.double(xtt), as.double(ytt), as.integer(nxtt), as.double(c.delta), as.integer(c.Cs), as.double(c.res), as.integer(c.type), mat=as.double(mat), as.integer(degree), as.integer(mm), PACKAGE="MICsplines")
   mat=matrix(ctmp[["mat"]], nrow=nn, byrow=F)
   aux.inf=NULL
  }

  #I spline
  if(type=="Is")
  {
    c.type=2
    c.Cs=0
    ctmp=.C("MIC_splines_basis_C", as.double(x), as.integer(c.ii), as.integer(nn), as.double(tt), as.double(c.min.tt), as.double(xtt), as.double(ytt), as.integer(nxtt), as.double(c.delta), as.integer(c.Cs), as.double(c.res), as.integer(c.type), mat=as.double(mat), as.integer(degree), as.integer(mm), PACKAGE="MICsplines")
   mat=matrix(ctmp[["mat"]], nrow=nn, byrow=F)
   aux.inf=colMeans(mat)
   mat=sweep(mat,2,colMeans(mat),"-")

  }
  #I spline not-normalized
  if(type=="IsN")  
  {
    c.type=2
    c.Cs=0
    ctmp=.C("MIC_splines_basis_C", as.double(x), as.integer(c.ii), as.integer(nn), as.double(tt), as.double(c.min.tt), as.double(xtt), as.double(ytt), as.integer(nxtt), as.double(c.delta), as.integer(c.Cs), as.double(c.res), as.integer(c.type), mat=as.double(mat), as.integer(degree), as.integer(mm), PACKAGE="MICsplines")
    mat=matrix(ctmp[["mat"]], nrow=nn, byrow=F)
    aux.inf=NULL
  }

  #C spines
  if(type=="Cs")
  {
    c.type=2
    c.Cs=1
    ctmp=.C("MIC_splines_basis_C", as.double(x), as.integer(c.ii), as.integer(nn), as.double(tt), as.double(c.min.tt), as.double(xtt), as.double(ytt), as.integer(nxtt), as.double(c.delta), as.integer(c.Cs), as.double(c.res), as.integer(c.type), mat=as.double(mat), as.integer(degree), as.integer(mm), PACKAGE="MICsplines")
   mat=matrix(ctmp[["mat"]], nrow=nn, byrow=F)

   x1=x-mean(x)
   x1=x1/sqrt(sum(x1*x1))
   mat=sweep(mat,2,colMeans(mat),"-")
   for(j in 1:mm)
   {
     mat[,j]=mat[,j]-x1*sum(mat[,j]*x1)
   }
   mat=cbind(x1,mat)
  
   aux.inf=NULL
  }
  ###return results
  res=list(mat=mat, x=x, df=df, knots=knots, boundary.knots=boundary.knots,
  type=type, degree=degree, delta=c.delta, aux.inf=aux.inf)
  attr(res,"class")="MICsplines"
  return(res)
}
