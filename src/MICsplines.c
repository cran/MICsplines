#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>


//******************************************************************************
double m_spline_x(double x, double *tt, int i, int k)
{
  double res, ti, tik, a1, a2;
  
  ti=tt[i];
  tik=tt[i+k];
  
  if((x<ti)|(x>=tik))
  {
     res=0.0;
     return res;
  }else{
        if(k==1)
        {
          res=1/(tik-ti);
          return res;
        }else{
              a1=m_spline_x(x,tt,i,k-1);
              a2=m_spline_x(x,tt,i+1,k-1);
              res=(k*((x-ti)*a1+(tik-x)*a2))/((k-1)*(tik-ti));
              return res;
             }
       }
}
//******************************************************************************
void i_spline_x(double *xx, int *ii, int *nxx, double *tt, double *min_tt, double *xtt, double *ytt, int *nxtt, int i, int k, double *delta, int *Cs, double *res)
{
  int j;

  for(j=0; j<(*nxtt); j++)
  {
    ytt[j]=m_spline_x(xtt[j],tt,i,k);
  }

  ytt[0]=ytt[0]*(*delta);
  
  for(j=1; j<(*nxtt); j++)
  {
    ytt[j]=ytt[j-1]+ytt[j]*(*delta); //?rieman integral??
  }
  
  if((*Cs)==1)
  {
    ytt[0]=ytt[0]*(*delta);

    for(j=1; j<(*nxtt); j++)
    {
      ytt[j]=ytt[j-1]+ytt[j]*(*delta);
    }
   }

  for(j=0; j<(*nxx); j++)
  {
    ii[j]=floor((xx[j]-(*min_tt))/(*delta));
    if(ii[j]>=(*nxtt))
    {
      ii[j]=(*nxtt)-1;
    }
  }

  for(j=0; j<(*nxx); j++)
  {
    res[j]=ytt[ii[j]];
  }
  
  return;
}
//******************************************************************************
void MIC_splines_basis_C(double *xx, int *ii, int *nxx, double *tt, double *min_tt, double *xtt, double *ytt, int *nxtt, double *delta, int *Cs, double *res, int *type, double *mat, int *degree, int *mm)
{
  int i, j;
  //xx: 
  
  if((*type)==1)
  {
   for(i=0; i<(*nxx); i++)
   {
    for(j=0; j<(*mm); j++)
    {
      mat[i+(*nxx)*j]=m_spline_x(xx[i],tt,j,(*degree)); ///the value of spline basis is calculated at each xx position
    }
   }
  }

  if((*type)==2)
  {
   for(j=0; j<(*mm); j++)
   {
     i_spline_x(xx, ii, nxx, tt, min_tt, xtt, ytt, nxtt, j, (*degree), delta, Cs, res); //xtt grid is used for integral????

     for(i=0; i<(*nxx); i++)
     {
       mat[i+(*nxx)*j]=res[i];
     }
   }
  }

}
