#include"Edo.h"


void jacY(double const x, double * y, double *fl, double *param)
{
/*    | df1/dy1 df1/dy2 df1/dx |
  J = | df2/dy1 df2/dy2 df2/dx |
      |     0.0    0.0     1.0 |
*/

/*... df1/dy1 df1/dy2*/
  fl[0] = 0.e0;
  fl[1] = 1.0; 
/*... df2/dy1  df2/dy2*/
  fl[2] = -(2.e0*y[0] * y[1] + 1.e0) / param[0];
  fl[3] = (1.e0 - y[0]*y[0])/param[0];
}

void jacX(double const x, double * y, double *fl, double *param)
{
/*... df1/x df2/dx*/
  fl[0] = 0.e0;
  fl[1] = 0.e0;
}


void rhs(DOUBLE const x,DOUBLE * y, DOUBLE *f, DOUBLE *param)
{
  f[0] = y[1];
  f[1] = ((1.e0 - y[0] * y[0])*y[1] - y[0]) / param[0];
}


void rhs2(DOUBLE const x,DOUBLE * y, DOUBLE *f, DOUBLE *param)
{
  f[0] = -0.013*y[0] - 1000.0*y[0]*y[2];
  f[1] = -2500.0*y[1]*y[2];
  f[2] = -0.013*y[0]-1000.0*y[0]*y[2]-2500.0*y[1]*y[2];
}

void jacX2(double const x, double * y, double *fl, double *param)
{
/*... df1/x df2/dx*/
  fl[0] = 0.e0;
  fl[1] = 0.e0;
  fl[2] = 0.e0;
}

void jacY2(double const x, double * y, double *fl, double *param)
{

/*
      | df1/dy1 df1/dy2 df1/dy3 |
  J = | df2/dy1 df2/dy2 df2/dy3 |
      | df3/dy1 df3/dy2 df3/dy3 |
*/

/*... df1/dy1 df1/dy2 df1/dy3 */
  fl[0] = -0.013-1000.0*y[2];
  fl[1] = 0.0;
  fl[2] = -1000.0*y[0];
/*... df2/dy1 df2/dy2 df2/dy3*/
  fl[3] = 0.0;
  fl[4] =-2500.0*y[2];
  fl[5] =-2500.0*y[1];
/*...df3/dy1 df3/dy2 df3/dy3*/
  fl[6] = -0.013-1000.0*y[2];
  fl[7] =-2500.0*y[2];
  fl[8] =-1000.0*y[0] - 2500.0*y[1];

}

int main(int argc, char *argv[])
{
  short nEdo = 2;
  unsigned INT it,maxInt;
  DOUBLE y[3], param[5], x1, x2, h, tol;

  param[0] = 1.e-07;
  y[0]     =   2.0;
  y[1]     =   0.0;
//y[2] = 0.0;
  x1      = 0.0;
  x2      = 2.0;
  h       = 1.0e-06;
  maxInt = 10000000;
  tol    = 1.e-06;
/*it = stepperDopr853(y    , param
                , x1       , x2
                , h        , tol 
                , tol      ,  nEdo
                , maxInt   , true  
                , FILE_OUT ,"dopr853.out"
                , &rhs);
  printf("%d\n",it); */
  
  y[0] =   2.0;
  y[1] =   0.0;
  it = stepperRoss(y     ,param
             ,x1         ,x2
             ,h          ,tol
             , tol       , nEdo
             , maxInt   , true   
             , FILE_OUT , "ross.out"
             , &rhs     , &jacY 
             , &jacX);  
  printf("%d\n",it);   

  y[0] =   2.0;
  y[1] =   0.0;
  it = StepperSie(y     , param
             ,x1        , x2
             ,h         , tol
             , tol      , nEdo
             , maxInt   , true   
             , FILE_OUT , "steppersie.out"
//           , &rhs     , &jacY); 
             , &rhs     , NULL);  
  printf("%d\n",it); 


  return EXIT_SUCCESS;

}