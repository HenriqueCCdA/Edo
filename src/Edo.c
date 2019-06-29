#include<Edo.h>

/*********************************************************************
 * Data de criacao    : 00/00/0000                                   *
 * Data de modificaco : 19/05/2018                                   *
 * ------------------------------------------------------------------*
 * openFile : abre um arquivo                                        *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * name  -> nome do arquivo                                          * 
 * mod   -> mode de abertura                                         *
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * FILE * ponterio par o arquivo de saida                            * 
 * ------------------------------------------------------------------*
 * *******************************************************************/
FILE* openFile(const char* const name,const char* const mod)
{
  FILE *aux;
  if(( aux =fopen(name,mod))==NULL)
  {
      fprintf(stderr,"Erro na abertura do arquivo %s.\n",name);
      exit(EXIT_FAILURE);
  }

  return aux;  
}
/********************************************************************/

/*********************************************************************
 * Data de criacao    : 19/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * writeOutput :                                                     *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------* 
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * *******************************************************************/
void writeOutput(int const it     , double const x
                , double const h  , double *y
                , short const nEdo, short const outCod
                , FILE *file )
{
  short i;

  switch (outCod)
  {
  case SCREEN_OUT:
    printf("%d %12lf %12.6lf", it, x, h);
    for (i = 0; i<nEdo; i++)
      printf(" %.14e ", y[i]);
    printf("\n");
  case FILE_OUT:
    fprintf(file, "%12d %.14e %.14e", it, x, h);
    for (i = 0; i<nEdo; i++)
      fprintf(file, " %.14e ", y[i]);
    fprintf(file, "\n");
    break;
  }
}
/*********************************************************************/

/*********************************************************************
 * Data de criacao    : 19/06/2019                                   *
 * Data de modificaco : 00/00/0000                                   *
 * ------------------------------------------------------------------*
 * edoError:                                                         *
 * ------------------------------------------------------------------*
 * Parametros de entrada:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * Paramanetros de saida:                                            *
 * ------------------------------------------------------------------*
 * ------------------------------------------------------------------*
 * OBS:                                                              *
 * *******************************************************************/
double edoError(DOUBLE *y        ,DOUBLE *yOut
               ,DOUBLE *yErr     ,DOUBLE const atol
               ,DOUBLE const rtol,short const nEdo)
{
  short i;
  double err = 0.e0, sk;
  for (i = 0; i<nEdo; i++) {
    sk = atol + rtol * max(fabs(y[i]), fabs(yOut[i]));
    err += (yErr[i] / sk)*(yErr[i] / sk);
  }
  return sqrt(err / nEdo);
}
/********************************************************************/

void testAlloc(void* pt,char *const name)
{
  if(pt==NULL)
  {
    printf("Erro na alocacao do vetor %s!!!",name);
    exit(EXIT_FAILURE);
  }
}

DOUBLE sqr(DOUBLE const x)
{
  return x*x;
}