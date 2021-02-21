# ifndef Gauss_h_
# define Gauss_h_
#include"utilitaire.h"
#include <math.h>
#define pi 3.14159265358979323846
#define N_pyr 100
#define alpha_u 15
#define  gaussDim 3
#define sigma 0.8

double GaussFiltre(int x,int y);

double** GaussMatrix();

double** Convlution(double** h, int dimX,int dimY,double** x);

double** ProduitMatTermeParTerme(double** terme1,double** terme2,int dimX,int dimY );

void CalcU_V_avecConv(double**U,double**V,double** ubar ,double ** vbar, double** FG , double** Ix , double** Iy,double** It,int lignes,int cols);
double **reduire_res(double** matIn,int dimX,int dimY);

double** augment_res(double** matIn,int dimX,int dimY);
# endif

