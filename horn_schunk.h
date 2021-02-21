#ifndef _Horn_Schunk_h_
#define _Horn_Schunk_h_
#define alpha 15
#define Nmax 50 //NOMBRE DES ITERATIONS

void gradient(double** image1,double** image2,double** Ix,double**Iy,double** It,int h,int w);


double R(double** x,int h,int w,int i, int j);


void calcul_u_v(double** u,double** v,double** Ix,double** Iy,double** It,double** uM,double** vM,int h,int w,int n);


#endif

