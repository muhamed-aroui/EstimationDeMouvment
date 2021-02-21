#include"Horn_Schunk.h"
#include "utilitaire.h"
void gradient(double** image1,double** image2,double** Ix,double**Iy,double** It,int h,int w)
{
    int i,j;
    for(i = 0; i < h-1; i++)
        {
            for(j = 0; j < w-1; j++)
            {
                if(j==0) Ix[i][j]=image1[i][1]-image1[i][0];
                else if(j==w-1) Ix[i][j]=image1[i][w-1]-image1[i][w-2];
                else Ix[i][j]=0.5*(image1[i][j+1]-image1[i][j-1]);

                if(i==0) Iy[i][j]=image1[1][j]-image1[0][j];
                else if(j==w-1) Iy[i][j]=image1[h-1][j]-image1[h-2][j];
                else Iy[i][j]=0.5*(image1[i+1][j]-image1[i-1][j]);

                It[i][j]=image2[i][j]-image1[i][j];
            }
        }
}


/** Methode de Recadrement de l'image ( la matrice )
*/
double R(double** x,int h,int w,int i, int j){  // Reencadrement de l'image  w= largeur h=hauteur

	if (i < 0) i = 0;
	if (j < 0) j = 0;
	if (i >= w) i = w - 1;
	if (j >= h) j = h - 1;
	return x[i][j];
}


/** Methode de calcul LES MATRICES : U , V
*/
void calcul_u_v(double** u,double** v,double** Ix,double** Iy,double** It,double** uM,double** vM,int h,int w,int n){  // procedure qui calcule u et v par la methode de horn & schrunk

   if(w%2!=0)w--;
   if(h%2!=0)h--;

    double fu,fv;
    int i=0,j=0;


n=0;
    while (n<=Nmax ) {

          //  compute_one_iteration(Ix,Iy,It,u,v,uM,vM,w,h);
        for (i = 0; i < h-1; i++){
            for (j = 0; j < w-1; j++){
            uM[i][j] = (1.0/8) * ( R(u,h,w,i-1,j-1)+ R(u,h,w,i-1,j)+ R(u,h,w,i-1,j-1)+
                                     R(u,h,w,i-1,j+1)+ R(u,h,w,i,j-1)+ R(u,h,w,i,j+1)+
                                    R(u,h,w,i+1,j-1)+ R(u,h,w,i+1,j+1));

            vM[i][j] = (1.0/8) * ( R(v,h,w,i-1,j-1)+ R(v,h,w,i-1,j)+ R(v,h,w,i-1,j-1)+
                                     R(v,h,w,i-1,j+1)+ R(v,h,w,i,j-1)+ R(v,h,w,i,j+1)+
                                    R(v,h,w,i+1,j-1)+ R(v,h,w,i+1,j+1));


            fu= (Ix[i][j]*uM[i][j]+Iy[i][j]*vM[i][j]+It[i][j])/(alpha+Ix[i][j]*Ix[i][j]+Iy[i][j]*Iy[i][j]);
            fv= (Ix[i][j]*uM[i][j]+Iy[i][j]*vM[i][j]+It[i][j])/(alpha+Ix[i][j]*Ix[i][j]+Iy[i][j]*Iy[i][j]);
            u[i][j]=uM[i][j]-Ix[i][j]*fu;
            v[i][j]=vM[i][j]-Iy[i][j]*fv;

            }
        }
        n++;
    }
}
void horn_shunck (double ** A, double ** B , int dimX, int dimY)
{
    double**Ix=CreerMatrice(dimX,dimY);
    double**Iy=CreerMatrice(dimX,dimY);
    double**It=CreerMatrice(dimX,dimY);
    gradient(A,B,Ix,Iy,It,dimX,dimY);
    double** U=CreerMatrice(dimX,dimY);
    double** uM=CreerMatrice(dimX,dimY);
    double** V=CreerMatrice(dimX,dimY);
    double** vM=CreerMatrice(dimX,dimY);

    calcul_u_v(U,V,Ix,Iy,It,uM,vM,dimX,dimY,Nmax);
    printf("calcul des matrices : uM,vM,U,V \t Nombre d'iterations : %d\n",Nmax);
    //ecriture des fichiers U,V:
    WriteTable(U,dimX,dimY,File_U);
    printf("Ecriture de la matrice U dans un fchier : U.txt\n");
    WriteTable(V,dimX,dimY,File_V);
    printf("Ecriture de la matrice V dans un fchier : V.txt\n\n");
    DetruireToutMatrices(Ix,Iy,It,U,V,uM,vM,dimX);

}
