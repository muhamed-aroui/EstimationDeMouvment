#include "Gauss.h"
#include "utilitaire.h"
#include "myBmpGris.h"
/** Methode de calcul des coifficients du filtre de GAUSS
*/
double GaussFiltre(int x,int y)
{
    double den = -2.*sigma*sigma;
    double num=((x*x)+(y*y))/den;
    double result=exp(num)*(1./den*pi);
    return result;
}
double ** GaussMatrix()
{
    int i,j;
    double somme=0;

    double ** matrice=CreerMatrice(gaussDim,gaussDim);
    for(i=0;i<gaussDim;i++)
    {
        for(j=0;j < gaussDim; j++)
        {
            matrice[i][j]=GaussFiltre(i-1,j-1);
            somme+=matrice[i][j];
        }
    }

    //normalisation de la matrice de gauss
    for(i=0;i<gaussDim;i++)
         for(j=0;j < gaussDim; j++)
            matrice[i][j]/=somme;

    return matrice;
}


/** Methode de produit de convolution
*/
double** Convlution(double** h, int dimX,int dimY,double** x)
{
    double** y=CreerMatrice(dimX,dimY);
    int i,j;
    for(i=1;i<dimX-1;i++)
    {
        for(j=1;j<dimY-1;j++)
        {
            y[i][j] =  h[i+1][j+1]*x[0][0] + h[i][j+1]*x[1][0] + h[i-1][j+1]*x[2][0]+
                            h[i+1][j]*x[0][1]   + h[i][j]*x[1][1]   + h[i-1][j]*x[2][1]+
                            h[i+1][j-1]*x[0][2] + h[i][j-1]*x[1][2] + h[i-1][j-1]*x[2][2];
        }
    }
    return y;
}


/** Methode de Produit des matrices terme par terme
*/
double** ProduitMatTermeParTerme(double** terme1,double** terme2,int dimX,int dimY )
{
    int i,j;
    double**result=CreerMatrice(dimX,dimY);
    for(i=0;i<dimX;i++)
        for(j=0;j<dimY;j++)
            result[i][j]=terme1[i][j]*terme2[i][j];

    return result;
}

/** Methode de calcul U et V
*/
void CalcU_V_avecConv(double**U,double**V,double** uM ,double ** vM, double** FG , double** Ix , double** Iy,double** It,int lignes,int cols)
{
    //recadrement de l'image
    if(lignes%2!=0)
        lignes--;
    if(cols%2!=0)
        cols--;

    double** IyIt=ProduitMatTermeParTerme(Iy,It,lignes,cols);
    double** IxIt=ProduitMatTermeParTerme(Ix,It,lignes,cols);
    double** IxIy=ProduitMatTermeParTerme(Ix,Iy,lignes,cols);
    double** IxIx=ProduitMatTermeParTerme(Ix,Ix,lignes,cols);
    double** IyIy=ProduitMatTermeParTerme(Iy,Iy,lignes,cols);

    double** temp=Convlution(IyIt,lignes,cols,FG);      // temp = Covolution [ (Iy*It) , FG ]
    double** temp1=Convlution(IxIt,lignes,cols,FG);     // temp1 = Covolution [ (Ix*It) , FG ]
    double** temp2=Convlution(IxIy,lignes,cols,FG);     // temp2 = Covolution [ (Ix*Iy) , FG ]
    double** temp3=Convlution(IxIx,lignes,cols,FG);        // temp3 = Covolution [ (Ix*Ix) , FG ]
    double** temp4=Convlution(IyIy,lignes,cols,FG);      // temp4 = Covolution [ (Iy*Iy) , FG ]

    int i,j;
    for(i=0;i<lignes;i++)
    {
        for(j=0;j<cols;j++)
        {
            //denomunateur
            //temp5 = alpha+Convolution [ (Ix*Ix) , FG]+ Convolution [ (Iy*Iy) , FG]
            double temp5=alpha_u+temp3[i][j]+temp4[i][j];
            //relation avec u
            //temp6= Covolution [ (Ix*It) , FG ]+ ( Covolution [ (Ix*Iy) , FG ] * vM)+ (Covolution [ (Ix*Ix) , FG ] * uM)
            double temp6=temp1[i][j]+(temp2[i][j]*vM[i][j])+(temp3[i][j]*uM[i][j]);
            //relation avec v
            //temp7= Covolution [ (Iy*It) , FG ]+ ( Covolution [ (Iy*Iy) , FG ] * vM)+ (Covolution [ (Ix*Iy) , FG ] * uM)
            double temp7=temp[i][j]+(temp4[i][j]*vM[i][j])+(temp2[i][j]*uM[i][j]);
            // formule de U
            U[i][j]=uM[i][j]-(temp6/temp5);
            // formule de V
            V[i][j]=vM[i][j]-(temp7/temp5);
        }
    }
}
/** Methode pour reduire des dimensions d'une matrice
*/
double **reduire_res(double** matIn,int dimX,int dimY)
{
    double** matOut=CreerMatrice(dimX/2,dimY/2);
    int i,j;
    for(i=0;i<dimX-2; i=i+2)
         for(j=0;j<dimY-2; j=j+2)
            matOut[i/2][j/2]=matIn[i][j];

    return matOut;
}

/** Methode de agrandissement des dimensions
*/
double** augment_res(double** matIn,int dimX,int dimY)
{
    double** matOut=CreerMatrice(dimX*2,dimY*2);
    int i,j;
    for(i=0;i<dimX-2;i=i+2)
        for(j=0;j<dimY-2;j=j+2)
            matOut[i*2][j*2]=matIn[i][j];

    return matOut;
}
multi_resolution (double **A,double **B, int dimX , int dimY,int dimX2,int dimY2)
{
    double** Gauss_M=GaussMatrix(); // FILTRE_GAUSSIEN_SIGMA_0.8
    double** img1_Lisse=A;
    img1_Lisse=Convlution(A,dimX,dimY,Gauss_M);
    printf("creation d'une Image 1 -SMOOTH - avec un filtre de GAUSS dim 3 - taille : %d x %d\n",dimX,dimY);

    double** img2_Lisse=B;
    img2_Lisse=Convlution(B,dimX2,dimY2,Gauss_M);
    printf("creation d'une Image 2 -SMOOTH - avec un filtre de GAUSS dim 3 - taille : %d x %d\n",dimX,dimY);

    // REDUIRE LA RESOLUTION  SCALE 2 :

    double** img1_LisseScale2=reduire_res(img1_Lisse,dimX, dimY);
    double** img2_LisseScale2=reduire_res(img2_Lisse, dimX2 , dimY2);

    double** img1_LisseScale2_Lisse=Convlution(img1_LisseScale2,dimX/2,dimY/2,Gauss_M);
    double** img2_LisseScale2_Lisse=Convlution(img2_LisseScale2,dimX2/2,dimY2/2,Gauss_M);

    // REDUIRE LA RESOLUTION  SCALE 4 :
    double** img1_LisseScale4=reduire_res(img1_LisseScale2_Lisse, dimX/2 ,dimY/2);
    double** img2_LisseScale4=reduire_res(img2_LisseScale2_Lisse , dimX2/2 , dimY2/2);

    // debut de calcul :
       double** u_4=CreerMatrice(dimX/4,dimY/4);
    double** uM_4=CreerMatrice(dimX/4,dimY/4);
    double** v_4=CreerMatrice(dimX/4,dimY/4);
    double** vM_4=CreerMatrice(dimX/4,dimY/4);
    double** Ix_4=CreerMatrice(dimX/4,dimY/4);
    double** Iy_4=CreerMatrice(dimX/4,dimY/4);
    double** It_4=CreerMatrice(dimX/4,dimY/4);

    gradient(img1_LisseScale4,img2_LisseScale4,Ix_4,Iy_4,It_4, dimX/4 ,dimY/4 );

    CalcU_V_avecConv(u_4,v_4,uM_4,vM_4,Gauss_M,Ix_4,Iy_4,It_4,dimX/4,dimY/4);
    double** u_2=augment_res(u_4,dimX/4,dimY/4);
    double** uM_2=augment_res(uM_4,dimX/4,dimY/4);
    double** v_2=augment_res(v_4,dimX/4,dimY/4);
    double** vM_2=augment_res(vM_4,dimX/4,dimY/4);
    double** Ix_2=augment_res(Ix_4,dimX/4,dimY/4);
    double** Iy_2=augment_res(Iy_4,dimX/4,dimY/4);
    double** It_2=augment_res(It_4,dimX/4,dimY/4);

    gradient(img1_LisseScale2,img2_LisseScale2, Ix_2 , Iy_2 , It_2 , dimX/2 , dimY/2);
    int i;
    for(i=0;i<N_pyr;i++)
        CalcU_V_avecConv( u_2, v_2, uM_2, vM_2, Gauss_M, Ix_2, Iy_2, It_2, dimX/2 , dimY/2);

    double** u_pyr=augment_res(u_2,dimX/2,dimY/2);
    double** uM_pyr=augment_res(uM_2,dimX/2,dimY/2);
    double** v_pyr=augment_res(v_2,dimX/2,dimY/2);
    double** vM_pyr=augment_res(vM_2,dimX/2,dimY/2);
    double** Ix_pyr=augment_res(Ix_2,dimX/2,dimY/2);
    double** Iy_pyr=augment_res(Iy_2,dimX/2,dimY/2);
    double** It_pyr=augment_res(It_2,dimX/2,dimY/2);

    gradient(img1_Lisse,img2_Lisse,Ix_pyr,Iy_pyr,It_pyr,dimX,dimY);

    for(i=0;i<N_pyr;i++){
        CalcU_V_avecConv(u_pyr,v_pyr,uM_pyr,vM_pyr,Gauss_M,Ix_pyr,Iy_pyr,It_pyr,dimX,dimY);
}

    WriteTable(u_pyr,dimX,dimY,File_U_pyr);
    printf("Ecriture de la matrice U dans un fchier : U_pyr.txt\n");
    WriteTable(v_pyr,dimX,dimY,File_V_pyr);
    printf("Ecriture de la matrice V dans un fchier : V_pyr.txt\n");
    printf("Liberation de la memoire\n ");
    DetruireMatrice(A,dimX);
    DetruireMatrice(B,dimX);
    DetruireToutMatrices(Ix_4,Iy_4,It_4,u_4,v_4,uM_4,vM_4,dimX/4);
    DetruireToutMatrices(Ix_2,Iy_2,It_2,u_2,v_2,uM_2,vM_2,dimX/2);
    DetruireToutMatrices(Ix_pyr,Iy_pyr,It_pyr,u_pyr,uM_pyr,v_pyr,vM_pyr,dimX);
}
