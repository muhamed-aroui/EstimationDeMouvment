#include "utilitaire.h"
/** Methode de creation d'une matrice des flottans
*/
double** CreerMatrice(int x,int y){
    int i;
    double** matrice=malloc(x * sizeof(double*));
    for(i=0;i<x;i++)
        matrice[i]=calloc(y, sizeof(double));

    return matrice;
}

/** Methode de destruction d'une matrice des flottans
*/
void DetruireMatrice(double** matrice, int x){
    int i;
    for(i=0;i<x;i++)
        free(matrice[i]);
    free(matrice); matrice=NULL;
}

double** bmpToMatrice(BmpImg bmpImg){
    int i, j;
    double** res = CreerMatrice(bmpImg.dimX, bmpImg.dimY); //initialisation résultat
    for(i = 0; i < bmpImg.dimX; i++){
        for(j = 0; j < bmpImg.dimY; j++){
            res[i][j] = 1.0*getPixel(bmpImg, i, j);         //recuperation valeur pixel aux coord correspondantes
        }
    }
    return res;
}


void WriteTable(double** MatToBeWritten,int height, int width,char* path)
{

    if(height%2!=0)
        height--;
    if(width%2!=0)
        width--;

    FILE * fTxt = fopen ( path , "w" );

    if( fTxt != NULL )
    {
        int i = 0,j = 0;
        for(i=height-1;i>=0;i--)
        {
            for(j=0;j<width;j++)
            {
                fprintf( fTxt , "%f ",MatToBeWritten [i][j]);
            }
            fprintf( fTxt , "\n");
        }



        fclose ( fTxt );
    }

    }
   void DetruireToutMatrices (double**mat1,double**mat2,double**mat3,double**mat4,double**mat5,double**mat6,double**mat7,int DIM)
   {
       DetruireMatrice(mat1,DIM);
       DetruireMatrice(mat2,DIM);
       DetruireMatrice(mat3,DIM);
       DetruireMatrice(mat4,DIM);
       DetruireMatrice(mat5,DIM);
       DetruireMatrice(mat6,DIM);
       DetruireMatrice(mat7,DIM);
   }
