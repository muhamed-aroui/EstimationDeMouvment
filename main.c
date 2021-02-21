#include <stdio.h>
#include <stdlib.h>
#include "myBmpGris.h"
#include "utilitaire.h"
#include "horn_schunk.h"
#include "Gauss.h"
#define IMAGE1 "image1.bmp"
#define IMAGE2 "image2.bmp"
int main()
{
    int i,j;
    double** A;
    double** B;
    //--------- lire image 1----------
    BmpImg img1;
    img1 = readBmpImage(IMAGE1);
    // ------- lire image 2  -----
    BmpImg img2;
    img2 = readBmpImage(IMAGE2);
    A= bmpToMatrice(img1);
    B= bmpToMatrice(img2);
    int dimX =img1.dimX;
    int dimY =img1.dimY;
    int dimX2 =img2.dimX;
    int dimY2 =img2.dimY;
    printf("/*\tP6 ESTIMATION DU MOUVEMENT\t*/ \n ");
    printf("1- METHODE 1 : \t Horn & Schunck\n"); // METHODE DE HORN & SCHUNCK:
    horn_shunck(A,B ,dimX,dimY);
    printf("2- METHODE 2 : \t Methode de multi-resolution\n\n"); // METHODE DE MULTI-RESOLUTION :
    multi_resolution(A,B,dimX,dimY,dimX2,dimY2);
    freeBmpImg(&img1); // détruire img1
    freeBmpImg(&img2); // détruire img2
    // NB: nous avons détruit tout les matrices utilisées dans les fonctions horn_shunck et multi_resolution
    return 0;
}
