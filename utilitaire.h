#ifndef _utilitaire_h_
#define _utilitaire_h_
#include <stdlib.h>
#include <stdio.h>
#define File_U "U.txt"
#define File_V "V.txt"
#define File_U_pyr "U_pyr.txt"
#define File_V_pyr "V_pyr.txt"
#include "myBmpGris.h"
double** CreerMatrice(int x,int y);

void DetruireMatrice(double** matrice, int x);

double** bmpToMatrice(BmpImg bmpImg);

void printMatrice(double** matrice, int dimX, int dimY);

void WriteTable(double** MatToBeWritten,int height, int width,char* path);

void DetruireToutMatrices (double**mat1,double**mat2,double**mat3,double**mat4,double**mat5,double**mat6,double**mat7,int DIM);
#endif
