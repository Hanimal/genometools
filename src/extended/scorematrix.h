#ifndef SCOREMATRIX_H
#define SCOREMATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <limits.h>
#define CHECK_READ(READ,N) if(READ != N){fprintf(stderr,"Fehler beim "\
                                         "Lesen der Datei\n");return(NULL);}
                                         
                                         
typedef struct{
    int **matrix;
    unsigned int dim;
    char *order;
} Scorematrix;

Scorematrix *read_score(FILE *fp);

void allocate_matrix(Scorematrix *smatrix);

Scorematrix *read_score(FILE *fp);

long access_scorematrix(Scorematrix *smatrix, char aminoacid1, char aminoacid2);

int show_scorematrix(Scorematrix *smatrix);

void delete_scorematrix(Scorematrix *smatrix);

#endif