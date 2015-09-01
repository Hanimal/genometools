#include "core/fileutils.h"
#include "core/ma.h"

#include "scorematrix.h"
#define SIZE 20

void allocate_matrix(Scorematrix *smatrix)
{
    GtUword i;
    smatrix->matrix = malloc(sizeof(GtUword*)*smatrix->dim);
    *(smatrix->matrix) = malloc(((smatrix->dim*(smatrix->dim+1))>>1)
                         *sizeof(GtUword));

    for (i=1; i<smatrix->dim; i++)
        smatrix->matrix[i]=smatrix->matrix[i-1]+(smatrix->dim-(i-1));
}

Scorematrix *read_score(FILE *fp)
{
    Scorematrix *smatrix = malloc(sizeof(Scorematrix));
    int tmp;
    smatrix->dim=0;
    GtUword j,
            i = 0,
            size = SIZE;
    char c,
         *s;

    s = malloc(sizeof(char)*SIZE);
    CHECK_READ(fread(&c, sizeof char,1, fp),1);
    while (c == '#')
    {
        while (c != '\n')
            CHECK_READ(fread(&c, sizeof char,1, fp),1);
        CHECK_READ(fread(&c, sizeof char,1, fp),1);
    }
    while (c != '\n')
    {
        if (i == size)
        {
            size = size+SIZE;
            s = realloc(s, size);
        }
        if (c != ' ')
        {
            s[i] = c;
            i++;
        }
        CHECK_READ(fread(&c, sizeof char,1, fp),1);
    }
    if (i == size)
    {
        size+=1;
        s=realloc(s, size);
    }
    s[i] = '\0';
    smatrix->order = s;
    smatrix->dim = i;
    allocate_matrix(smatrix);

    for (i = 0; i < smatrix->dim; i++)
    {
        CHECK_READ(fscanf(fp,"%c",&c),1);
        if (c != smatrix->order[i])
        {
            fprintf(stderr, "Datei hat falsches Format: "
                    "Die Reihenfolge der Aminosaueren stimmt nicht, "
                    "Matrix ist nicht symmetrisch\n");
            return(NULL);
        }
        for (j = 0;j<i;j++)
        {
            CHECK_READ(fscanf(fp," %d", &tmp),1);
            if (tmp != smatrix->matrix[j][i-j])
            {
                fprintf(stderr, "Datei hat falsches Format: "
                                "Matrix ist nicht symmetrisch\n");
                return(NULL);
            }
        }
        for (; j < smatrix->dim; j++)
        {
            if (fscanf(fp," %d", &smatrix->matrix[i][j-i]) != 1)
            {
                fprintf(stderr, "Datei hat falsches Format: Die Zeilen "
                                "sind nicht gleich lang\n");
                return(NULL);
            }
        }
        CHECK_READ(fscanf(fp,"%c",&c),1);
        while (c != '\n' || feof(fp) == 1)
        {
            if (c != ' ')
            {
                fprintf(stderr, "Datei hat falsches Format: Die Zeilen "
                                "sind nicht gleich lang\n");
                return(NULL);
            }
            CHECK_READ(fscanf(fp,"%c",&c),1);
        }
    }
    return(smatrix);
}

GtWord access_scorematrix(Scorematrix *smatrix,
                          char aminoacid1,
                          char aminoacid2)
{
    int i,
        pos1 = smatrix->dim,
        pos2 = smatrix->dim;
    GtWord value;
    gt_assert(smatrix != NULL);

    for (i = 0; i < smatrix->dim; i++)
    {
        if (smatrix->order[i] == aminoacid1)
            pos1 = i;
        if (smatrix->order[i] == aminoacid2)
            pos2 = i;
        if (pos1 != smatrix->dim && pos2 != smatrix->dim)
            break;
    }
    if (pos1 == smatrix->dim)
    {
        printf("Aminosaeure: %c existiert in uebergegebner "
                        "Scorematrix nicht\n", aminoacid1);
        return(LONG_MAX);
    }
    else if (pos2 == smatrix->dim)
    {
        printf("Aminosaeure: %c existiert in uebergegebner "
                        "Scorematrix nicht\n", aminoacid2);
        return(LONG_MAX);
    }
    if (pos1 < pos2)
        value = smatrix->matrix[pos1][pos2-pos1];
    else
        value = smatrix->matrix[pos2][pos1-pos2];

    return(value);
}

int show_scorematrix(Scorematrix *smatrix)
{
    unsigned i,j;
    GtWord score;
    gt_assert(smatrix != NULL);
    printf("Die gespeicherte Scorematrix enthaelt folgende Werte:\n");
    printf("  ");
    for (i = 0; i < smatrix->dim; i++)
        printf("%4c", smatrix->order[i]);
    printf("\n");
    for (i = 0; i < smatrix->dim; i++)
    {
        printf("%c ", smatrix->order[i]);
        for (j = 0; j < smatrix->dim; j++)
        {
            score=access_scorematrix(smatrix,
                                     smatrix->order[i],
                                     smatrix->order[j]);
            if (score == LONG_MAX)
                return(EXIT_FAILURE);
            printf("%4ld", score);
        }
        printf(" \n");
    }
    return(EXIT_SUCCESS);
}

/*
 * gibt allokierten Speicher wieder frei
 */
void delete_scorematrix(Scorematrix *smatrix)
{
    if (smatrix != NULL)
    {
        free(*(smatrix->matrix));
        free(smatrix->matrix);
        free(smatrix->order);
        free(smatrix);
    }
}
