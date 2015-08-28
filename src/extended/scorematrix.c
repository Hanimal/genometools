
#include "scorematrix.h"
#define SIZE 20


/* allokiert Speicherplatz fuer eine obere Dreiecksmatrix*/
void allocate_matrix(Scorematrix *smatrix)
{
    unsigned int i;
    smatrix->matrix=malloc(sizeof(int*)*smatrix->dim);
    *(smatrix->matrix)=malloc(((smatrix->dim*(smatrix->dim+1))>>1)*sizeof(int));

    for(i=1; i<smatrix->dim; i++)
    {
        smatrix->matrix[i]=smatrix->matrix[i-1]+(smatrix->dim-(i-1));
    }
}

/* liest die Scorematrix aus einer Datei ein */
Scorematrix *read_score(FILE *fp)
{
    Scorematrix *smatrix=malloc(sizeof(Scorematrix));
    int temp;
    smatrix->dim=0;
    unsigned int j,i=0, size=SIZE;
    char c;
    char *s;
    
    s=malloc(sizeof(char)*SIZE); 
    CHECK_READ(fread(&c, sizeof(char),1, fp),1);
    while (c=='#')/*Kommentarzeile der Datei lesen*/
    {
        while (c!='\n')
        {
            CHECK_READ(fread(&c, sizeof(char),1, fp),1);
        }
        CHECK_READ(fread(&c, sizeof(char),1, fp),1);
    }
    while(c!='\n')/*Aminsoaeuren-Reihenfolge erfahren*/
    {
        if(i==size)
        {
            size=size+SIZE;
            s=realloc(s, size);
        }
        if(c!=' ')
        {
            s[i]=c;
            i++;
        }
        CHECK_READ(fread(&c, sizeof(char),1, fp),1);
    }
    if(i==size)
    {
        size+=1;
        s=realloc(s, size);
    }
    s[i]='\0';
    smatrix->order=s;
    smatrix->dim=i;
    allocate_matrix(smatrix);
    
    /*Score-Werte einlesen*/
    for(i=0; i<smatrix->dim; i++)
    {
        CHECK_READ(fscanf(fp,"%c",&c),1);
        if(c != smatrix->order[i])
        {
            fprintf(stderr, "Datei hat falsches Format: "
                    "Die Reihenfolge der Aminosaueren stimmt nicht, "
                    "Matrix ist nicht symmetrisch\n");
            return(NULL);
        }
        for(j=0;j<i;j++)
        {
            CHECK_READ(fscanf(fp," %d", &temp),1);
            if(temp != smatrix->matrix[j][i-j])
            {
                fprintf(stderr, "Datei hat falsches Format: "
                                "Matrix ist nicht symmetrisch\n");
                return(NULL);
            }
        } 
        for(; j<smatrix->dim; j++)
        {
            if(fscanf(fp," %d", &smatrix->matrix[i][j-i])!=1)
            {
                fprintf(stderr, "Datei hat falsches Format: Die Zeilen "
                                "sind nicht gleich lang\n");
                return(NULL);
            }
        }
        CHECK_READ(fscanf(fp,"%c",&c),1);
        while(c!='\n' || feof(fp)==1)
        {
            if(c!=' ')
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

/* 
 * Zugriff auf den Score fuer eine Ersetzungsoperation von zwei 
 * uebergebenen Aminosaeuren
 * 
 * bei Fehlern wie nicht existierender Aminsoaeuren wird LONG_MAX 
 * zurueck gegeben
 */
long access_scorematrix(Scorematrix *smatrix, char aminoacid1, char aminoacid2)
{
    int i, pos1=smatrix->dim, pos2=smatrix->dim;
    long value;
    assert(smatrix != NULL);

    /*Position der Aminosaeuren bestimmen*/
    for(i=0; i<smatrix->dim; i++)
    {
        if(smatrix->order[i]==aminoacid1)
            pos1=i;
        if(smatrix->order[i]==aminoacid2)
            pos2=i;
        if(pos1!=smatrix->dim && pos2!=smatrix->dim)
            break;
    }
    if(pos1==smatrix->dim)
    {
        fprintf(stderr, "Aminosaeure: %c existiert in uebergegebner "
                        "Scorematrix nicht\n", aminoacid1);
        return(LONG_MAX);
    }
    else if(pos2==smatrix->dim)
    {
        fprintf(stderr, "Aminosaeure: %c existiert in uebergegebner "
                        "Scorematrix nicht\n", aminoacid2);
        return(LONG_MAX);
    }
    /* Scorewert auslesen*/
    if(pos1<pos2)
        value=smatrix->matrix[pos1][pos2-pos1];
    else
        value=smatrix->matrix[pos2][pos1-pos2];
        
    return(value);
}

/* 
 * gibt eine die gesamte Matrix einer uebergebenen Scorematrix, die als
 * obere Dreiecksmatrix vorliegt, auf stdout aus
 * 
 * Anmerkung: wird hier nur zum testen vom access_scorematrix verwendet, 
 * daher diese Implementierung und nicht einfach eine doppelte 
 * for-Schleife die ueber die smatrix->matrix laeuft
 */
int show_scorematrix(Scorematrix *smatrix)
{
    unsigned i,j;
    long score;
    assert(smatrix != NULL);
    printf("Die gespeicherte Scorematrix enthaelt folgende Werte:\n");
    printf("  ");
    for(i=0; i<smatrix->dim; i++)
    {
        printf("%4c", smatrix->order[i]);
    }
    printf("\n");
    for(i=0; i<smatrix->dim; i++)
    {
        printf("%c ", smatrix->order[i]);
        for(j=0; j<smatrix->dim; j++)
        {
            score=access_scorematrix(smatrix, smatrix->order[i], smatrix->order[j]);
            if(score == LONG_MAX)
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
    if(smatrix!=NULL)
    {
        free(*(smatrix->matrix));
        free(smatrix->matrix);
        free(smatrix->order);
        free(smatrix);
    }
}

//int main(int argc, char **argv)
//{
    //FILE *fp;
    //Scorematrix *smatrix;
    //int control;
    
    //if(argc==2)
    //{
        //fp=fopen(argv[1], "r");
        //if(fp==NULL)
            //{
                //fprintf(stderr, "Dateiname '%s' nicht bekannt\n", argv[1]);
                //return(EXIT_FAILURE);
            //}
        //smatrix=read_score(fp);
        //fclose(fp);
        //if(smatrix == NULL)
            //return(EXIT_FAILURE);
        //control=show_scorematrix(smatrix);/*zum testen von access-Funktion*/
        //delete_scorematrix(smatrix);
        //if(control == EXIT_FAILURE)
            //return(EXIT_FAILURE);
    //}
    //else
    //{
        //fprintf(stderr, "Falsche Anzahl an Parametern\n"
                        //"USAGE: %s <dateiname>\n", argv[0]);
        //return(EXIT_FAILURE);
    //}
    //return(EXIT_SUCCESS);
//}
