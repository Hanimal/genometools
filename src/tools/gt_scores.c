#include "core/encseq_api.h"
#include "core/str_api.h"
#include "assert.h"
#include "core/error.h"
#include "stdbool.h"
#include "core/warning_api.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#include "extended/fscore.h"
#include "gt_scores.h"

typedef struct{
    GtUword score,
            seqnum_u,
            seqnum_v;
}Fscore;


static GtEncseq *gt_encseq_get_encseq(const char *seqfile,
                                             GtError *err)
{
    GtEncseqLoader *encseq_loader;
    GtEncseq *encseq;
    int had_err = 0;
    gt_error_check(err);
    gt_assert(seqfile);
  
    encseq_loader = gt_encseq_loader_new();
    if (!(encseq = gt_encseq_loader_load(encseq_loader, seqfile, err)))
      had_err = -1;

    gt_encseq_loader_delete(encseq_loader);
  
    if(!had_err) 
    {
      if (!gt_encseq_has_description_support(encseq))
        gt_warning("Missing description support for file %s", seqfile);
        return encseq;
    } 
    else
      return NULL;
}


int gt_fscore(int argc, const char **argv, GtError *err);
{
  if(argc != 4)
  {
    printf("USAGE: %s <alphabetfile> <sequencefile> <k>", argv[0]);
    return(EXIT_FAILURE);
  }
  else
  {
    GtEncseq *encseq;
    GtAlphabet *alpha; 
    Fscore *fscore;
    GtUword *alpha_code;
    GtUword k, 
      i,
      j,
      file_num,
      r,
      idx,
      struct_size = 0,
      size = 10;
  
    sscanf(argv[3], GT_WU, &k);
    alpha_code = malloc(sizeof(GtUword)*UCHAR_MAX);
    fscore = malloc(sizeof(Fscore)*size);
  
    //err = gt_error_new();
    //gt_error_set_progname(err, argv[0]);
    encseq = gt_encseq_get_encseq(argv[1], err);
    alpha = gt_alphabet_new_from_file_no_suffix(argv[1], err);
    gt_error_check(err);
    r = gt_alphabet_size(alpha);
    match_alphabetcode(alpha, alpha_code, err);
    gt_error_check(err);
  
    file_num = gt_encseq_num_of_files(encseq);
    idx = 0;
    for(i = 0; i < file_num; i++)
    {
      for(j = i+1; j < file_num; j++)
      {
        /*GtUword score = f_score(encseq, 
                    i, 
                    j,
                    r, 
                    k, 
                    alpha_code, 
                    err);
        if(struct_size == size)
        {
          size += 10;
          fscore = realloc(fscore, size*sizeof(Fscore));
        }
        fscore[idx].score = score;
        fscore[idx].seqnum_u = i;
        fscore[idx].seqnum_v = j;
        idx++;
        struct_size++;*/
      }
    }
    gt_encseq_delete(encseq);
    return(EXIT_SUCCESS);
  }
}
