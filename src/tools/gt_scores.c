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


int gt_scores(int argc, const char **argv, GtError *err)
{
  if(argc != 5)
  {
    printf("USAGE: %s <sequencefile> <sequencefile> <alphabetsize> <k>\n", argv[0]);
    return(EXIT_FAILURE);
  }
  else
  {
	gt_error_check(err);
    GtEncseq *encseq_first;
    GtEncseq *encseq_second;
    FScore *score;
    GtUword k, 
			r,
			i;
  
    if((sscanf(argv[4], GT_WU, &k) != 1) || k < 1)
    {
      printf("%s is no valid input for k\n", argv[4]);
      return(EXIT_FAILURE);
	}
	if((sscanf(argv[3], GT_WU, &r) != 1) || r < 1)
    {
      printf("%s is no valid input for the alphabetsize\n", argv[3]);
      return(EXIT_FAILURE);
	}
  
    encseq_first = gt_encseq_get_encseq(argv[1], err);
    gt_error_check(err);
    encseq_second = gt_encseq_get_encseq(argv[2], err);
    gt_error_check(err);
  
    score = fscore(encseq_first, 
                    encseq_second, 
					r, 
					k, 
					err);
	
	for( i = 0; i < score->pos; i++)
	{
		printf("Fscore between sequence "GT_WU" and "GT_WU" is %.3f.\n", 
				score[i].seqnum_u, score[i].seqnum_v, score[i].dist);
	}
	
    gt_encseq_delete(encseq_first);
    gt_encseq_delete(encseq_second);
    free(score);
    return(EXIT_SUCCESS);
  }
}
