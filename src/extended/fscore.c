#include "core/encseq_api.h"
#include "core/str_api.h"
#include "core/codetype.h"
#include "core/error.h"
#include "core/warning_api.h"
#include "core/alphabet_api.h"

#include "match/sfx-mappedstr.h"

#include "stdbool.h"
#include "assert.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#include "fscore.h"
#define SIZE 5

void gt_get_kmercodes(const GtEncseq *encseq, 
                         GtUword kmerlen,
                         GtUword **tau)
{
  GtKmercodeiterator *kc_iter;
  const GtKmercode *kmercode;
  GtUword seqnum,
		  pos;
  gt_assert(encseq != NULL);
  kc_iter = gt_kmercodeiterator_encseq_new(encseq, 
                                           GT_READMODE_FORWARD, 
                                           kmerlen,
                                           0);                                         
  while ((kmercode = gt_kmercodeiterator_encseq_next(kc_iter)) != NULL) 
  {
	
	pos = gt_kmercodeiterator_encseq_get_currentpos(kc_iter);
	if(gt_encseq_total_length(encseq) <= pos)
	  break;
    if (!(kmercode->definedspecialposition))
    {  
	  //if(gt_encseq_position_is_separator(encseq, pos, GT_READMODE_FORWARD))
		//break;
	  seqnum = gt_encseq_seqnum(encseq, pos);
      tau[seqnum][kmercode->code]++;
	}
  }
  gt_kmercodeiterator_delete(kc_iter);
}

GtUword** new_tau(GtUword numofsequences, GtUword r, GtUword k)
{
  GtUword **tau,
          i,
          range;
  range = (pow(r, k)+1);
  tau = malloc((numofsequences)*sizeof(GtUword*));
  *tau = calloc((numofsequences),range*sizeof(GtUword));  
  for(i=1; i < numofsequences; i++)
    tau[i] = tau[i-1]+range;
  return tau;
}

void delete_tau(GtUword **tau)
{
	if(tau != NULL)
	{
		free(*tau);
		free(tau);
	}
}

GtUword min(GtUword a, GtUword b)
{
  if(a <= b)
    return (a);
  else
    return (b);
}

FScore *fscore(GtEncseq *encseq_first, 
                GtEncseq *encseq_second, 
                GtUword r, 
                GtUword k, 
                GtError *err)
{
  gt_error_check(err);
  assert(encseq_first != NULL && encseq_second != NULL);
  assert(r > 0 && k > 0);
    
  GtUword numofseqfirst;
  GtUword numofseqsecond;
  GtUword **tu,
          **tv,
          i, j, l,
          tmp = 0;
  float dist = 0;
  FScore *score;
    
  numofseqfirst = gt_encseq_num_of_sequences(encseq_first);
  numofseqsecond = gt_encseq_num_of_sequences(encseq_second);
  tu = new_tau(numofseqfirst, r, k);
  tv = new_tau(numofseqsecond, r, k);
  
  score = malloc(sizeof(FScore)*((numofseqfirst*numofseqsecond)+1));
  score->pos = 0;
  
  gt_get_kmercodes(encseq_first, k, tu);
  gt_get_kmercodes(encseq_second, k, tv);
  
  for(i = 0; i < numofseqfirst; i++)
  {
    for(j = 0; j < numofseqsecond; j++)
    {
      tmp = 0;
      for(l = 0; l < pow(r,k); l++)
      {
        if(tu[i][l] == tv[j][l] && tu[i][l] > 0 && tv[j][l] > 0)
          tmp += min(tu[i][l],tv[j][l]);
      }
      GtUword length_u = gt_encseq_seqlength(encseq_first, i);
      GtUword length_v = gt_encseq_seqlength(encseq_second, j);
      dist = (float)tmp/(float)(min(length_u, length_v)-k+1);
      score[score->pos].dist = dist;
      score[score->pos].seqnum_u = i;
      score[score->pos].seqnum_v = j;
      score->pos++;
    }
  }
  return score;
}
