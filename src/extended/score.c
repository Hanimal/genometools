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

#include "score.h"
#define SIZE 5

GtUword** new_tau(GtUword numofsequences, GtUword r, GtUword length)
{
  GtUword **tau,
          i,
          range;
  range = (pow(r, length)+1);
  tau = malloc((numofsequences)*sizeof(GtUword*));
  *tau = calloc((numofsequences),range*sizeof(GtUword));  
  for(i = 1; i < numofsequences; i++)
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

void gt_get_codes(const GtEncseq *encseq, 
                  GtUword len,
                  GtUword **tau)
{
  GtUword i;
  /*printf("\n Sequenz: "GT_WU"\n", gt_encseq_seqnum(encseq, 0));
  for(i = 0; i < (gt_encseq_total_length(encseq)); i++)
  {
    if(gt_encseq_position_is_separator(encseq, i, GT_READMODE_FORWARD))
    {
      printf("\n");
      printf("Sequenz: "GT_WU"\n", gt_encseq_seqnum(encseq, i+1));
      continue;
    }
    printf("%c",gt_encseq_get_decoded_char(encseq, i, GT_READMODE_FORWARD));
  }
  printf("\n");*/
  GtKmercodeiterator *kc_iter;
  const GtKmercode *kmercode;
  GtUword seqnum,
          pos = 0;
  gt_assert(encseq != NULL);
  for(i = 0; i < len; i++)
  {
    if(gt_encseq_total_length(encseq) <= pos)
      return;
    if(gt_encseq_position_is_separator(encseq, pos, GT_READMODE_FORWARD))
    {
      gt_warning("Sequence "GT_WU" is shorter than given k.\n", gt_encseq_seqnum(encseq, pos));
      i = 0;
      pos++;
      continue;
    }
    pos++;
  }
  kc_iter = gt_kmercodeiterator_encseq_new(encseq, 
                                           GT_READMODE_FORWARD, 
                                           len,
                                           pos-len);                                         
  while ((kmercode = gt_kmercodeiterator_encseq_next(kc_iter)) != NULL) 
  {
	  pos = (gt_kmercodeiterator_encseq_get_currentpos(kc_iter)-len); 
	  if(gt_encseq_total_length(encseq) <= pos)
	    break;
    if (!(kmercode->definedspecialposition))
    {  
      seqnum = gt_encseq_seqnum(encseq, pos);
      tau[seqnum][kmercode->code]++;
    }
  }
  gt_kmercodeiterator_delete(kc_iter);
}

Score *calc_qgram(GtEncseq *encseq_first, 
                           GtEncseq *encseq_second, 
                           GtUword r, 
                           GtUword q, 
                           GtError *err)
{
  gt_error_check(err);
  assert(encseq_first != NULL);
  assert(r > 0 && q > 0);
    
  GtUword numofseqfirst,
          numofseqsecond;
  GtUword **tu,
          **tv,
          i, j, l;
  float dist = 0;
  Score *score;
  
  numofseqfirst = gt_encseq_num_of_sequences(encseq_first);
  tu = new_tau(numofseqfirst, r, q);
  gt_get_codes(encseq_first, q, tu);
  
  if(encseq_first && encseq_second)
  {  
    numofseqsecond = gt_encseq_num_of_sequences(encseq_second);
    tv = new_tau(numofseqsecond, r, q);
    gt_get_codes(encseq_second, q, tv);
    
    score = malloc(sizeof(Score)*((numofseqfirst*numofseqsecond)+1));
    score->pos = 0;
    
    for(i = 0; i < numofseqfirst; i++)
    {
      for(j = 0; j < numofseqsecond; j++)
      {
        dist = 0;
        for(l = 0; l < pow(r,q); l++)
        {
          if(tu[i][l] > 0 || tv[j][l] > 0)
            dist += abs(tu[i][l] - tv[j][l]);
        }
        score[score->pos].dist = dist;
        score[score->pos].seqnum_u = i;
        score[score->pos].seqnum_v = j;
        score->pos++;
      }
    }
    delete_tau(tv);
  }
  else
  {
    score = malloc(sizeof(Score)*((numofseqfirst*numofseqfirst)+1));
    score->pos = 0;
  
    for(i = 0; i < numofseqfirst; i++)
    {
      for(j = (i+1); j < numofseqfirst; j++)
      {
        dist = 0;
        for(l = 0; l < pow(r,q); l++)
        {
          if(tu[i][l] > 0 || tu[j][l] > 0)
            dist += abs(tu[i][l] - tu[j][l]);
        }
        score[score->pos].dist = dist;
        score[score->pos].seqnum_u = i;
        score[score->pos].seqnum_v = j;
        score->pos++;
      }
    }
  }
  delete_tau(tu);
  return score;
}


Score *calc_fscore(GtEncseq *encseq_first, 
                   GtEncseq *encseq_second, 
                   GtUword r, 
                   GtUword k, 
                   GtError *err)
{
  gt_error_check(err);
  assert(encseq_first != NULL);
  assert(r > 0 && k > 0);
    
  GtUword numofseqfirst,
          numofseqsecond;
  GtUword **tu,
          **tv,
          i, j, l,
          tmp;
  float dist = 0;
  Score *score;

  numofseqfirst = gt_encseq_num_of_sequences(encseq_first);
  tu = new_tau(numofseqfirst, r, k);
  gt_get_codes(encseq_first, k, tu);
  
  if(encseq_first && encseq_second)
  {  
    numofseqsecond = gt_encseq_num_of_sequences(encseq_second);
    tv = new_tau(numofseqsecond, r, k);
    gt_get_codes(encseq_second, k, tv);
    
    score = malloc(sizeof(Score)*((numofseqfirst*numofseqsecond)+1));
    score->pos = 0;
    
    for(i = 0; i < numofseqfirst; i++)
    {
      for(j = 0; j < numofseqsecond; j++)
      {
        tmp = 0;
        for(l = 0; l < pow(r,k); l++)
        {
          if(tu[i][l] > 0 && tv[j][l] > 0)
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
    delete_tau(tv);
  }
  else
  {
    score = malloc(sizeof(Score)*((numofseqfirst*numofseqfirst)+1));
    score->pos = 0;
  
    for(i = 0; i < numofseqfirst; i++)
    {
      for(j = (i+1); j < numofseqfirst; j++)
      {
        tmp = 0;
        for(l = 0; l < pow(r,k); l++)
        {
          if(tu[i][l] > 0 && tu[j][l] > 0)
            tmp += min(tu[i][l],tu[j][l]);
        }
        GtUword length_u = gt_encseq_seqlength(encseq_first, i);
        GtUword length_v = gt_encseq_seqlength(encseq_first, j);
        dist = (float)tmp/(float)(min(length_u, length_v)-k+1);
        score[score->pos].dist = dist;
        score[score->pos].seqnum_u = i;
        score[score->pos].seqnum_v = j;
        score->pos++;
      }
    }
  }
  delete_tau(tu);
  return score;
}
