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

GtUword *def_factor(GtUword r, GtUword q)
{
    GtUword i, 
            *factor;
    factor=calloc(sizeof(GtUword),(q+1));
    factor[0]=1;
    for(i = 1; i < q; i++)
    {
        factor[i]=factor[i-1]*r;
    }
    return(factor);
}

GtUword *alphabetcode(GtAlphabet *alpha, GtError *err)
{
  gt_error_check(err);
  assert(alpha != NULL);
  GtUword i;
  GtUword *alpha_tab;
  GtUword r;
  const GtUchar* alphabet;
  bool haserr = false;
  
  alpha_tab = malloc(sizeof(GtUword)*UCHAR_MAX);
  r = gt_alphabet_size(alpha);
  alphabet = gt_alphabet_characters(alpha);

  for(i=0; i<UCHAR_MAX; i++)
    alpha_tab[i] = -1;
  for(i=0; i < (r-1); i++)
  {
    if(alpha_tab[(GtUword)alphabet[i]] == -1)
      alpha_tab[(GtUword)alphabet[i]] = i;      
    else
    {
      gt_error_set(err,"The same symbol %c occured more than" 
        "once in the given alphabet\n", alphabet[i]);
      haserr = true;
    }
  }
  GtUchar wildcard = gt_alphabet_wildcard_show(alpha);
  alpha_tab[wildcard] = i;
  return((haserr)? NULL : alpha_tab);
}

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
      seqnum = gt_encseq_seqnum(encseq, pos);
      tau[seqnum][kmercode->code]++;
    }
  }
  gt_kmercodeiterator_delete(kc_iter);
}
//Noch betrachtet, dass sequenz kleiner sein kann als q.
void gt_get_qgramcodes(const GtEncseq *encseq, 
                       GtUword r, 
                       GtUword q,
                       GtUword **tau,
                       GtError *err)
{
  gt_assert(encseq != NULL);
  GtUword *factor,
          *alphatab,
          seqnum,
          code = 0,
          i, j,
          tmp,
          delete;
          
  alphatab = alphabetcode(gt_encseq_alphabet(encseq), err);
  gt_error_check(err);
  factor = def_factor(r, q);
                                            
  for(i = 0; i < q; i++)
  {
    tmp = (GtUword) gt_encseq_get_decoded_char(encseq, i, 
                                               GT_READMODE_FORWARD);
    code += (alphatab[tmp]*factor[q-1-i]);
  }
  seqnum = gt_encseq_seqnum(encseq, i);
  tau[seqnum][code]+= 1;
  
  for(i = 1; i <= (gt_encseq_total_length(encseq)-q); i++)
  {
    if(gt_encseq_position_is_separator(encseq, i+q-1, GT_READMODE_FORWARD))
    {
      code = 0;
      i = (i + q);
      for(j = 0; j < q; j++)
      { 
        tmp = (GtUword) gt_encseq_get_decoded_char(encseq, i, 
                                                   GT_READMODE_FORWARD);
        code += (alphatab[tmp]*factor[q-1-j]);
        i++;
      }
      seqnum = gt_encseq_seqnum(encseq, i-1);
      tau[seqnum][code]+= 1;
      i = (i - (q - 1));
    }
    

    tmp = (GtUword) gt_encseq_get_decoded_char(encseq, i+q-1, 
                                               GT_READMODE_FORWARD);
    delete = (GtUword) gt_encseq_get_decoded_char(encseq, i-1, 
                                                  GT_READMODE_FORWARD);
    code = (code - alphatab[delete]*factor[q-1])*r + alphatab[tmp];
    seqnum = gt_encseq_seqnum(encseq, i);
    tau[seqnum][code]+= 1;
  }
  free(factor);
  free(alphatab);
}       

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
          i, j, l,
          dist = 0;
  Score *score;
  
  numofseqfirst = gt_encseq_num_of_sequences(encseq_first);
  tu = new_tau(numofseqfirst, r, q);
  gt_get_qgramcodes(encseq_first, r, q, tu, err);
  
  if(encseq_first && encseq_second)
  {
    numofseqsecond = gt_encseq_num_of_sequences(encseq_second);
    tv = new_tau(numofseqsecond, r, q);
    gt_get_qgramcodes(encseq_second, r, q, tv, err);

    score = malloc(sizeof(Score)*((numofseqfirst*numofseqsecond)));
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
    score = malloc(sizeof(Score)*((numofseqfirst*numofseqfirst)));
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
  gt_get_kmercodes(encseq_first, k, tu);
  
  if(encseq_first && encseq_second)
  {  
    numofseqsecond = gt_encseq_num_of_sequences(encseq_second);
    tv = new_tau(numofseqsecond, r, k);
    gt_get_kmercodes(encseq_second, k, tv);
    
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
