#include "core/encseq_api.h"
#include "core/str_api.h"
#include "core/codetype.h"
#include "core/error.h"
#include "core/warning_api.h"
#include "core/alphabet_api.h"
#include "core/unused_api.h"

#include "match/sfx-mappedstr.h"

#include "stdbool.h"
#include "core/assert_api.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#include "score.h"
#include "scorematrix.h"

#define SIZE 5

GtUword** new_tau(GtUword numofsequences, GtUword r, GtUword length)
{
  GtUword **tau,
          i,
          range;
  range = (pow(r, length)+1);
  tau = malloc((numofsequences)*sizeof(GtUword*));
  *tau = calloc((numofsequences),range*sizeof(GtUword));
  for (i = 1; i < numofsequences; i++)
    tau[i] = tau[i-1]+range;
  return tau;
}

void delete_tau(GtUword **tau)
{
  if (tau != NULL)
  {
    free(*tau);
    free(tau);
  }
}

GtUword min(GtUword a, GtUword b)
{
  if (a <= b)
    return (a);
  else
    return (b);
}

long max(long a, long b)
{
  if (a <= b)
    return (b);
  else
    return (a);
}

void gt_get_codes(const GtEncseq *encseq,
                  GtUword len,
                  GtUword **tau)
{
  GtUword i;
  GtKmercodeiterator *kc_iter;
  const GtKmercode *kmercode;
  GtUword seqnum,
          pos = 0;
  gt_assert(encseq != NULL);
  for (i = 0; i < len; i++)
  {
    if (gt_encseq_total_length(encseq) <= pos)
      return;
    if (gt_encseq_position_is_separator(encseq, pos, GT_READMODE_FORWARD))
    {
      gt_warning("Sequence "GT_WU" is shorter than given k.\n",
                 gt_encseq_seqnum(encseq, pos));
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
    if (gt_encseq_total_length(encseq) <= pos)
      break;
    if (!(kmercode->definedspecialposition))
    {
      seqnum = gt_encseq_seqnum(encseq, pos);
      tau[seqnum][kmercode->code]++;
    }
  }
  gt_kmercodeiterator_delete(kc_iter);
}

GtWord** new_table(GtUword seqlength_first, GtUword seqlength_second)
{
  GtUword i;
  GtWord **table;

  table = malloc((seqlength_first+1)*sizeof(GtWord*));
  *table = malloc((seqlength_first+1)*(seqlength_second+1)*sizeof(GtWord));
  for (i=1; i <= seqlength_first; i++)
      table[i] = table[i-1]+(seqlength_second+1);

  return table;
}

void delete_table(GtWord **table)
{
  if (table != NULL)
  {
    free(*table);
    free(table);
  }
}

Score *calc_qgram(GtEncseq *encseq_first,
                           GtEncseq *encseq_second,
                           GtUword r,
                           GtUword q,
                           GtError *err)
{
  gt_error_check(err);
  gt_assert(encseq_first != NULL);
  gt_assert(r > 0 && q > 0);

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

  if (encseq_first && encseq_second)
  {
    numofseqsecond = gt_encseq_num_of_sequences(encseq_second);
    tv = new_tau(numofseqsecond, r, q);
    gt_get_codes(encseq_second, q, tv);

    score = malloc(sizeof(Score)*((numofseqfirst*numofseqsecond)+1));
    score->pos = 0;

    for (i = 0; i < numofseqfirst; i++)
    {
      for (j = 0; j < numofseqsecond; j++)
      {
        dist = 0;
        for (l = 0; l < pow(r,q); l++)
        {
          if (tu[i][l] > 0 || tv[j][l] > 0)
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

    for (i = 0; i < numofseqfirst; i++)
    {
      for (j = (i+1); j < numofseqfirst; j++)
      {
        dist = 0;
        for (l = 0; l < pow(r,q); l++)
        {
          if (tu[i][l] > 0 || tu[j][l] > 0)
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
  gt_assert(encseq_first != NULL);
  gt_assert(r > 0 && k > 0);

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

  if (encseq_first && encseq_second)
  {
    numofseqsecond = gt_encseq_num_of_sequences(encseq_second);
    tv = new_tau(numofseqsecond, r, k);
    gt_get_codes(encseq_second, k, tv);

    score = malloc(sizeof(Score)*((numofseqfirst*numofseqsecond)+1));
    score->pos = 0;

    for (i = 0; i < numofseqfirst; i++)
    {
      for (j = 0; j < numofseqsecond; j++)
      {
        tmp = 0;
        for (l = 0; l < pow(r,k); l++)
        {
          if (tu[i][l] > 0 && tv[j][l] > 0)
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

    for (i = 0; i < numofseqfirst; i++)
    {
      for (j = (i+1); j < numofseqfirst; j++)
      {
        tmp = 0;
        for (l = 0; l < pow(r,k); l++)
        {
          if (tu[i][l] > 0 && tu[j][l] > 0)
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

Score *calc_edist(GtEncseq *encseq_first,
                  GtEncseq *encseq_second,
                  GtStr *scorematrix,
                  GT_UNUSED GtError *err)
{
  Scorematrix *smatrix;
  FILE *fp;
  GtUword numofseqfirst,
          numofseqsecond,
          m,
          n;
  GtUword i, j, k, l,
          tmp,
          maximum,
          startone,
          starttwo;
  GtWord ins, del, rep,
         **table;
  Score *score;

  fp = fopen(gt_str_get(scorematrix), "r");
  smatrix = read_score(fp);
  fclose(fp);

  numofseqfirst = gt_encseq_num_of_sequences(encseq_first);
  numofseqsecond = gt_encseq_num_of_sequences(encseq_second);

  score = malloc(sizeof(Score)*((numofseqfirst*numofseqsecond)+1));
  score->pos = 0;

  for (i = 0; i < numofseqfirst; i++)
  {
    for (j = 0; j < numofseqsecond; j++)
    {
      startone = gt_encseq_seqstartpos(encseq_first,i);
      starttwo = gt_encseq_seqstartpos(encseq_second,j);
      m = gt_encseq_seqlength(encseq_first, i);
      n = gt_encseq_seqlength(encseq_second, j);

      table = new_table(m, n);
      table[0][0] = 0;
      for (k = 1; k <= m; k++)
        table[k][0] = table[k-1][0] + 1;

      for (l = 1; l <= n; l++)
      {
        table[0][l] = table[0][l-1] + 1;
        for (k = 1; k <= m; k++)
        {
          ins = table[k][l-1] + 1;
          del = table[k-1][l] + 1;
          char one = gt_encseq_get_decoded_char(encseq_first,
                                                k+startone-1,
                                                GT_READMODE_FORWARD);
          char two = gt_encseq_get_decoded_char(encseq_second,
                                                l+starttwo-1,
                                                GT_READMODE_FORWARD);
          rep = table[k-1][l-1] + access_scorematrix(smatrix, one, two);
          tmp = max(ins, del);
          maximum = max(tmp, rep);
          table[k][l] = maximum;
        }
      }
      score[score->pos].dist = table[m][n];
      score[score->pos].seqnum_u = i;
      score[score->pos].seqnum_v = j;
      score->pos++;
    }
  }
  delete_scorematrix(smatrix);
  return score;
}