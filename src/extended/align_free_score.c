/*printf("\n Sequenz: "GT_WU"\n", gt_encseq_seqnum(encseq, 0));
for(i = 0; i < (gt_encseq_total_length(encseq)); i++)
{
  if (gt_encseq_position_is_separator(encseq, i, GT_READMODE_FORWARD))
  {
    printf("\n");
    printf("Sequenz: "GT_WU"\n", gt_encseq_seqnum(encseq, i+1));
    continue;
  }
  printf("%c",gt_encseq_get_decoded_char(encseq, i, GT_READMODE_FORWARD));
}
printf("\n");*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include "match/esa-splititv.h"
#include "match/sfx-mappedstr.h"
#include "match/esa-map.h"
#include "core/encseq_api.h"
#include "core/str_api.h"
#include "core/codetype.h"
#include "core/error.h"
#include "core/warning_api.h"
#include "core/alphabet_api.h"
#include "core/unused_api.h"
#include "core/types_api.h"
#include "core/seq_iterator_sequence_buffer_api.h"
#include "core/assert_api.h"
#include "core/minmax.h"
#include "core/array2dim_api.h"

#include "align_free_score.h"
#include "scorematrix.h"

typedef uint32_t GtKmercount;

GtKmercount** kmer_profile_table_new(GtUword numofsequences,
                                     GtUword r,
                                     GtUword length)
{
  GtKmercount **tau;
  GtUword i, range;
  range = pow(r, length);
  tau = gt_malloc((numofsequences)*sizeof(GtKmercount*));
  *tau = gt_calloc((numofsequences),range*sizeof(GtKmercount));
  for (i = 1; i < numofsequences; i++)
  {
    tau[i] = tau[i-1]+range;
  }
  return tau;
}

void kmer_profile_table_delete(GtKmercount **tau)
{
  if (tau != NULL)
  {
    gt_free(*tau);
    gt_free(tau);
  }
}

void gt_get_codes(const GtEncseq *encseq,
                  GtUword len,
                  GtKmercount **tau)
{
  GtKmercodeiterator *kc_iter;
  const GtKmercode *kmercode;
  GtUword seqnum,
          pos = 0;
  gt_assert(encseq != NULL);
  kc_iter = gt_kmercodeiterator_encseq_new(encseq,
                                           GT_READMODE_FORWARD,
                                           len,
                                           pos);
  while ((kmercode = gt_kmercodeiterator_encseq_next(kc_iter)) != NULL)
  {
    pos = (gt_kmercodeiterator_encseq_get_currentpos(kc_iter)-len);
    if (!(kmercode->definedspecialposition))
    {
      seqnum = gt_encseq_seqnum(encseq, pos);
      tau[seqnum][kmercode->code]++;
    }
  }
  gt_kmercodeiterator_delete(kc_iter);
}

Score *calc_qgram(const GtEncseq *encseq_first,
                  const GtEncseq *encseq_second,
                  GtUword r,
                  GtUword q,
                  GtError *err)
{
  GtUword numofseqfirst,
          numofseqsecond;
  GtKmercount **tu,
              **tv;
  GtUword i, j, l, range;
  Score *score;
  bool compare = true;

  gt_error_check(err);
  gt_assert(encseq_first != NULL && encseq_second != NULL && r > 0 && q > 0);

  if (encseq_first == encseq_second)
  {
    compare = false;
  }
  range = pow(r,q);
  numofseqfirst = gt_encseq_num_of_sequences(encseq_first);
  tu = kmer_profile_table_new(numofseqfirst, r, q);
  gt_get_codes(encseq_first, q, tu);
  if (!compare)
  {
    numofseqsecond = numofseqfirst;
    tv = tu;
  }
  else
  {
    numofseqsecond = gt_encseq_num_of_sequences(encseq_second);
    tv = kmer_profile_table_new(numofseqsecond, r, q);
    gt_get_codes(encseq_second, q, tv);
  }
  score = gt_malloc(sizeof(Score)*(numofseqfirst*numofseqsecond));
  score->pos = 0;
  for (i = 0; i < numofseqfirst; i++)
  {
    GtUword startidx = (encseq_first == encseq_second)? i+1 : 0;
    for (j = startidx; j < numofseqsecond; j++)
    {
      double dist = 0;
      for (l = 0; l < range; l++)
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
  if (!compare)
  {
    kmer_profile_table_delete(tu);
  }
  else
  {
    kmer_profile_table_delete(tv);
    kmer_profile_table_delete(tu);
  }
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
  GtKmercount **tu,
              **tv;
  GtUword i, j, l,
          tmp;
  bool compare = true;
  Score *score;

  if (encseq_first == encseq_second)
  {
    compare = false;
  }
  numofseqfirst = gt_encseq_num_of_sequences(encseq_first);
  tu = kmer_profile_table_new(numofseqfirst, r, k);
  gt_get_codes(encseq_first, k, tu);
  if (!compare)
  {
    numofseqsecond = numofseqfirst;
    tv = tu;
  }
  else
  {
    numofseqsecond = gt_encseq_num_of_sequences(encseq_second);
    tv = kmer_profile_table_new(numofseqsecond, r, k);
    gt_get_codes(encseq_second, k, tv);
  }
  score = gt_malloc(sizeof(Score)*(numofseqfirst*numofseqsecond));
  score->pos = 0;

  for (i = 0; i < numofseqfirst; i++)
  {
    GtUword startidx = (encseq_first == encseq_second)? i+1 : 0;
    for (j = startidx; j < numofseqsecond; j++)
    {
      tmp = 0;
      for (l = 0; l < pow(r,k); l++)
      {
        if (tu[i][l] > 0 && tv[j][l] > 0)
           tmp += MIN(tu[i][l],tv[j][l]);
      }
      GtUword length_u = gt_encseq_seqlength(encseq_first, i);
      GtUword length_v = gt_encseq_seqlength(encseq_second, j);
      double dist = (double)tmp/(double)(MIN(length_u, length_v)-k+1);
      score[score->pos].dist = dist;
      score[score->pos].seqnum_u = i;
      score[score->pos].seqnum_v = j;
      score->pos++;
    }
  }
  if (!compare)
  {
    kmer_profile_table_delete(tu);
  }
  else
  {
    kmer_profile_table_delete(tv);
    kmer_profile_table_delete(tu);
  }
  return score;
}

Score *calc_edist(GtEncseq *encseq_first,
                  GtEncseq *encseq_second,
                  GtStr *scorematrix,
                  int indelscore,
                  GtError *err)
{
  gt_error_check(err);
  gt_assert(encseq_first != NULL);

  Scorematrix *smatrix;
  FILE *fp;
  GtUword numofseqfirst,
          numofseqsecond,
          m,
          n;
  GtUword i, j, u, v,
          tmp,
          maximum,
          startone,
          starttwo;
  GtWord ins, del, rep,
         **table;
  Score *score;
  bool haserr = false;

  fp = fopen(gt_str_get(scorematrix), "r");
  if (fp == NULL)
  {
    gt_error_set(err,"Scorematrixfile %s does not exist.\n",
                 gt_str_get(scorematrix));
    haserr = true;
  }
  if (!haserr)
  {
    smatrix = read_score(fp);
    gt_assert(smatrix != NULL);
    fclose(fp);

    numofseqfirst = gt_encseq_num_of_sequences(encseq_first);
    if (encseq_first && encseq_second)
    {
      numofseqsecond = gt_encseq_num_of_sequences(encseq_second);

      score = malloc(sizeof(Score)*((numofseqfirst*numofseqsecond)+1));
      score->pos = 0;

      for (u = 0; u < numofseqfirst; u++)
      {
        for (v = 0; v < numofseqsecond; v++)
        {
          startone = gt_encseq_seqstartpos(encseq_first,u);
          starttwo = gt_encseq_seqstartpos(encseq_second,v);
          m = gt_encseq_seqlength(encseq_first,u);
          n = gt_encseq_seqlength(encseq_second,v);

          gt_array2dim_malloc(table, m, n);
          table[0][0] = 0;
          for (i = 1; i <= m; i++)
            table[i][0] = table[i-1][0] + indelscore;

          for (j = 1; j <= n; j++)
          {
            table[0][j] = table[0][j-1] + indelscore;
            for (i = 1; i <= m; i++)
            {
              ins = table[i][j-1] + indelscore;
              del = table[i-1][j] + indelscore;
              char one = gt_encseq_get_decoded_char(encseq_first,
                                                    i+startone-1,
                                                    GT_READMODE_FORWARD);
              char two = gt_encseq_get_decoded_char(encseq_second,
                                                    j+starttwo-1,
                                                    GT_READMODE_FORWARD);
              rep = table[i-1][j-1] + access_scorematrix(smatrix, one, two);
              tmp = MAX(ins, del);
              maximum = MAX(tmp, rep);
              if (maximum == LONG_MAX)
                haserr = true;
              table[i][j] = maximum;
            }
          }
          score[score->pos].dist = table[m][n];
          score[score->pos].seqnum_u = u;
          score[score->pos].seqnum_v = v;
          score->pos++;
          gt_array2dim_delete(table);
        }
      }
    }
    else
    {
      score = malloc(sizeof(Score)*((numofseqfirst*numofseqfirst)+1));
      score->pos = 0;

      for (u = 0; u < numofseqfirst; u++)
      {
        for (v = (u+1); v < numofseqfirst; v++)
        {
          startone = gt_encseq_seqstartpos(encseq_first,u);
          starttwo = gt_encseq_seqstartpos(encseq_first,v);
          m = gt_encseq_seqlength(encseq_first,u);
          n = gt_encseq_seqlength(encseq_first,v);

          gt_array2dim_malloc(table, m, n);
          table[0][0] = 0;
          for (i = 1; i <= m; i++)
            table[i][0] = table[i-1][0] + indelscore;

          for (j = 1; j <= n; j++)
          {
            table[0][j] = table[0][j-1] + indelscore;
            for (i = 1; i <= m; i++)
            {
              ins = table[i][j-1] + indelscore;
              del = table[i-1][j] + indelscore;
              char one = gt_encseq_get_decoded_char(encseq_first,
                                                    i+startone-1,
                                                    GT_READMODE_FORWARD);
              char two = gt_encseq_get_decoded_char(encseq_first,
                                                    j+starttwo-1,
                                                    GT_READMODE_FORWARD);
              rep = table[i-1][j-1] + access_scorematrix(smatrix, one, two);
              tmp = MAX(ins, del);
              maximum = MAX(tmp, rep);
              if (maximum == LONG_MAX)
                haserr = true;
              table[i][j] = maximum;
            }
          }
          score[score->pos].dist = table[m][n];
          score[score->pos].seqnum_u = u;
          score[score->pos].seqnum_v = v;
          score->pos++;
          gt_array2dim_delete(table);
        }
      }
    }
  delete_scorematrix(smatrix);
  }
  return (haserr)? NULL : score;
}

void calc_maxmatches(GtStrArray *seq,
                     Suffixarray *suffixarray,
                     GtError *err)
{
  GtSeqIterator *seqit;
  bool haserr = false;
  GtUword totallength,
          maxpreflength = 0,
          queryunitnum;
  int retval;
  const GtUchar *query;
  GtUword querylen;
  char *desc = NULL;
  GtUword score = 0;

  seqit = gt_seq_iterator_sequence_buffer_new(seq, err);
  if (seqit == NULL)
    haserr = true;
  else
  {
    gt_assert(suffixarray);
    gt_seq_iterator_set_symbolmap(seqit,
                                  gt_alphabet_symbolmap(gt_encseq_alphabet(
                                                        suffixarray->encseq)));

    for (queryunitnum = 0; /* Nothing */; queryunitnum++)
    {
      retval = gt_seq_iterator_next(seqit,
                                    &query,
                                    &querylen,
                                    &desc,
                                    err);
      if (retval < 0)
      {
        haserr = true;
        break;
      }
      if (retval == 0)
        break;
      if (!haserr)
      {
        GtUword i = 0;
        score = 0;
        while (i < querylen)
        {
          totallength = gt_encseq_total_length(suffixarray->encseq);
          maxpreflength = gt_findmaximalprefixinESA(suffixarray->encseq,
                                                    suffixarray->readmode,
                                                    totallength,
                                                    suffixarray->suftab,
                                                    &query[i],
                                                    querylen-i);
          /* i == i + maxpreflength <=> maxpreflength == 0 ? */
          if (maxpreflength == 0 || maxpreflength == querylen)
          {
            break;
          }
          score++;
          i = i + maxpreflength + 1; /* + 1? */
        }
        printf("Sequenz: "GT_WU"\tScore: "GT_WU"\n", queryunitnum, score);
      }
    }
    gt_seq_iterator_delete(seqit);
  }
}
