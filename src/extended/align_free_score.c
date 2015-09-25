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
#include "core/encseq_api.h"
#include "core/error.h"
#include "core/alphabet_api.h"
#include "core/unused_api.h"
#include "core/types_api.h"
#include "core/seq_iterator_sequence_buffer_api.h"
#include "core/assert_api.h"
#include "core/minmax.h"
#include "core/array2dim_api.h"
#include "core/score_matrix.h"

#include "align_free_score.h"

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
  if (tau)
  {
    gt_free(*tau);
    gt_free(tau);
  }
}

void get_codes(const GtEncseq *encseq,
               GtUword len,
               GtKmercount **tau)
{
  GtKmercodeiterator *kc_iter;
  const GtKmercode *kmercode;
  GtUword seqnum,
          pos = 0;

  gt_assert(encseq);
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

double **calc_qgram(const GtEncseq *encseq_first,
                     const GtEncseq *encseq_second,
                     GtUword r,
                     GtUword q,
                     GtError *err)
{
  GtUword numofseqfirst,
          numofseqsecond,
          i, j, l, range;
  double **score;
  GtKmercount **tu,
              **tv;
  bool compare = true;

  gt_error_check(err);
  gt_assert(encseq_first && encseq_second && r > 0 && q > 0);

  if (encseq_first == encseq_second)
  {
    compare = false;
  }
  numofseqfirst = gt_encseq_num_of_sequences(encseq_first);
  tu = kmer_profile_table_new(numofseqfirst, r, q);
  get_codes(encseq_first, q, tu);
  if (!compare)
  {
    numofseqsecond = numofseqfirst;
    tv = tu;
  }
  else
  {
    numofseqsecond = gt_encseq_num_of_sequences(encseq_second);
    tv = kmer_profile_table_new(numofseqsecond, r, q);
    get_codes(encseq_second, q, tv);
  }

  gt_array2dim_malloc(score, numofseqfirst, numofseqsecond);
  range = pow(r,q);
  for (i = 0; i < numofseqfirst; i++)
  {
    GtUword startidx = (!compare)? i+1 : 0;
    for (j = startidx; j < numofseqsecond; j++)
    {
      GtUword dist = 0;
      for (l = 0; l < range; l++)
      {
        if (tu[i][l] > 0 || tv[j][l] > 0)
        {
          dist += abs(tu[i][l] - tv[j][l]);
        }
      }
      GtUword length_u = gt_encseq_seqlength(encseq_first, i);
      GtUword length_v = gt_encseq_seqlength(encseq_second, j);
      GtUword qmax = (length_u-q+1)+(length_v-q+1);
      score[i][j] = (double)(qmax-dist)/(double)qmax;
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

double **calc_fscore(const GtEncseq *encseq_first,
                     const GtEncseq *encseq_second,
                     GtUword r,
                     GtUword k,
                     GtError *err)
{
  GtUword numofseqfirst,
          numofseqsecond,
          i, j, l,
          tmp, range;
  GtKmercount **tu,
              **tv;
  bool compare = true;
  double **score;

  gt_error_check(err);
  gt_assert(encseq_first && encseq_second && r > 0 && k > 0);

  if (encseq_first == encseq_second)
  {
    compare = false;
  }
  numofseqfirst = gt_encseq_num_of_sequences(encseq_first);
  tu = kmer_profile_table_new(numofseqfirst, r, k);
  get_codes(encseq_first, k, tu);
  if (!compare)
  {
    numofseqsecond = numofseqfirst;
    tv = tu;
  }
  else
  {
    numofseqsecond = gt_encseq_num_of_sequences(encseq_second);
    tv = kmer_profile_table_new(numofseqsecond, r, k);
    get_codes(encseq_second, k, tv);
  }

  gt_array2dim_malloc(score, numofseqfirst, numofseqsecond);
  range = pow(r,k);
  for (i = 0; i < numofseqfirst; i++)
  {
    GtUword startidx = (!compare)? i+1 : 0;
    for (j = startidx; j < numofseqsecond; j++)
    {
      tmp = 0;
      for (l = 0; l < range; l++)
      {
        if (tu[i][l] > 0 && tv[j][l] > 0)
        {
           tmp += MIN(tu[i][l],tv[j][l]);
        }
      }
      GtUword length_u = gt_encseq_seqlength(encseq_first, i);
      GtUword length_v = gt_encseq_seqlength(encseq_second, j);
      double dist = (double)tmp/(double)(MIN(length_u, length_v)-k+1);
      score[i][j] = dist;
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

GtWord **calc_edist(const GtEncseq *encseq_first,
                    const GtEncseq *encseq_second,
                    GtStr *scorematrix,
                    int indelscore,
                    GtError *err)
{
  GtScoreMatrix *smatrix;
  GtUword numofseqfirst,
          numofseqsecond,
          m, n, i, j, u, v,
          startone,
          starttwo;
  GtWord ins, del, rep, tmp,
         maximum,
         **table,
         **score;
  bool haserr = false,
       compare = true;

  gt_error_check(err);
  gt_assert(encseq_first && encseq_second && scorematrix);

  GtAlphabet *alpha = gt_encseq_alphabet(encseq_first);
  if (!(smatrix = gt_score_matrix_new_read(gt_str_get(scorematrix),
                                           alpha, err)))
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (encseq_first == encseq_second)
    {
      compare = false;
    }

    numofseqfirst = gt_encseq_num_of_sequences(encseq_first);
    if (!compare)
    {
      numofseqsecond = numofseqfirst;
    }
    else
    {
      numofseqsecond = gt_encseq_num_of_sequences(encseq_second);
    }
    gt_array2dim_malloc(score, numofseqfirst, numofseqsecond);

    for (u = 0; u < numofseqfirst; u++)
    {
      GtUword startidx = (!compare)? u+1 : 0;
      for (v = startidx; v < numofseqsecond; v++)
      {
        startone = gt_encseq_seqstartpos(encseq_first,u);
        starttwo = gt_encseq_seqstartpos(encseq_second,v);
        m = gt_encseq_seqlength(encseq_first,u);
        n = gt_encseq_seqlength(encseq_second,v);
        gt_array2dim_malloc(table, m+1, n+1);
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
            GtUchar code1 = gt_alphabet_encode(alpha, one);
            GtUchar code2 = gt_alphabet_encode(alpha, two);
            rep = table[i-1][j-1] + gt_score_matrix_get_score(smatrix,
                                                              code1, code2);
            tmp = MAX(ins, del);
            maximum = MAX(tmp, rep);
            if (maximum == LONG_MAX)
              haserr = true;
            table[i][j] = maximum;
          }
        }
        score[u][v] = table[m][n];
        gt_array2dim_delete(table);
      }
    }
    gt_score_matrix_delete(smatrix);
  }
  return (haserr)? NULL : score;
}

Maxmatch *calc_maxmatches(const GtStrArray *seq,
                          const Suffixarray *suffixarray,
                          GtError *err)
{
  GtSeqIterator *seqit;
  bool haserr = false;
  GtUword totallength,
          maxpreflength = 0,
          queryunitnum,
          match,
          size = 10,
          querylen;
  Maxmatch *score;
  int retval;
  const GtUchar *query;
  char *desc = NULL;

  seqit = gt_seq_iterator_sequence_buffer_new(seq, err);
  if (seqit == NULL)
  {
    haserr = true;
  }
  else
  {
    gt_assert(suffixarray);
    gt_seq_iterator_set_symbolmap(seqit,
                                  gt_alphabet_symbolmap(gt_encseq_alphabet(
                                                        suffixarray->encseq)));
    score = gt_malloc(sizeof(Maxmatch));
    score->dist = gt_malloc(sizeof(double)*size);
    for (queryunitnum = 0; /* Nothing */; queryunitnum++)
    {
      if (queryunitnum == size)
      {
        size += 10;
        score->dist = gt_realloc(score->dist,sizeof(double)*size);
      }
      retval = gt_seq_iterator_next(seqit, &query, &querylen, &desc, err);
      if (retval < 0)
      {
        haserr = true;
        break;
      }
      if (retval == 0)
      {
        break;
      }
      if (!haserr)
      {
        GtUword i = 0;
        match = 0;
        GtUword maxmm = querylen;
        while (i < querylen)
        {
          totallength = gt_encseq_total_length(suffixarray->encseq);
          maxpreflength = gt_findmaximalprefixinESA(suffixarray->encseq,
                                                    suffixarray->readmode,
                                                    totallength,
                                                    suffixarray->suftab,
                                                    &query[i],
                                                    querylen-i);
          if (maxpreflength == querylen)
          {
            break;
          }
          match++;
          i = i + maxpreflength + 1; /* + 1? */
        }
        score->dist[queryunitnum] = (double)(maxmm-match)/(double)maxmm;
      }
    }
    gt_seq_iterator_delete(seqit);
  }
  score->numofseq = queryunitnum;
  return (haserr)? NULL : score;
}
