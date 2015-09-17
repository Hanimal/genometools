#ifndef ALIGN_FREE_SCORE_H
#define ALIGN_FREE_SCORE_H

typedef struct{
  GtUword *dist,
          numofseq;
}Maxmatch;

double **calc_fscore(const GtEncseq *encseq_first,
                     const GtEncseq *encseq_second,
                     GtUword r,
                     GtUword k,
                     GtError *err);

GtUword **calc_qgram(const GtEncseq *encseq_first,
                     const GtEncseq *encseq_second,
                     GtUword r,
                     GtUword q,
                     GtError *err);

GtWord **calc_edist(const GtEncseq *encseq_first,
                    const GtEncseq *encseq_second,
                    GtStr *scorematrix,
                    int indelscore,
                    GtError *err);

Maxmatch *calc_maxmatches(const GtStrArray *seq,
                          const Suffixarray *suffixarray,
                          GtError *err);
#endif
