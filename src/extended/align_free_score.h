#ifndef ALIGN_FREE_SCORE_H
#define ALIGN_FREE_SCORE_H

typedef struct{
  GtUword seqnum_u,
          seqnum_v,
          pos;
  float dist;
}Score;

Score* calc_fscore(GtEncseq *encseq_first,
                   GtEncseq *encseq_second,
                   GtUword r,
                   GtUword k,
                   GtError *err);

Score *calc_qgram(const GtEncseq *encseq_first,
                  const GtEncseq *encseq_second,
                  GtUword r,
                  GtUword q,
                  GtError *err);

Score *calc_edist(GtEncseq *encseq_first,
                  GtEncseq *encseq_second,
                  GtStr *scorematrix,
                  int indelscore,
                  GT_UNUSED GtError *err);

void calc_maxmatches(GtStrArray *seq,
                     Suffixarray *suffixarray,
                     GtError *err);
#endif
