#ifndef FSCORE_H
#define FSCORE_H

typedef struct{
  GtUword seqnum_u,
          seqnum_v,
          dist;
}Score;

Score* f_score(GtEncseq *encseq_first, 
                GtEncseq *encseq_second, 
                GtUword r, 
                GtUword k, 
                GtError *err);

#endif
