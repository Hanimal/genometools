#ifndef FSCORE_H
#define FSCORE_H

typedef struct{
  GtUword seqnum_u,
          seqnum_v,
          pos;
  float dist;
}FScore;

FScore* fscore(GtEncseq *encseq_first, 
                GtEncseq *encseq_second, 
                GtUword r, 
                GtUword k, 
                GtError *err);

#endif
