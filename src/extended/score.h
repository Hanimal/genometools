#ifndef SCORE_H
#define SCORE_H

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
                
Score* calc_qgram(GtEncseq *encseq_first, 
                  GtEncseq *encseq_second,
                  GtUword r, 
                  GtUword q,
                  GtError *err);

Score* calc_edist(GtEncseq *encseq_first, 
                  GtEncseq *encseq_second,
                  GtStr *scorematrix,
                  GtError *err);
#endif
