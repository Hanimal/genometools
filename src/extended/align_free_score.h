#ifndef ALIGN_FREE_SCORE_H
#define ALIGN_FREE_SCORE_H

typedef struct{
  double *dist;
  GtUword numofseq;
}Maxmatch;

typedef struct{
    double value;
    bool minRep;
    bool minIns;
    bool minDel;
}Entry;

typedef struct{
    Entry del;
    Entry ins;
    Entry rep;
}AffineNode;

double **calc_fscore(const GtEncseq *encseq_first,
                     const GtEncseq *encseq_second,
                     GtUword r,
                     GtUword k,
                     GtError *err);

double **calc_qgram(const GtEncseq *encseq_first,
                     const GtEncseq *encseq_second,
                     GtUword r,
                     GtUword q,
                     bool distance,
                     GtError *err);

GtWord **calc_edist(const GtEncseq *encseq_first,
                    const GtEncseq *encseq_second,
                    GtStr *scorematrix,
                    int indelscore,
                    GtError *err);

Maxmatch *calc_maxmatches(const GtStrArray *seq,
                          const Suffixarray *suffixarray,
                          bool distance,
                          GtError *err);
                        
GtWord **calc_edist_affine(const GtEncseq *encseq_first,
                           const GtEncseq *encseq_second,
                           GtStr *scorematrix,
                           int gapopen,
                           int gapextend,
                           GtError *err);
#endif
