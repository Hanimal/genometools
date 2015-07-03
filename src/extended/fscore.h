#ifndef FSCORE_H
#define FSCORE_H

GtUword f_score(GtEncseq *encseq, 
				GtUword seqnum_u, 
				GtUword seqnum_v,
				GtUword r, 
				GtUword k, 
				GtUword *alpha_code, 
				GtError *err);
				
void match_alphabetcode(GtAlphabet *alpha, 
						GtUword *alpha_code,
						GtError *err);

#endif
