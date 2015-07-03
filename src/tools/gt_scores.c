#include "core/encseq_api.h"
#include "core/str_api.h"
#include "assert.h"
#include "core/error.h"
#include "stdbool.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#include "extended/fscore.h"
#include "gt_scores.h"

typedef struct{
    GtUword score,
			seqnum_u,
			seqnum_v;
}Fscore;

int calc_fscore(int argc, char *argv[]) 
{
	if(argc != 4)
	{
		printf("USAGE: %s <alphabetfile> <sequencefile> <k>", argv[0]);
		return(EXIT_FAILURE);
	}
	else
	{
		GtEncseq *encseq;
		GtAlphabet *alpha; 
		Fscore *fscore;
		GtError* err;
		GtUword *alpha_code;
		GtUword k, 
				i,
				j,
				file_num,
				r,
				idx,
				struct_size = 0,
				size = 10;
		
		sscanf(argv[3], GT_WU, &k);
		alpha_code = malloc(sizeof(GtUword)*UCHAR_MAX);
		fscore = malloc(sizeof(Fscore)*size);
		
		err = gt_error_new();
		gt_error_set_progname(err, argv[0]);
		alpha = gt_alphabet_new_from_file_no_suffix(argv[1], err);
		gt_error_check(err);
		r = gt_alphabet_size(alpha);
		match_alphabetcode(alpha, alpha_code, err);
		gt_error_check(err);
		
		/*encseq = blabla*/
		file_num = gt_encseq_num_of_files(encseq);
		idx = 0;
		for(i = 0; i < file_num; i++)
		{
			for(j = i+1; j < file_num; j++)
			{
				GtUword score = f_score(encseq, 
										i, 
										j,
										r, 
										k, 
										alpha_code, 
										err);
				if(struct_size == size)
				{
					size += 10;
					fscore = realloc(fscore, size*sizeof(Fscore));
				}
				fscore[idx].score = score;
				fscore[idx].seqnum_u = i;
				fscore[idx].seqnum_v = j;
				idx++;
				struct_size++;
			}
		}
		gt_encseq_delete(encseq);
		return(EXIT_SUCCESS);
	}	
}
