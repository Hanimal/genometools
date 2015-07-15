#include "core/encseq_api.h"
#include "core/str_api.h"
#include "core/codetype.h"
#include "core/error.h"
#include "core/warning_api.h"
#include "core/alphabet_api.h"

#include "match/sfx-mappedstr.h"

#include "stdbool.h"
#include "assert.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

//#include "core/str_api.h"
//#include "core/str_array_api.h"

#include "fscore.h"

typedef struct{
    GtUword *tu, 
			*tv,
			*C;
}Score;

void add_to_C(set_of_qword *sq, int icode)
{
    if(sq->j==sq->size)
    {
        sq->size=sq->size+SIZE;
        sq->C=(int*)realloc(sq->C, sq->size*sizeof(int));
    }
    sq->C[sq->j]=icode;
    sq->j++;
}

GtUword gt_get_kmercodes(const GtEncseq *encseq, 
								unsigned int kmerlen)
{
  GtKmercodeiterator *kc_iter;
  const GtKmercode *kmercode;
  GtUword numberofkmerscollected = 0;

  gt_assert(encseq != NULL);
  kc_iter = gt_kmercodeiterator_encseq_new(encseq, GT_READMODE_FORWARD, kmerlen,
                                           0);
  while ((kmercode = gt_kmercodeiterator_encseq_next(kc_iter)) != NULL) 
  {
    if (!kmercode->definedspecialposition) 
    {
      /* store (kmercode, seqnum, endpos) in array */

      //process kmercode->code;

    }
  }
  gt_kmercodeiterator_delete(kc_iter);
  return numberofkmerscollected;
}

static GtEncseq *gt_encseq_get_encseq(const char *seqfile,
                                      GtError *err)
{
  GtEncseqLoader *encseq_loader;
  GtEncseq *encseq;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(seqfile);
  
  encseq_loader = gt_encseq_loader_new();
  if (!(encseq = gt_encseq_loader_load(encseq_loader, seqfile, err)))
    had_err = -1;

  gt_encseq_loader_delete(encseq_loader);
  if (!had_err) 
  {
    if (!gt_encseq_has_description_support(encseq))
      gt_warning("Missing description support for file %s", seqfile);
    return encseq;
  } 
  else
    return NULL;
}

/*GtUword min(GtUword fst, GtUword snd)
{
	if(fst < snd)
		return (fst);
	else
		return (snd);
}*/

/*void match_alphabetcode(GtAlphabet *alpha, 
							GtUword *alpha_code,
							GtError *err)
{
	gt_error_check(err);
	assert(alpha != NULL);
    GtUword i;
    GtUword r;
    const GtUchar* alphabet;
    
    r = gt_alphabet_size(alpha);
    alphabet = gt_alphabet_characters(alpha);

    for(i=0; i<UCHAR_MAX; i++)
        alpha_code[i] = -1;
    for(i=0; i<r; i++)
    {
        if(alpha_code[(GtUword)alphabet[i]] ==-1)
            alpha_code[(GtUword)alphabet[i]] = i;
        else
        {
            gt_error_set(err,"The same symbol %c occured more than" 
						"once in the given alphabet\n", alphabet[i]);
		}
    }
}*/
/*
 * encseq : contains the sequence u and v
 * r : length of alphabet
 * k : length of kmer
 * alpha_code : integercode for each symbol of the alphabet
 */
/*GtUword f_score(GtEncseq *encseq, 
				GtUword seqnum_u, 
				GtUword seqnum_v,
				GtUword r, 
				GtUword k, 
				GtUword *alpha_code, 
				GtError *err)
{
	gt_error_check(err);
    assert(encseq != NULL && alpha_code != NULL);
    assert(r > 0 && k > 0);
    
    GtUword i, 
			tmp, 
			dist,
			pos_u,
			pos_v,
			int_code = 0, 
			idx_C = 0,
			size = 10;
			
    GtUword *tu, 
			*tv, 
			*factor, 
			*C;

    GtEncseqReader* reader;
    
    GtUchar symbol_old = NULL, 
			symbol_new = NULL;
			
	bool haserr = false;
	
	pos_u = gt_encseq_seqstartpos(encseq, seqnum_u);
	pos_v = gt_encseq_seqstartpos(encseq, seqnum_v);
    
    tu = calloc(pow(r, k), sizeof(GtUword));
    tv = calloc(pow(r, k), sizeof(GtUword));
    C = malloc(sizeof(GtUword)*size);
    reader = gt_encseq_create_reader_with_readmode(encseq, 
                                                   GT_READMODE_FORWARD,
                                                   pos_u);
                                                   
    //define factor for position in sequence
    factor=malloc(sizeof(GtUword)*k);
    factor[0]=1;
    for(i=1; i<k; i++)
        factor[i]=factor[i-1]*r;
        
    //integer_code for first kmer in u
    for(i = 0; i < k; i++)
    {   
		if(gt_encseq_position_is_separator(encseq, pos_u, 
											GT_READMODE_FORWARD))
			break;         
		symbol_old = gt_encseq_reader_next_decoded_char(reader);		
        tmp = alpha_code[(GtUword)symbol_old];
        if(tmp == -1)
            return(INT_MAX);
        int_code += tmp*factor[k-1-i];
        pos_u++;
    }
    tu[int_code]++;

	//integer_code for the rest of u
    while(!(gt_encseq_position_is_separator(encseq, pos_u, 
									GT_READMODE_FORWARD)) && !haserr)
	{
		symbol_new = gt_encseq_reader_next_decoded_char(reader);
		tmp = alpha_code[(GtUword)symbol_new];
		if(tmp == -1)
        {
            gt_error_set(err,"The same symbol "GT_WU" does not occur"
							"in the given alphabet\n", tmp);
			haserr = true;
		}
        int_code = (int_code-alpha_code[(GtUword)symbol_old]*factor[k-1])*r+tmp; 
        tu[int_code]++; 
        symbol_old = symbol_new;
        pos_u++;
	}
    
	//integer_code for first kmer in v
    int_code=0;
    gt_encseq_reader_reinit_with_readmode(reader, encseq,
                                          GT_READMODE_FORWARD, pos_v);
    for(i = 0; i < k; i++)
    {  		
		if(gt_encseq_position_is_separator(encseq, pos_v, 
										GT_READMODE_FORWARD) || haserr)
			break;         
		symbol_old = gt_encseq_reader_next_decoded_char(reader);		
        tmp = alpha_code[(GtUword)symbol_old];
        if(tmp == -1)
        {
            gt_error_set(err,"The same symbol "GT_WU" does not occur"
							"in the given alphabet\n", tmp);
			haserr = true;
		}
        int_code += tmp*factor[k-1-i];
        pos_v++;
    }
    //if kmer occurs in both sequences, then add to C
    tv[int_code]++;
    if(tu[int_code] > 0)
    {
		C[idx_C] = int_code;
		idx_C++;
	}
        
	//integer_code for te rest of v
    while(!(gt_encseq_position_is_separator(encseq, pos_v, 
									GT_READMODE_FORWARD)) && !haserr)
	{
		symbol_new = gt_encseq_reader_next_decoded_char(reader);
		tmp = alpha_code[(GtUword)symbol_new];
		if(tmp == -1)
            return(INT_MAX);
        int_code = (int_code-alpha_code[(GtUword)symbol_old]*factor[k-1])*r+tmp; 
        tu[int_code]++; 
        if(tu[int_code] > 0)
		{
			if(idx_C == size)
			{
				size += 10;
				C= (GtUword*) realloc(C, size*sizeof(GtUword));
			}
			C[idx_C] = int_code;
			idx_C++;
		}
        symbol_old = symbol_new;
        pos_v++;
	}
   
    GtUword length_u = gt_encseq_seqlength(encseq, seqnum_u);
    GtUword length_v = gt_encseq_seqlength(encseq, seqnum_v);
    
    tmp = 0;
    for(i = 0; i < idx_C; i++)
		tmp += min(tu[C[i]],tv[C[i]]);
	dist = tmp/(min(length_u, length_v)-k+1);
	
	free(C);
    free(tu);
    free(tv);
    gt_encseq_reader_delete(reader);
    return((haserr)? UINT_MAX : dist);
}*/
