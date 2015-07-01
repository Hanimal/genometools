#include "core/encseq_api.h"
#include "core/str_api.h"
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "assert.h"
#include "stdbool.h"
#include "f_score.h"

GtUword min(GtUword fst, GtUword scd)
{
	if(fst < scd)
		return (fst);
	else
		return (scd);
}

GtUword *match_alphabetcode(GtAlphabet *alpha, GtError *err)
{
	gt_error_check(err);
	assert(alpha != NULL);
    GtUword i;
    GtUword *alpha_tab;
    GtUword r;
    const GtUchar* alphabet;
    bool haserr = false;
    
    alpha_tab = malloc(sizeof(GtUword)*UCHAR_MAX);
    r = gt_alphabet_size(alpha);
    alphabet = gt_alphabet_characters(alpha);

    for(i=0; i<UCHAR_MAX; i++)
        alpha_tab[i] = -1;
    for(i=0; i<r; i++)
    {
        if(alpha_tab[(GtUword)alphabet[i]] ==-1)
            alpha_tab[(GtUword)alphabet[i]] = i;
        else
        {
            gt_error_set(err,"The same symbol %c occured more than" 
						"once in the given alphabet\n", alphabet[i]);
			haserr = true;
		}
    }
    return((haserr)? NULL : alpha_tab);
}
/*
 * encseq : contains the sequence u and v
 * r : length of alphabet
 * k : length of kmer
 * alpha_code : integercode for each symbol of the alphabet
 */
GtUword f_score(GtEncseq *encseq, GtUword r, GtUword k, GtUword *alpha_code, GtError *err)
{
	gt_error_check(err);
    assert(encseq != NULL && alpha_code != NULL);
    assert(r > 0 && k > 0);
    
    GtUword read_at = 0, 
			i = 0, 
			tmp = 0, 
			int_code=0, 
			size = 10, 
			idx_C = 0,
			dist = 0;
			
    GtUword *tu, 
			*tv, 
			*factor, 
			*C;

    GtEncseqReader* reader;
    
    GtUchar symbol_old, 
			symbol_new;
    
    tu = calloc(pow(r, k), sizeof(GtUword));
    tv = calloc(pow(r, k), sizeof(GtUword));
    C = malloc(sizeof(GtUword)*size);
    reader = gt_encseq_create_reader_with_readmode(encseq, 
                                                   GT_READMODE_FORWARD,
                                                   0);
                                                   
    /*define factor for position in sequence*/
    factor=malloc(sizeof(GtUword)*k);
    factor[0]=1;
    for(i=1; i<k; i++)
        factor[i]=factor[i-1]*r;
        
    /*integer_code for first kmer in u*/
    for(i = 0; i < k; i++)
    {   
		if(gt_encseq_position_is_separator(encseq, read_at, 
											GT_READMODE_FORWARD))
			break;         
		symbol_old = gt_encseq_reader_next_decoded_char(reader);		
        tmp = alpha_code[(GtUword)symbol_old];
        if(tmp == -1)
            return(INT_MAX);
        int_code += tmp*factor[k-1-i];
        read_at++;
    }
    tu[int_code]++;

	/*integer_code for the rest of u*/
    while(!(gt_encseq_position_is_separator(encseq, read_at, 
											GT_READMODE_FORWARD)))
	{
		GtUchar symbol_new = gt_encseq_reader_next_decoded_char(reader);
		tmp = alpha_code[(GtUword)symbol_new];
		if(tmp == -1)
            return(INT_MAX);
        int_code = (int_code-alpha_code[(GtUword)symbol_old]*factor[k-1])*r+tmp; 
        tu[int_code]++; 
        symbol_old = symbol_new;
        read_at++;
	}
    
	/*integer_code for first kmer in v*/
    int_code=0;
    read_at++;
    for(i = 0; i < k; i++)
    {  		
		if(gt_encseq_position_is_separator(encseq, read_at, 
											GT_READMODE_FORWARD))
			break;         
		symbol_old = gt_encseq_reader_next_decoded_char(reader);		
        tmp = alpha_code[(GtUword)symbol_old];
        if(tmp == -1)
            return(INT_MAX);
        int_code += tmp*factor[k-1-i];
        read_at++;
    }
    /*if kmer occurs in both sequences, then add to C*/
    tv[int_code]++;
    if(tu[int_code] > 0)
    {
		C[idx_C] = int_code;
		idx_C++;
	}
        
	/*integer_code for te rest of v*/
    while(!(gt_encseq_position_is_separator(encseq, read_at, 
											GT_READMODE_FORWARD)))
	{
		GtUchar symbol_new = gt_encseq_reader_next_decoded_char(reader);
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
        read_at++;
	}

    
    GtUword length_u = gt_encseq_seqlength(encseq, 0);
    GtUword length_v = gt_encseq_seqlength(encseq, 1);
    
    tmp = 0;
    for(i = 0; i < idx_C; i++)
		tmp += min(tu[C[i]],tv[C[i]]);
	dist = tmp/(min(length_u, length_v)-k+1);
	
	free(C);
    free(tu);
    free(tv);
    gt_encseq_reader_delete(reader);
    return(dist);
}
