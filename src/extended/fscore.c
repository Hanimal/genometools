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
#define SIZE 5

typedef struct{
    GtUword *C,
            pos,
            size;
}Kmers;

void add_to_kmers_table(Kmers *kmers, GtUword kmercode)
{
    if(kmers->pos == kmers->size)
    {
        kmers->size = kmers->size+SIZE;
        kmers->C = (GtUword*)realloc(kmers->C, kmers->size*sizeof(GtUword));
    }
    kmers->C[kmers->pos] = kmercode;
    kmers->pos++;
}

Kmers* new_kmers_table()
{
  Kmers *kmers;
  kmers = malloc(sizeof(Kmers));
  kmers->C = malloc(SIZE*sizeof(GtUword));
  kmers->size = SIZE;
  kmers->pos = 0;
  return kmers;
}

void delete_kmers_table(Kmers *kmers)
{
  free(kmers->C);
  free(kmers);
}

void gt_get_kmercodes(const GtEncseq *encseq, 
                         GtUword kmerlen,
                         GtUword **tau)
{
  GtKmercodeiterator *kc_iter;
  const GtKmercode *kmercode;
  GtUword seqnum;
  gt_assert(encseq != NULL);
  kc_iter = gt_kmercodeiterator_encseq_new(encseq, 
                                           GT_READMODE_FORWARD, 
                                           kmerlen,
                                           0);
                                           
  while ((kmercode = gt_kmercodeiterator_encseq_next(kc_iter)) != NULL) 
  {
    
    seqnum = gt_encseq_seqnum(encseq,
                            gt_kmercodeiterator_encseq_get_currentpos(kc_iter));
    if (!(kmercode->definedspecialposition)) 
    {
      tau[seqnum][kmercode->code]++;
      /* store (kmercode, seqnum, endpos) in array */
    }
  }
  gt_kmercodeiterator_delete(kc_iter);
}

GtUword** new_tau(GtUword numofsequences, GtUword r, GtUword k)
{
  GtUword **tau;
  tau = malloc((numofsequences)*sizeof(GtUword*));
  *tau = calloc((numofsequences),(pow(r, k))*sizeof(GtUword));
  return tau;
}

GtUword min(GtUword fst, GtUword snd)
{
  if(fst < snd)
    return (fst);
  else
    return (snd);
}

Score *f_score(GtEncseq *encseq_first, 
                GtEncseq *encseq_second, 
                GtUword r, 
                GtUword k, 
                GtError *err)
{
  gt_error_check(err);
  assert(encseq_first != NULL && encseq_second != NULL);
  assert(r > 0 && k > 0);
    
  GtUword numofsequences_first;
  GtUword numofsequences_second;
  GtUword **tu,
          **tv,
          i,
          j,
          l,
          p,
          dist = 0,
          tmp = 0,
          score_pos = 0;
  Kmers *kmers;
  Score *score;
  
  
  
  numofsequences_first = gt_encseq_num_of_sequences(encseq_first);
  tu = new_tau(numofsequences_first, r, k);
  numofsequences_second = gt_encseq_num_of_sequences(encseq_second);
  tv = new_tau(numofsequences_second, r, k);
  score = malloc(sizeof(Score)*(numofsequences_first*numofsequences_second));
  
  gt_get_kmercodes(encseq_first, 
                   k,
                   tu);
  gt_get_kmercodes(encseq_second, 
                   k,
                   tv);
  for(i = 0; i < numofsequences_first; i++)
  {
    for(j = 0; j < numofsequences_second; j++)
    {
      kmers = new_kmers_table();
      for(l = 0; l < pow(r,k); l++)
      {
        if(tu[i][l] == tv[j][l])
          add_to_kmers_table(kmers, l);
      }
      GtUword length_u = gt_encseq_seqlength(encseq_first, i);
      GtUword length_v = gt_encseq_seqlength(encseq_second, j);
      for(p = 0; p < kmers->pos; p++)
        tmp += min(tu[i][kmers->C[p]],tv[j][kmers->C[p]]);
      dist = tmp/(min(length_u, length_v)-k+1);
      score[score_pos].dist = dist;
      score[score_pos].seqnum_u = i;
      score[score_pos].seqnum_v = j;
      score_pos++;
      delete_kmers_table(kmers);
    }
  }
  return score;
}

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
/*GtUword f_score2(GtEncseq *encseq, 
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
