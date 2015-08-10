/*
  Copyright (c) 2015 Hannah <hannah@rauterberg.eu>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include "core/ma.h"
#include "core/unused_api.h"
#include "core/alphabet.h"
#include "core/encseq.h"
#include "core/encseq_options.h"
#include "core/fileutils.h"
#include "core/str_array_api.h"
#include "tools/gt_encseq_encode.h"
#include "core/warning_api.h"
#include "assert.h"

#include "gt_sequencescorer.h"
#include "score.h"


typedef struct 
{
  unsigned int k;
  unsigned int q;
  bool fscore;
  bool qgram;
  GtStrArray *queryfiles;
} GtSequencescorerArguments;

static void* gt_sequencescorer_arguments_new(void)
{
  GtSequencescorerArguments *arguments;
  
  arguments = gt_calloc(1, sizeof(*arguments));
  arguments->queryfiles = gt_str_array_new();
  return arguments;
}

static void gt_sequencescorer_arguments_delete(void *tool_arguments)
{
  GtSequencescorerArguments *arguments = tool_arguments;
  
  if (arguments != NULL) 
  {
    gt_str_array_delete(arguments->queryfiles);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_sequencescorer_option_parser_new(void *tool_arguments)
{
  GtSequencescorerArguments *arguments = tool_arguments;
  
  GtOptionParser *op;
  GtOption *k, *q, *fscore, *queryoption, *qgram;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option] -s sequencefile [sequencefile] [...]",
                            "Computes scores.");
                            
  gt_option_parser_set_mail_address(op,"<2rauterb@informatik.uni-hamburg.de>");

  k = gt_option_new_uint_min("k","length of kmer", &arguments->k, 6, 1);
  gt_option_parser_add_option(op, k);
  
  q = gt_option_new_uint_min("q","length of qgram", &arguments->q, 6, 1);
  gt_option_parser_add_option(op, q);
                                  
  fscore = gt_option_new_bool("fscore", "computes fscore", 
                              &arguments->fscore, false);
  gt_option_parser_add_option(op, fscore);
  
  qgram = gt_option_new_bool("qgram", "computes qgram distance", 
                              &arguments->qgram, false);
  gt_option_parser_add_option(op, qgram);
   
  queryoption = gt_option_new_filename_array("s", "Specify query files",
                                             arguments->queryfiles);
  gt_option_parser_add_option(op, queryoption);
  gt_option_is_mandatory(queryoption);
 
  gt_option_imply(k, fscore);
  gt_option_imply(q, qgram);
  return op;
}

static int gt_sequencescorer_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtSequencescorerArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  return had_err;
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
  
    if(!had_err) 
    {
      if (!gt_encseq_has_description_support(encseq))
        gt_warning("Missing description support for file %s", seqfile);
        return encseq;
    } 
    else
      return NULL;
}

static int gt_sequencescorer_runner(GT_UNUSED int argc, 
                                    GT_UNUSED const char **argv, 
                                    GT_UNUSED int parsed_args,
                                    void *tool_arguments, 
                                    GT_UNUSED GtError *err)
{
  GtSequencescorerArguments *arguments = tool_arguments;
  int haserr = 0;

  gt_error_check(err);
  gt_assert(arguments);

  GtEncseq *encseq_first = NULL;
  GtEncseq *encseq_second = NULL;
  GtAlphabet* alpha = NULL;
  
  if((gt_str_array_size(arguments->queryfiles) == 0)|| 
      (gt_str_array_size(arguments->queryfiles) > 2))
  {
    gt_error_set(err,"At least one and at most two sequencefiles allowed.\n");
    haserr = true;
  }
  else if(gt_str_array_size(arguments->queryfiles) == 1)
  {
    encseq_first = gt_encseq_get_encseq(gt_str_array_get(arguments->queryfiles,0), err);
    alpha = gt_encseq_alphabet(encseq_first);
    gt_error_check(err);
  }
  else
  {
    encseq_first = gt_encseq_get_encseq(gt_str_array_get(arguments->queryfiles,0), err);
    gt_error_check(err);
    alpha = gt_encseq_alphabet(encseq_first);
    encseq_second = gt_encseq_get_encseq(gt_str_array_get(arguments->queryfiles,1), err);
    gt_error_check(err);
    GtAlphabet* tmp = gt_encseq_alphabet(encseq_second);
    assert(gt_alphabet_equals(alpha, tmp));
  }
    
  if(arguments->fscore == true)
  {
    Score *score;
    GtUword r,
            i;
            
    if(!haserr)
    {
      r = gt_alphabet_size(alpha);
      
      score = calc_fscore(encseq_first, 
                          encseq_second, 
                          r, 
                          arguments->k, 
                          err);
      
      for( i = 0; i < score->pos; i++)
      {
        printf("Fscore between sequence "GT_WU" and "GT_WU" is %.3f.\n", 
            score[i].seqnum_u, score[i].seqnum_v, score[i].dist);
      }
      free(score);
    }
  }
  if(arguments->qgram == true)
  {
    Score *score;
    GtUword r,
            i;
            
    if(!haserr)
    {
      r = gt_alphabet_size(alpha);
      
      score = calc_qgram(encseq_first, 
                         encseq_second, 
                         r, 
                         arguments->q, 
                         err);
      
      for( i = 0; i < score->pos; i++)
      {
        printf("Qgramdistance between sequence "GT_WU" and "GT_WU" is %.0f.\n", 
            score[i].seqnum_u, score[i].seqnum_v, score[i].dist);
      }
      free(score);
    }
  }
  gt_encseq_delete(encseq_first);
  gt_encseq_delete(encseq_second);
  return haserr;
}

GtTool* gt_sequencescorer(void)
{
  return gt_tool_new(gt_sequencescorer_arguments_new,
                     gt_sequencescorer_arguments_delete,
                     gt_sequencescorer_option_parser_new,
                     gt_sequencescorer_arguments_check,
                     gt_sequencescorer_runner);
}
