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

#include "gt_sequencescorer.h"
#include "fscore.h"


typedef struct 
{
  GtEncseqOptions *eopts;
  unsigned int k;
  bool fscore;
  GtStr *indexname;
  GtStrArray *queryfiles;
} GtSequencescorerArguments;

static void* gt_sequencescorer_arguments_new(void)
{
  GtSequencescorerArguments *arguments;
  
  arguments = gt_calloc(1, sizeof(*arguments));
  arguments->indexname = gt_str_new();
  arguments->queryfiles = gt_str_array_new();
  return arguments;
}

static void gt_sequencescorer_arguments_delete(void *tool_arguments)
{
  GtSequencescorerArguments *arguments = tool_arguments;
  
  if (arguments != NULL) 
  {
    gt_str_delete(arguments->indexname);
    gt_str_array_delete(arguments->queryfiles);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_sequencescorer_option_parser_new(void *tool_arguments)
{
  GtSequencescorerArguments *arguments = tool_arguments;
  
  GtOptionParser *op;
  GtOption *option, *fscore, *queryoption;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option] sequencefile [sequencefile ...] "
                            "-ii indexname [...]",
                            "Computes scores.");
  gt_option_parser_set_mail_address(op,"<hannah@rauterberg.eu>");

  option = gt_option_new_uint_min("k","length of kmer",
                                  &arguments->k,
                                  6,
                                  1);
  gt_option_parser_add_option(op, option);
                                  
  fscore = gt_option_new_bool("fscore", "computes fscore",
                              &arguments->fscore, false);
  gt_option_parser_add_option(op, fscore);
  
  
  queryoption = gt_option_new_filename_array("s",
                                             "Specify query files",
                                             arguments->queryfiles);
  gt_option_parser_add_option(op, queryoption);

  arguments->eopts = gt_encseq_options_register_encoding(op,
                                                         arguments->indexname,
                                                         NULL);
  gt_option_imply(option, fscore);
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

static int gt_sequencescorer_runner(GT_UNUSED int argc, GT_UNUSED const char **argv, GT_UNUSED int parsed_args,
                              void *tool_arguments, GT_UNUSED GtError *err)
{
  GtSequencescorerArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  if(arguments->fscore == true)
  {
    gt_error_check(err);
    GtEncseq *encseq_first;
    GtEncseq *encseq_second;
    FScore *score;
    GtAlphabet* alpha;
    GtUword r,
            i;

  
    encseq_first = gt_encseq_get_encseq(gt_str_array_get(arguments->queryfiles,0), err);
    gt_error_check(err);
    encseq_second = gt_encseq_get_encseq(gt_str_array_get(arguments->queryfiles,1), err);
    gt_error_check(err);
    
    //Schoener waere es, noch die alphabets miteinander zu vergleichen
    //bool gt_alphabet_equals(const GtAlphabet *a, const GtAlphabet *b);
    alpha = gt_encseq_alphabet(encseq_first);
    r = gt_alphabet_size(alpha);
    
    score = fscore(encseq_first, 
                   encseq_second, 
                   r, 
                   arguments->k, 
                   err);
    
    for( i = 0; i < score->pos; i++)
    {
      printf("Fscore between sequence "GT_WU" and "GT_WU" is %.3f.\n", 
          score[i].seqnum_u, score[i].seqnum_v, score[i].dist);
    }
    
    gt_encseq_delete(encseq_first);
    gt_encseq_delete(encseq_second);
    free(score);
    return(EXIT_SUCCESS);
  }

  return had_err;
}

GtTool* gt_sequencescorer(void)
{
  return gt_tool_new(gt_sequencescorer_arguments_new,
                     gt_sequencescorer_arguments_delete,
                     gt_sequencescorer_option_parser_new,
                     gt_sequencescorer_arguments_check,
                     gt_sequencescorer_runner);
}
