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
#include "core/assert_api.h"
#include "match/esa-splititv.h"
#include "match/sfx-mappedstr.h"
#include "match/esa-map.h"

#include "gt_sequencescorer.h"
#include "score.h"

typedef struct
{
  unsigned int k;
  unsigned int q;
  int indelscore;
  bool edist;
  bool fscore;
  bool qgram;
  bool maxmatches;
  GtStr *scorematrix;
  GtStrArray *queryfiles;
  GtStrArray *seq;
} GtSequencescorerArguments;

static void* gt_sequencescorer_arguments_new(void)
{
  GtSequencescorerArguments *arguments;

  arguments = gt_calloc(1, sizeof (*arguments));

  arguments->queryfiles = gt_str_array_new();
  arguments->seq = gt_str_array_new();
  arguments->scorematrix = gt_str_new();
  return arguments;
}

static void gt_sequencescorer_arguments_delete(void *tool_arguments)
{
  GtSequencescorerArguments *arguments = tool_arguments;

  if (arguments != NULL)
  {
    gt_str_array_delete(arguments->queryfiles);
    gt_str_array_delete(arguments->seq);
    gt_str_delete(arguments->scorematrix);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_sequencescorer_option_parser_new(void *tool_arguments)
{
  GtSequencescorerArguments *arguments = tool_arguments;

  GtOptionParser *op;
  GtOption *k, *q, *fscore, *queryoption, *qgram, *edist, *scorematrix,
           *indelscore, *maxmatches, *seq;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option] -ii queryfile [queryfile] [...]",
                            "Computes scores.");

  gt_option_parser_set_mail_address(op,"<hannah@rauterberg.eu>");
  
  queryoption = gt_option_new_filename_array("ii", "Specify query files.\n"
                                            "Either Encseq or Suffixarray",
                                             arguments->queryfiles);
  gt_option_parser_add_option(op, queryoption);
  gt_option_is_mandatory(queryoption);
  
  /*fscore*/
  fscore = gt_option_new_bool("fscore", "computes fscore",
                              &arguments->fscore, false);
  gt_option_parser_add_option(op, fscore);
  k = gt_option_new_uint_min("k","length of kmer", &arguments->k, 6, 1);
  gt_option_parser_add_option(op, k);

  /*edist*/
  edist = gt_option_new_bool("edist", "computes editdistance",
                              &arguments->edist, false);
  gt_option_parser_add_option(op, edist);
  indelscore = gt_option_new_int_min("indelscore",
                                     "set score for inserstion or deletion",
                                     &arguments->indelscore, -1, INT_MIN);
  gt_option_parser_add_option(op, indelscore);
  scorematrix = gt_option_new_filename("smatrix", "Specify Scorematrix",
                                       arguments->scorematrix);
  gt_option_parser_add_option(op, scorematrix);
  
  /*qgram*/
  qgram = gt_option_new_bool("qgram", "computes qgram distance",
                              &arguments->qgram, false);
  gt_option_parser_add_option(op, qgram);
  q = gt_option_new_uint_min("q","length of qgram", &arguments->q, 6, 1);
  gt_option_parser_add_option(op, q);

  /*maxmatches*/
  maxmatches = gt_option_new_bool("maxmatches", "computes maximalmatches from\n"
                                  " Sequencefiles with respect to given "
                                  "Suffixtab",
                                  &arguments->maxmatches, false);
  gt_option_parser_add_option(op, maxmatches);
  seq = gt_option_new_filename_array("seq", "Specify Sequencefiles",
                                    arguments->seq);
  gt_option_parser_add_option(op, seq);

  gt_option_imply(fscore, k);
  gt_option_imply(edist, scorematrix);
  gt_option_imply(edist, indelscore);
  gt_option_imply(qgram, q);
  gt_option_imply(maxmatches, seq);
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
    bool haserr = false;
    gt_error_check(err);
    gt_assert(seqfile);

    encseq_loader = gt_encseq_loader_new();
    if (!(encseq = gt_encseq_loader_load(encseq_loader, seqfile, err)))
      haserr = true;
    gt_encseq_loader_delete(encseq_loader);
    if (!haserr)
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
  if(arguments->maxmatches == false)
  {
    if ((gt_str_array_size(arguments->queryfiles) == 0)||
        (gt_str_array_size(arguments->queryfiles) > 2))
    {
      gt_error_set(err,"At least one and at most two sequencefiles allowed.\n");
      haserr = true;
    }
    else if (gt_str_array_size(arguments->queryfiles) == 1)
    {
      encseq_first = gt_encseq_get_encseq(gt_str_array_get(arguments->queryfiles,
                                                           0), err);
      if (!encseq_first)
      {
        gt_error_set(err,"Sequencefile does not exist.\n");
        haserr = true;
      }
    }
    else
    {
      encseq_first = gt_encseq_get_encseq(gt_str_array_get(arguments->queryfiles,
                                                           0), err);
      if (encseq_first == NULL)
      {
        gt_error_set(err,"Sequencefile %s does not exist.\n",
                          gt_str_array_get(arguments->queryfiles,0));
        haserr = true;
      }
      encseq_second = gt_encseq_get_encseq(gt_str_array_get(arguments->queryfiles,
                                                            1), err);
      if (encseq_second == NULL)
      {
        gt_error_set(err,"Sequencefile %s does not exist.\n",
                          gt_str_array_get(arguments->queryfiles,1));
        haserr = true;
      }
      if(!haserr)
      {
        if(!gt_alphabet_equals(gt_encseq_alphabet(encseq_first),
                               gt_encseq_alphabet(encseq_second)))
        {
          gt_error_set(err,"Files encoded with different alphabets.");
          haserr = true;
        }
      }
    }
  }
  if (arguments->fscore == true && !haserr)
  {
    Score *score;
    GtUword r,
            i;
    gt_assert(encseq_first);
    r = gt_alphabet_size(gt_encseq_alphabet(encseq_first));

    score = calc_fscore(encseq_first,
                        encseq_second,
                        r,
                        arguments->k,
                        err);

    for (i = 0; i < score->pos; i++)
    {
      printf("Fscore between sequence "GT_WU" and "GT_WU" is %.3f.\n",
          score[i].seqnum_u, score[i].seqnum_v, score[i].dist);
    }
    free(score);
  }
  if (arguments->qgram == true && !haserr)
  {
    Score *score;
    GtUword r,
            i;
    gt_assert(encseq_first);
    r = gt_alphabet_size(gt_encseq_alphabet(encseq_first));

    score = calc_qgram(encseq_first,
                       encseq_second,
                       r,
                       arguments->q,
                       err);

    for (i = 0; i < score->pos; i++)
    {
      printf("Qgramdistance between sequence "GT_WU" and "GT_WU" is %.0f.\n",
          score[i].seqnum_u, score[i].seqnum_v, score[i].dist);
    }
    free(score);
  }
  if (arguments->edist == true && !haserr)
  {
    GtUword i;
    Score *score;
    gt_assert(encseq_first);
    score = calc_edist(encseq_first,
                       encseq_second,
                       arguments->scorematrix,
                       arguments->indelscore,
                       err);
    if(score == NULL)
      gt_error_set(err,"Error computing Editdistance\n");
    else
    {
      for (i = 0; i < score->pos; i++)
      {
        printf("Editdistance between sequence "GT_WU" and "GT_WU" is %.0f.\n",
            score[i].seqnum_u, score[i].seqnum_v, score[i].dist);
      }
      free(score);
    }
  }
  if (arguments->maxmatches == true && !haserr)
  {
    Suffixarray suffixarray;
    GtLogger *logger;
    logger = gt_logger_new(true, "# ", stderr); 
    gt_mapsuffixarray(&suffixarray,
                      SARR_SUFTAB | SARR_ESQTAB, 
                      gt_str_array_get(arguments->queryfiles,0),
                      logger,
                      err);
    gt_error_check(err);
    calc_maxmatches(arguments->seq,
                    &suffixarray,
                    err);
    gt_freesuffixarray(&suffixarray);
    gt_logger_delete(logger);
  }
  if(encseq_first)
    gt_encseq_delete(encseq_first);
  if(encseq_second)
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
