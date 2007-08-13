/*
** Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
**  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**  
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**  
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
**  
*/
/**
 * \file checkblockcomp.c
 * \brief Experimentally try methods for block-compressed sequences.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif /* HAVE_STDINT_H */
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif /* HAVE_INTTYPES_H */
#include <stdlib.h>
#include <string.h>

#include <time.h>
#include <sys/time.h>

#include <libgtcore/env.h>
#include <libgtcore/ensure.h>
#include <libgtcore/option.h>
#include <libgtcore/str.h>
#include <libgtcore/versionfunc.h>
#include <encidxseq.h>
#include <warts.h>

enum {
  BLOCKSIZE = 8,
};

static OPrval
parseOptions(int *parsed_args, int argc, char **argv,
             Env *env)
{
  OptionParser *op;
  Option *randSeedOption;
  OPrval oprval;
  long seedVal;
  {
    struct timeval seed;
    gettimeofday(&seed, NULL);
    seedVal = seed.tv_sec + seed.tv_usec;
  }

  env_error_check(env);
  op = option_parser_new("indexname",
                         "Map <indexname> and build block composition index.",
                         env);
  randSeedOption = option_new_long("random-seed", "specify start value"
                                   " for random number generator", &seedVal,
                                   seedVal, env);
  option_parser_add_option(op, randSeedOption, env);

  oprval = option_parser_parse_min_max_args(op, parsed_args, argc,
                                            (const char **)argv,
                                            versionfunc, 1, 1, env);
  fprintf(stderr, "seedval = %lu\n", seedVal);
  srandom(seedVal);
  option_parser_delete(op, env);
  return oprval;
}


int
main(int argc, char *argv[])
{
//  MRAEnc *alphabet;
  struct encIdxSeq *seq;
  int dnasymcount = 4;
//  uint8_t charMapping[256];
  enum rangeStoreMode modes[] = { BLOCK_COMPOSITION_INCLUDE,
                                  REGIONS_LIST };
  Str *inputProject;
  int parsedArgs;
  int had_err = 0;
  Env *env = env_new();
  env_error_check(env);
  switch (parseOptions(&parsedArgs, argc, argv, env))
  {
    case OPTIONPARSER_OK:
      break;
    case OPTIONPARSER_ERROR:
      return EXIT_FAILURE;
    case OPTIONPARSER_REQUESTS_EXIT:
      return EXIT_SUCCESS;
  }

//  memset(charMapping, 0, sizeof(charMapping));
//  charMapping[(uint8_t)'A'] = 0;
//  charMapping[(uint8_t)'C'] = 1;
//  charMapping[(uint8_t)'T'] = 2;
//  charMapping[(uint8_t)'G'] = 3;
  inputProject = str_new_cstr(argv[parsedArgs], env);
//  alphabet = MRAEncUInt8New(1, &dnasymcount, charMapping);
  seq = newBlockEncIdxSeq(modes, inputProject,
                          /* int blockSize */BLOCKSIZE, env);
  ensure(had_err, seq);
  if(had_err)
    return EXIT_FAILURE;
  {
    int i;
    Symbol exampleBlock[BLOCKSIZE];
    size_t indices[2];
    for(i = 0; i < BLOCKSIZE; ++i)
    {
      exampleBlock[i] = ((unsigned long)random())%dnasymcount;
    }
    searchBlock2IndexPair(seq, exampleBlock, indices, env);
    exampleBlock[0] = 3;
    exampleBlock[1] = 0;
    exampleBlock[2] = 0;
    exampleBlock[3] = 2;
    exampleBlock[4] = 2;
    exampleBlock[5] = 2;
    exampleBlock[6] = 2;
    exampleBlock[7] = 2;
    searchBlock2IndexPair(seq, exampleBlock, indices, env);
    fprintf(stderr, "indices: %lu, %lu\n", (unsigned long)indices[0],
            (unsigned long)indices[1]);
  }
  {
    int integrity = verifyIntegrity(seq, inputProject, 100000, stderr, env);
    if((integrity))
    {
      if(integrity == -1)
        perror("I/O error when checking index integrity");
      else
        fputs("Integrity check failed for index.\n", stderr);
      return EXIT_FAILURE;
    }
  }
  deleteEncIdxSeq(seq, env);
//  MRAEncDelete(alphabet);
  str_delete(inputProject, env);
  env_delete(env);
  return EXIT_SUCCESS;
}
