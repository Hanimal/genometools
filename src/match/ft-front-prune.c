#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include "core/assert_api.h"
#include "core/unused_api.h"
#include "core/divmodmul.h"
#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/intbits.h"
#include "core/encseq.h"
#include "match/extend-offset.h"
#include "core/ma_api.h"
#include "core/types_api.h"
#include "ft-front-prune.h"
#include "ft-trimstat.h"
#include "core/minmax.h"
#include "ft-polish.h"
#include "ft-front-generation.h"

#define UPDATE_MATCH_HISTORY(MC,MH)\
        if ((MH) & mask)\
        {\
          gt_assert((MC) > 0);\
          (MC)--;\
        }\
        (MH) <<= 1

typedef unsigned int Rowvaluetype;
typedef uint8_t Matchcounttype;
typedef uint8_t Backreferencetype;

typedef struct
{
  uint64_t matchhistory;
  Rowvaluetype row,
               localmatch_count;
  Matchcounttype matchhistory_count;
  Backreferencetype backreference;
} Frontvalue;

#ifndef OUTSIDE_OF_GT
typedef struct
{
  const GtTwobitencoding *twobitencoding;
  const GtEncseq *encseq;
  const GtUchar *bytesequenceptr;
  GtEncseqReader *encseqreader;
  GtUchar *cache_ptr;
  GtAllocatedMemory *sequence_cache;
  GtUword substringlength,
          totallength,
          min_access_pos, /* no position accessed will be smaller than this */
          cache_num_positions; /* number of positions in cache */
  GtUword offset,
          seqstartpos;
  bool read_seq_left2right,
       dir_is_complement;
} Sequenceobject;

static void ft_sequenceobject_init(Sequenceobject *seq,
                                   GtExtendCharAccess extend_char_access_mode,
                                   const GtEncseq *encseq,
                                   bool rightextension,
                                   GtReadmode readmode,
                                   GtUword seqstartpos,
                                   GtUword startpos,
                                   GtUword len,
                                   GtEncseqReader *encseq_r,
                                   GtAllocatedMemory *sequence_cache,
                                   const GtUchar *bytesequence,
                                   GtUword totallength
                                   )
{
  gt_assert(seq != NULL);
  seq->encseq = NULL;
  seq->encseqreader = NULL;
  seq->twobitencoding = NULL;
  seq->cache_ptr = NULL;
  seq->sequence_cache = NULL;
  seq->bytesequenceptr = NULL;
  seq->seqstartpos = seqstartpos;
  gt_assert(seqstartpos <= startpos);
  seq->offset = GT_EXTEND_OFFSET(rightextension,
                                 readmode,
                                 totallength,
                                 seqstartpos,
                                 startpos,
                                 len);
  seq->read_seq_left2right = GT_EXTEND_READ_SEQ_LEFT2RIGHT(rightextension,
                                                           readmode);
  if (encseq != NULL && extend_char_access_mode == GT_EXTEND_CHAR_ACCESS_ANY &&
      gt_encseq_has_twobitencoding(encseq) && gt_encseq_wildcards(encseq) == 0)
  {
    seq->twobitencoding = gt_encseq_twobitencoding_export(encseq);
  }
  if (encseq != NULL && seq->twobitencoding == NULL &&
      (extend_char_access_mode == GT_EXTEND_CHAR_ACCESS_ANY ||
       extend_char_access_mode == GT_EXTEND_CHAR_ACCESS_ENCSEQ_READER))
  {
    GtUword full_totallength = gt_encseq_total_length(encseq);
    gt_encseq_reader_reinit_with_readmode(encseq_r, encseq,
                                          seq->read_seq_left2right
                                            ? GT_READMODE_FORWARD
                                            : GT_READMODE_REVERSE,
                                          seq->read_seq_left2right
                                            ?  seq->offset
                                            : GT_REVERSEPOS(full_totallength,
                                                            seq->offset));
    seq->encseqreader = encseq_r;
    gt_assert(seq->encseqreader != NULL && sequence_cache != NULL);
    seq->sequence_cache = sequence_cache;
    seq->cache_ptr = (GtUchar *) sequence_cache->space;
    seq->min_access_pos = GT_UWORD_MAX; /* undefined */
    seq->cache_num_positions = 0;
  }
  if (encseq != NULL && seq->twobitencoding == NULL &&
      seq->encseqreader == NULL &&
      (extend_char_access_mode == GT_EXTEND_CHAR_ACCESS_ANY ||
       extend_char_access_mode == GT_EXTEND_CHAR_ACCESS_ENCSEQ))
  {
    seq->encseq = encseq;
  }
  if (extend_char_access_mode == GT_EXTEND_CHAR_ACCESS_DIRECT)
  {
    gt_assert(seq->twobitencoding == NULL && seq->encseqreader == NULL &&
              seq->encseq == NULL);
    seq->bytesequenceptr = bytesequence;
  }
  seq->substringlength = len;
  seq->totallength = totallength;
  seq->dir_is_complement = GT_ISDIRCOMPLEMENT(readmode) ? true : false;
  gt_assert(seq->twobitencoding != NULL || seq->encseqreader != NULL ||
            seq->encseq != NULL || seq->bytesequenceptr != NULL);
}

static GtUchar gt_twobitencoding_char_at_pos(
                                      const GtTwobitencoding *twobitencoding,
                                      GtUword pos)
{
  return (twobitencoding[GT_DIVBYUNITSIN2BITENC(pos)] >>
          GT_MULT2(GT_UNITSIN2BITENC - 1 - GT_MODBYUNITSIN2BITENC(pos))) & 3;
}

static GtUchar ft_sequenceobject_get_char(Sequenceobject *seq,GtUword idx)
{
  GtUchar cc;
  GtUword accesspos;

  if (seq->twobitencoding != NULL)
  {
    gt_assert (seq->read_seq_left2right || seq->offset >= idx);
    accesspos = seq->read_seq_left2right ? seq->offset + idx
                                         : seq->offset - idx;
    gt_assert(accesspos < seq->seqstartpos + seq->totallength);
    cc = gt_twobitencoding_char_at_pos(seq->twobitencoding, accesspos);
    return seq->dir_is_complement ? GT_COMPLEMENTBASE(cc) : cc;
  }
  if (seq->encseqreader != NULL)
  {
    gt_assert(idx < seq->substringlength);
    if (idx >= seq->cache_num_positions)
    {
      GtUword idx, tostore;
      const GtUword addamount = 256UL;

      tostore = MIN(seq->cache_num_positions + addamount,
                    seq->substringlength);
      if (tostore > seq->sequence_cache->allocated)
      {
        seq->sequence_cache->allocated += addamount;
        seq->sequence_cache->space
          = gt_realloc(seq->sequence_cache->space,
                       sizeof (GtUchar) * seq->sequence_cache->allocated);
        seq->cache_ptr = (GtUchar *) seq->sequence_cache->space;
      }
      for (idx = seq->cache_num_positions; idx < tostore; idx++)
      {
        seq->cache_ptr[idx]
          = gt_encseq_reader_next_encoded_char(seq->encseqreader);
      }
      seq->cache_num_positions = tostore;
    }
    gt_assert(seq->cache_ptr != NULL && idx < seq->cache_num_positions);
    return seq->cache_ptr[idx];
  }
  accesspos = seq->read_seq_left2right ? seq->offset + idx : seq->offset - idx;
  if (seq->encseq != NULL)
  {
    cc = gt_encseq_get_encoded_char(seq->encseq,accesspos,
                                    GT_READMODE_FORWARD);
  } else
  {
    gt_assert(seq->bytesequenceptr != NULL);
    cc = seq->bytesequenceptr[accesspos];
  }
  if (seq->dir_is_complement && !ISSPECIAL(cc))
  {
    return GT_COMPLEMENTBASE(cc);
  }
  return cc;
}

#undef SKDEBUG
#ifdef SKDEBUG
static char *gt_ft_sequencebject_get(Sequenceobject *seq)
{
  GtUword idx;
  char *buffer;
  char *map = "acgt";

  gt_assert(seq != NULL);
  buffer = gt_malloc(sizeof *buffer * (seq->substringlength+1));
  for (idx = 0; idx < seq->substringlength; idx++)
  {
    GtUchar cc = ft_sequenceobject_get_char(seq,idx);

    if (cc == WILDCARD)
    {
      buffer[idx] = '#';
    } else
    {
      if (cc == SEPARATOR)
      {
        buffer[idx] = '$';
      } else
      {
        gt_assert(cc < 4);
        buffer[idx] = map[cc];
      }
    }
  }
  buffer[seq->substringlength] = '\0';
  return buffer;
}

static void gt_greedy_show_context(bool rightextension,
                                   Sequenceobject *useq,Sequenceobject *vseq)
{
  char *uptr = gt_ft_sequencebject_get(useq);
  char *vptr = gt_ft_sequencebject_get(vseq);
  printf(">%sextension:\n>%s\n>%s\n",rightextension ? "right" : "left",
         uptr,vptr);
  gt_free(uptr);
  gt_free(vptr);
}
#endif

#else
typedef struct
{
  const GtUchar *sequence_ptr;
  GtUword substringlength;
} Sequenceobject;

static void ft_sequenceobject_init(Sequenceobject *seq,
                                   const GtUchar *ptr,
                                   GtUword seqstartpos,
                                   GtUword startpos,
                                   GtUword len)
{
  gt_assert(seq != NULL);
  seq->sequence_ptr = ptr + seqstartpos + startpos;
  seq->substringlength = len;
}
#endif

#define FRONT_DIAGONAL(FRONTPTR) (GtWord) ((FRONTPTR) - midfront)

static bool ft_sequenceobject_symbol_match(Sequenceobject *useq,
                                           GtUword upos,
                                           Sequenceobject *vseq,
                                           GtUword vpos)
{
#ifndef OUTSIDE_OF_GT
  GtUchar cu = ft_sequenceobject_get_char(useq,upos);
  if (ISSPECIAL(cu))
  {
    return false;
  }
  return cu == ft_sequenceobject_get_char(vseq,vpos) ? true : false;
#else
  GtUchar cu = useq->sequence_ptr[upos];
  return cu == vseq->sequence_ptr[vpos] ? true : false;
#endif
}

static void inline front_prune_add_matches(Frontvalue *midfront,
                                           Frontvalue *fv,
                                           uint64_t mask,
                                           Sequenceobject *useq,
                                           Sequenceobject *vseq)
{
  GtUword upos, vpos;

  fv->localmatch_count = 0;
  for (upos = fv->row, vpos = fv->row + FRONT_DIAGONAL(fv);
       upos < useq->substringlength && vpos < vseq->substringlength &&
       ft_sequenceobject_symbol_match(useq,upos,vseq,vpos);
       upos++, vpos++)
  {
    fv->localmatch_count++;
    if (!(fv->matchhistory & mask))
    {
      gt_assert(fv->matchhistory_count < INT8_MAX);
      fv->matchhistory_count++;
    }
    fv->matchhistory = (fv->matchhistory << 1) | (uint64_t) 1;
  }
  fv->row += fv->localmatch_count;
}

static GtUword front_next_inplace(Frontvalue *midfront,
                                  Frontvalue *lowfront,
                                  Frontvalue *highfront,
                                  GtUword history,
                                  Sequenceobject *useq,
                                  Sequenceobject *vseq)
{
  GtUword alignedlen, maxalignedlen;
  const uint64_t mask = ((uint64_t) 1) << (history-1);
  Frontvalue bestfront, insertion_value, replacement_value, *frontptr;

  insertion_value = *lowfront; /* from previous diag -(d-1) => -d => DELETION */
  bestfront = insertion_value;
  bestfront.row++;
  UPDATE_MATCH_HISTORY(bestfront.matchhistory_count,bestfront.matchhistory);
  *lowfront = bestfront;
  lowfront->backreference = FT_EOP_DELETION;
  front_prune_add_matches(midfront,lowfront,mask,useq,vseq);
  maxalignedlen = GT_MULT2(lowfront->row) + FRONT_DIAGONAL(lowfront);

  replacement_value = *(lowfront+1);
  if (bestfront.row < replacement_value.row + 1)
  {
    bestfront = replacement_value;
    bestfront.backreference = FT_EOP_DELETION;
    bestfront.row++;
    UPDATE_MATCH_HISTORY(bestfront.matchhistory_count,bestfront.matchhistory);
  } else
  {
    bestfront.backreference = FT_EOP_REPLACEMENT;
    if (bestfront.row == replacement_value.row + 1)
    {
      bestfront.backreference |= FT_EOP_DELETION;
    }
  }
  *(lowfront+1) = bestfront;
  front_prune_add_matches(midfront,lowfront + 1,mask,useq,vseq);
  alignedlen = GT_MULT2((lowfront+1)->row) + FRONT_DIAGONAL(lowfront + 1);
  if (maxalignedlen < alignedlen)
  {
    maxalignedlen = alignedlen;
  }
  for (frontptr = lowfront+2; frontptr <= highfront; frontptr++)
  {
    bestfront = insertion_value;
    bestfront.backreference = FT_EOP_INSERTION;
    if (frontptr <= highfront - 1)
    {
      if (bestfront.row < replacement_value.row + 1)
      {
        bestfront = replacement_value;
        bestfront.backreference = FT_EOP_REPLACEMENT;
        bestfront.row++;
      } else
      {
        if (bestfront.row == replacement_value.row + 1)
        {
          bestfront.backreference |= FT_EOP_REPLACEMENT;
        }
      }
    }
    if (frontptr <= highfront - 2)
    {
      if (bestfront.row < frontptr->row + 1)
      {
        bestfront = *frontptr;
        bestfront.backreference = FT_EOP_DELETION;
        bestfront.row++;
      } else
      {
        if (bestfront.row == frontptr->row + 1)
        {
          bestfront.backreference |= FT_EOP_DELETION;
        }
      }
    }
    UPDATE_MATCH_HISTORY(bestfront.matchhistory_count,bestfront.matchhistory);
    if (frontptr < highfront)
    {
      insertion_value = replacement_value;
      replacement_value = *frontptr;
    }
    *frontptr = bestfront;
    front_prune_add_matches(midfront,frontptr,mask,useq,vseq);
    alignedlen = GT_MULT2(frontptr->row) + FRONT_DIAGONAL(frontptr);
    if (maxalignedlen < alignedlen)
    {
      maxalignedlen = alignedlen;
    }
  }
  return maxalignedlen;
}

static GtUword front_second_inplace(Frontvalue *midfront,
                                    Frontvalue *lowfront,
                                    GtUword history,
                                    Sequenceobject *useq,
                                    Sequenceobject *vseq)
{
  GtUword alignedlen, maxalignedlen;
  const uint64_t mask = ((uint64_t) 1) << (history-1);

  *(lowfront+1) = *(lowfront+2) = *lowfront;
  lowfront->row++;
  lowfront->backreference = FT_EOP_DELETION;
  UPDATE_MATCH_HISTORY(lowfront->matchhistory_count,lowfront->matchhistory);
  front_prune_add_matches(midfront,lowfront,mask,useq,vseq);
  maxalignedlen = GT_MULT2(lowfront->row) + FRONT_DIAGONAL(lowfront);

  (lowfront+1)->row++;
  (lowfront+1)->backreference = FT_EOP_REPLACEMENT;
  UPDATE_MATCH_HISTORY((lowfront+1)->matchhistory_count,
                       (lowfront+1)->matchhistory);
  front_prune_add_matches(midfront,lowfront + 1,mask,useq,vseq);
  alignedlen = GT_MULT2((lowfront+1)->row) + FRONT_DIAGONAL(lowfront + 1);
  if (maxalignedlen < alignedlen)
  {
    maxalignedlen = alignedlen;
  }

  (lowfront+2)->backreference = FT_EOP_INSERTION;
  UPDATE_MATCH_HISTORY((lowfront+2)->matchhistory_count,
                       (lowfront+2)->matchhistory);
  front_prune_add_matches(midfront,lowfront + 2,mask,useq,vseq);
  alignedlen = GT_MULT2((lowfront+2)->row) + FRONT_DIAGONAL(lowfront + 2);
  if (maxalignedlen < alignedlen)
  {
    maxalignedlen = alignedlen;
  }
  return maxalignedlen;
}

#undef TRIM_INFO_OUT
#ifdef TRIM_INFO_OUT
static bool trimthisentry(GtUword distance,
                          Rowvaluetype row,
                          GtWord diagonal,
                          GtUword minlenforhistorycheck,
                          Matchcounttype matchhistory_count,
                          GtUword minmatchnum,
                          GtUword minlenfrommaxdiff)
{
  GtUword alignedlen = GT_MULT2(row) + diagonal;

  if (alignedlen >= minlenforhistorycheck && matchhistory_count < minmatchnum)
  {
    printf(GT_WD "&" GT_WU "&%u&1: matches=%d < " GT_WU "=minmatches\n",
              diagonal,distance,row,(int) matchhistory_count,minmatchnum);
    return true;
  }
  if (alignedlen < minlenfrommaxdiff)
  {
    printf(GT_WD "&" GT_WU "&%u&2: i'+j'=" GT_WU "<" GT_WU "=i+j-lag\n",
              diagonal,distance,row,alignedlen,minlenfrommaxdiff);
    return true;
  }
  printf(GT_WD "&" GT_WU "&%u\n", diagonal,distance,row);
  return false;
}
#else
static bool trimthisentry(Rowvaluetype row,
                          GtWord diagonal,
                          GtUword minlenforhistorycheck,
                          Matchcounttype matchhistory_count,
                          GtUword minmatchnum,
                          GtUword minlenfrommaxdiff)
{
  GtUword alignedlen = GT_MULT2(row) + diagonal;

  if (alignedlen >= minlenforhistorycheck && matchhistory_count < minmatchnum)
  {
    return true;
  }
  if (alignedlen < minlenfrommaxdiff)
  {
    return true;
  }
  return false;
}
#endif

static GtUword trim_front(bool upward,
#ifdef TRIM_INFO_OUT
                          GtUword distance,
#endif
                          GtUword ulen,
                          GtUword vlen,
                          GtUword minmatchnum,
                          GtUword minlenforhistorycheck,
                          GtUword minlenfrommaxdiff,
                          const Frontvalue *midfront,
                          const Frontvalue *from,
                          const Frontvalue *stop)
{
  const Frontvalue *frontptr;
  GtUword trim = 0;

  gt_assert ((upward && from < stop) || (!upward && stop < from));
  for (frontptr = from; frontptr != stop; frontptr = upward ? (frontptr + 1)
                                                            : (frontptr - 1))
  {
    if (frontptr->row > ulen ||
        frontptr->row + FRONT_DIAGONAL(frontptr) > vlen ||
        trimthisentry(
#ifdef TRIM_INFO_OUT
                        distance,
#endif
                        frontptr->row,
                        FRONT_DIAGONAL(frontptr),
                        minlenforhistorycheck,
                        frontptr->matchhistory_count,
                        minmatchnum,
                        minlenfrommaxdiff))
    {
      trim++;
    } else
    {
      break;
    }
  }
  return trim;
}

static void frontspace_check(GT_UNUSED const Frontvalue *from,
                             GT_UNUSED const Frontvalue *to,
                             GT_UNUSED const Frontvalue *ptr)
{
  gt_assert (ptr >= from && ptr <= to);
}

static Frontvalue *frontspace_allocate(GtUword minsizeforshift,
                                       GtUword trimleft,
                                       GtUword valid,
                                       GtAllocatedMemory *fs)
{
  if (trimleft - fs->offset + valid >= fs->allocated)
  {
    fs->allocated = 255UL + MAX(fs->allocated * 1.2,
                                           trimleft - fs->offset + valid);
    gt_assert(fs->allocated > trimleft - fs->offset + valid);
    fs->space = gt_realloc(fs->space,sizeof (Frontvalue) * fs->allocated);
    gt_assert(fs->space != NULL);
  }
  gt_assert(trimleft >= fs->offset);
  if (trimleft - fs->offset > MAX(valid,minsizeforshift))
  {
    memcpy(fs->space,((Frontvalue *) fs->space) + trimleft - fs->offset,
           sizeof (Frontvalue) * valid);
    fs->offset = trimleft;
  }
  return ((Frontvalue *) fs->space) - fs->offset;
}

static void update_trace_and_polished(Polished_point *best_polished_point,
#ifndef OUTSIDE_OF_GT
                                      GtUword *minrow,
                                      GtUword *mincol,
#endif
                                      Fronttrace *front_trace,
                                      const Polishing_info *pol_info,
                                      GtUword distance,
                                      GtUword trimleft,
                                      Frontvalue *midfront,
                                      Frontvalue *lowfront,
                                      Frontvalue *highfront)
{
  const Frontvalue *frontptr;
  uint64_t lsb;

#ifndef OUTSIDE_OF_GT
  *minrow = GT_UWORD_MAX;
  *mincol = GT_UWORD_MAX;
#endif
  for (frontptr = lowfront; frontptr <= highfront; frontptr++)
  {
    GtUword alignedlen = GT_MULT2(frontptr->row) + FRONT_DIAGONAL(frontptr);

#ifndef OUTSIDE_OF_GT
    GtUword currentcol;

    if (*minrow > frontptr->row)
    {
      *minrow = frontptr->row;
    }
    gt_assert(FRONT_DIAGONAL(frontptr) >= 0 ||
              frontptr->row >= -FRONT_DIAGONAL(frontptr));
    currentcol = frontptr->row + FRONT_DIAGONAL(frontptr);
    if (*mincol > currentcol)
    {
      *mincol = currentcol;
    }
#endif
    lsb = frontptr->matchhistory & pol_info->mask;
    if (HISTORY_IS_POLISHED(pol_info,frontptr->matchhistory,lsb) &&
        alignedlen > best_polished_point->alignedlen)
    {
      best_polished_point->alignedlen = alignedlen;
      best_polished_point->row = frontptr->row;
      best_polished_point->distance = distance;
      best_polished_point->trimleft = trimleft;
    }
    if (front_trace != NULL)
    {
      front_trace_add_trace(front_trace,frontptr->backreference,
                            frontptr->localmatch_count);
    }
  }
}

GtUword front_prune_edist_inplace(
#ifndef OUTSIDE_OF_GT
                         bool rightextension,
                         GtAllocatedMemory *frontspace,
#endif
                         Trimstat *trimstat,
                         Polished_point *best_polished_point,
                         Fronttrace *front_trace,
                         const Polishing_info *pol_info,
                         GtUword history,
                         GtUword minmatchnum,
                         GtUword maxalignedlendifference,
                         FTsequenceResources *ufsr,
                         GtUword ustart,
                         GtUword ulen,
                         GtUword vseqstartpos,
                         FTsequenceResources *vfsr,
                         GtUword vstart,
                         GtUword vlen)
{
  const GtUword sumseqlength = ulen + vlen,
                minsizeforshift = sumseqlength/1000,
                minlenforhistorycheck = GT_MULT2(history);
  /* so the space for allocating the fronts is
     sizeof (Frontvalue) * ((m+n)/1000 + maxvalid), where maxvalid is a small
     constant. */
  GtUword distance, trimleft = 0, valid = 1UL, maxvalid = 0, sumvalid = 0;
  const uint64_t mask = ((uint64_t) 1) << (history-1);
  Frontvalue *validbasefront;
  bool diedout = false;
  Sequenceobject useq, vseq;

#ifdef OUTSIDE_OF_GT
  GtAllocatedMemory *frontspace = gt_malloc(sizeof *frontspace);
  frontspace->space = NULL;
  frontspace->allocated = 0;
  frontspace->offset = 0;
  ft_sequenceobject_init(&useq,useqptr,0,ustart,ulen);
  ft_sequenceobject_init(&vseq,vseqptr,vseqstartpos,vstart,vlen);
#else
  /*
  printf("%sextension:useq->readmode=%s,vseq->readmode=%s\n",
          rightextension ? "right" : "left",
          gt_readmode_show(ufsr->readmode),gt_readmode_show(ufsr->readmode));
  */
  ft_sequenceobject_init(&useq,
                         ufsr->extend_char_access,
                         ufsr->encseq,
                         rightextension,
                         ufsr->readmode,
                         0,
                         ustart,
                         ulen,
                         ufsr->encseq_r,
                         ufsr->sequence_cache,
                         NULL,
                         ufsr->totallength);
  ft_sequenceobject_init(&vseq,
                         vfsr->extend_char_access,
                         vfsr->encseq,
                         rightextension,
                         vfsr->readmode,
                         vseqstartpos,
                         vstart,
                         vlen,
                         vfsr->encseq_r,
                         vfsr->sequence_cache,
                         vfsr->bytesequence,
                         vfsr->totallength);
#ifdef SKDEBUG
  gt_greedy_show_context(rightextension,&useq,&vseq);
#endif
  frontspace->offset = 0;
#endif
#ifdef TRIM_INFO_OUT
  printf("regionalquality(minmatchnum)=" GT_WU "\n",minmatchnum);
#endif
  for (distance = 0, valid = 1UL; /* Nothing */; distance++, valid += 2)
  {
    GtUword trim, maxalignedlen, minlenfrommaxdiff;

#ifdef TRIM_INFO_OUT
    printf("distance=" GT_WU ",full=" GT_WU ",trimleft=" GT_WU
           ",valid=" GT_WU "\n",distance,
                  GT_MULT2(distance) + 1,
                  trimleft,valid);
#endif
    gt_assert(valid <= GT_MULT2(distance) + 1);
    sumvalid += valid;
    if (maxvalid < valid)
    {
      maxvalid = valid;
    }
    validbasefront = frontspace_allocate(minsizeforshift,trimleft,valid,
                                         frontspace);
    if (distance == 0)
    {
      validbasefront->row = 0;
      validbasefront->matchhistory = 0;
      validbasefront->matchhistory_count = 0;
      validbasefront->backreference = 0; /* No back reference */
      front_prune_add_matches(validbasefront + distance,validbasefront,mask,
                              &useq,&vseq);
      maxalignedlen = GT_MULT2(validbasefront->row);
    } else
    {
      gt_assert(valid >= 3UL);
      frontspace_check((const Frontvalue *) frontspace->space,
                       ((const Frontvalue *) frontspace->space)
                        + frontspace->allocated - 1,
                       validbasefront + trimleft);
      frontspace_check((const Frontvalue *) frontspace->space,
                       ((const Frontvalue *) frontspace->space)
                         + frontspace->allocated - 1,
                       validbasefront + trimleft + valid - 1);
      if (valid == 3UL)
      {
        maxalignedlen
          = front_second_inplace(validbasefront + distance,
                                 validbasefront + trimleft,
                                 history,
                                 &useq,
                                 &vseq);
      } else
      {
        maxalignedlen
          = front_next_inplace(validbasefront + distance,
                               validbasefront + trimleft,
                               validbasefront + trimleft + valid - 1,
                               history,
                               &useq,
                               &vseq);
      }
    }
    gt_assert(valid > 0);
    minlenfrommaxdiff = maxalignedlen >= maxalignedlendifference
                          ? maxalignedlen - maxalignedlendifference
                          : 0;
#ifdef TRIM_INFO_OUT
    printf("maxalignedlen=" GT_WU ",maxlenfrommaxdiff=" GT_WU "\n",
           maxalignedlen,minlenfrommaxdiff);
#endif
    trim = trim_front(true,
#ifdef TRIM_INFO_OUT
                      distance,
#endif
                      ulen,
                      vlen,
                      minmatchnum,
                      minlenforhistorycheck,
                      minlenfrommaxdiff,
                      validbasefront + distance,
                      validbasefront + trimleft,
                      validbasefront + trimleft + valid);
#ifdef TRIM_INFO_OUT
    printf("trim on left=" GT_WU "\n",trim);
#endif
    if (trim > 0)
    {
      trimleft += trim;
      gt_assert(valid >= trim);
      valid -= trim;
    }
    if (valid > 0)
    {
      trim = trim_front(false,
#ifdef TRIM_INFO_OUT
                        distance,
#endif
                        ulen,
                        vlen,
                        minmatchnum,
                        minlenforhistorycheck,
                        minlenfrommaxdiff,
                        validbasefront + distance,
                        validbasefront + trimleft + valid - 1,
                        validbasefront + trimleft - 1);
#ifdef TRIM_INFO_OUT
      printf("trim on right=" GT_WU "\n",trim);
#endif
      gt_assert(trim < valid);
      if (trim > 0)
      {
        gt_assert(valid >= trim);
        valid -= trim;
      }
    }
    if (valid == 0)
    {
      diedout = true;
      break;
    }
    if (front_trace != NULL)
    {
      front_trace_add_gen(front_trace,trimleft,valid);
    }
    update_trace_and_polished(best_polished_point,
#ifndef OUTSIDE_OF_GT
                              &useq.min_access_pos,
                              &vseq.min_access_pos,
#endif
                              front_trace,
                              pol_info,
                              distance,
                              trimleft,
                              validbasefront + distance,
                              validbasefront + trimleft,
                              validbasefront + trimleft + valid - 1);
    if ((vlen > ulen && vlen - ulen <= distance) ||
        (vlen <= ulen && ulen - vlen <= distance))
    {
      if (distance + vlen - ulen >= trimleft &&
          distance + vlen - ulen <= trimleft + valid - 1 &&
          validbasefront[distance + vlen - ulen].row == ulen)
      {
        break;
      }
    }
    if (distance >= sumseqlength)
    {
      break;
    }
  }
  trimstat_add(trimstat,diedout,sumvalid,maxvalid,distance,
               sizeof (Frontvalue) * frontspace->allocated,
#ifndef OUTSIDE_OF_GT
               useq.sequence_cache != NULL &&
               vseq.sequence_cache != NULL ? MAX(useq.sequence_cache->allocated,
                                                 vseq.sequence_cache->allocated)
                                           : 0
#else
               0
#endif
              );
  return diedout ? sumseqlength + 1 : distance;
}
