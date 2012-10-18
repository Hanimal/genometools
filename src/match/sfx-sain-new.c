/*
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include <limits.h>
#include "core/minmax.h"
#include "core/unused_api.h"
#include "core/timer_api.h"
#include "sfx-linlcp.h"
#include "sfx-sain.h"

#define GT_SSTARLENGTH_MAX 50

#define GT_SAIN_SHOWTIMER(DESC)\
        if (timer != NULL)\
        {\
          gt_timer_show_progress(timer,DESC,stdout);\
        }

typedef enum
{
  GT_SAIN_PLAINSEQ,
  GT_SAIN_INTSEQ,
  GT_SAIN_ENCSEQ
} GtSainSeqtype;

typedef struct
{
  unsigned long totallength,
                numofchars,
                startoccupied,
                *bucketsize,
                *bucketfillptr;
  union
  {
    const GtEncseq *encseq;
    const GtUchar *plainseq;
    const unsigned long *array;
  } seq;
  GtSainSeqtype seqtype;
  bool bucketfillptrpoints2suftab,
       bucketsizepoints2suftab;
} GtSainseq;

static GtSainseq *gt_sain_seq_new_from_encseq(const GtEncseq *encseq)
{
  unsigned long idx;
  GtSainseq *sainseq = gt_malloc(sizeof *sainseq);

  sainseq->seqtype = GT_SAIN_ENCSEQ;
  sainseq->seq.encseq = encseq;
  sainseq->totallength = gt_encseq_total_length(encseq);
  sainseq->numofchars = (unsigned long) gt_encseq_alphabetnumofchars(encseq);
  sainseq->bucketsize = gt_malloc(sizeof (*sainseq->bucketsize) *
                                  sainseq->numofchars);
  sainseq->bucketfillptr = gt_malloc(sizeof (*sainseq->bucketfillptr) *
                                     sainseq->numofchars);
  sainseq->bucketfillptrpoints2suftab = false;
  sainseq->bucketsizepoints2suftab = false;
  for (idx = 0; idx<sainseq->numofchars; idx++)
  {
    sainseq->bucketsize[idx] = gt_encseq_charcount(encseq,(GtUchar) idx);
  }
  return sainseq;
}

static GtSainseq *gt_sain_seq_new_from_plainseq(const GtUchar *plainseq,
                                                unsigned long len)
{
  unsigned long idx;
  GtSainseq *sainseq = gt_malloc(sizeof *sainseq);

  sainseq->seqtype = GT_SAIN_PLAINSEQ;
  sainseq->seq.plainseq = plainseq;
  sainseq->totallength = len;
  sainseq->numofchars = UCHAR_MAX+1;
  sainseq->bucketsize = gt_malloc(sizeof (*sainseq->bucketsize) *
                                  sainseq->numofchars);
  sainseq->bucketfillptr = gt_malloc(sizeof (*sainseq->bucketfillptr) *
                                     sainseq->numofchars);
  sainseq->bucketfillptrpoints2suftab = false;
  sainseq->bucketsizepoints2suftab = false;
  for (idx = 0; idx<sainseq->numofchars; idx++)
  {
    sainseq->bucketsize[idx] = 0;
  }
  for (idx = 0; idx<sainseq->totallength; idx++)
  {
    sainseq->bucketsize[sainseq->seq.plainseq[idx]]++;
  }
  return sainseq;
}

static void gt_sain_determinebucketsize(GtSainseq *sainseq)
{
  unsigned long idx;

  for (idx = 0; idx < sainseq->numofchars; idx++)
  {
    sainseq->bucketsize[idx] = 0;
  }
  for (idx = 0; idx<sainseq->totallength; idx++)
  {
    gt_assert(sainseq->seq.array[idx] < sainseq->numofchars);
    sainseq->bucketsize[sainseq->seq.array[idx]]++;
  }
}

static GtSainseq *gt_sain_seq_new_from_array(unsigned long *arr,
                                             unsigned long len,
                                             unsigned long numofchars,
                                             unsigned long *suftab,
                                             unsigned long suftabentries)
{
  unsigned long maxused = GT_MULT2(len);
  GtSainseq *sainseq = gt_malloc(sizeof *sainseq);

  sainseq->seqtype = GT_SAIN_INTSEQ;
  sainseq->seq.array = arr;
  sainseq->totallength = len;
  sainseq->numofchars = numofchars;
  sainseq->bucketfillptrpoints2suftab = false;
  sainseq->bucketsizepoints2suftab = false;
  sainseq->startoccupied = suftabentries;
  if (maxused < suftabentries)
  {
    if (suftabentries - maxused >= numofchars)
    {
      /*printf("bucketsize: reclaim suftab at %lu for %lu elements\n",
              suftabentries - numofchars,numofchars);*/
      sainseq->bucketsize = suftab + suftabentries - numofchars;
      sainseq->bucketsizepoints2suftab = true;
      sainseq->startoccupied = suftabentries - numofchars;
    }
    if (suftabentries - maxused >= GT_MULT2(numofchars))
    {
      /*printf("bucketfillptr: reclaim suftab at %lu for %lu elements\n",
              suftabentries - GT_MULT2(numofchars),GT_MULT2(numofchars));*/
      sainseq->bucketfillptr = suftab + suftabentries - GT_MULT2(numofchars);
      sainseq->bucketfillptrpoints2suftab = true;
      sainseq->startoccupied = suftabentries - GT_MULT2(numofchars);
    }
  }
  if (!sainseq->bucketfillptrpoints2suftab)
  {
    sainseq->bucketfillptr = gt_malloc(sizeof (*sainseq->bucketfillptr) *
                                       numofchars);
    printf("bucketfillptr requires %lu bytes\n",
           (unsigned long) sizeof (*sainseq->bucketfillptr) * numofchars);
  }
  if (!sainseq->bucketsizepoints2suftab)
  {
    sainseq->bucketsize = gt_malloc(sizeof (*sainseq->bucketsize) * numofchars);
    printf("bucketsize requires %lu bytes\n",
           (unsigned long) sizeof (*sainseq->bucketsize) * numofchars);
  }
  gt_sain_determinebucketsize(sainseq);
  return sainseq;
}

static unsigned long gt_sain_seq_delete(GtSainseq *sainseq)
{
  unsigned long ret = ULONG_MAX;

  if (sainseq != NULL)
  {
    if (!sainseq->bucketfillptrpoints2suftab)
    {
      gt_free(sainseq->bucketfillptr);
    }
    if (!sainseq->bucketsizepoints2suftab)
    {
      gt_free(sainseq->bucketsize);
    }
    if (sainseq->bucketsizepoints2suftab ||
        sainseq->bucketsizepoints2suftab)
    {
      ret = sainseq->startoccupied;
    }
    gt_free(sainseq);
  }
  return ret;
}

static unsigned long gt_sain_countcharaccess = 0;

static unsigned long gt_sain_seq_getchar(const GtSainseq *sainseq,
                                         unsigned long position)
{
  gt_assert(position < sainseq->totallength);
  gt_sain_countcharaccess++;
  switch (sainseq->seqtype)
  {
    case GT_SAIN_PLAINSEQ:
      return (unsigned long) sainseq->seq.plainseq[position];
    case GT_SAIN_INTSEQ:
      return sainseq->seq.array[position];
    case GT_SAIN_ENCSEQ:
      {
        GtUchar cc = gt_encseq_get_encoded_char(sainseq->seq.encseq,
                                                position,
                                                GT_READMODE_FORWARD);
        return ISSPECIAL(cc) ? GT_UNIQUEINT(position) : (unsigned long) cc;
      }
  }
  /*@ignore@*/
  return 0;
  /*@end@*/
}

static unsigned long gt_sain_seq_getchar_nospecial(const GtSainseq *sainseq,
                                                   unsigned long position)
{
  gt_assert(position < sainseq->totallength);
  gt_sain_countcharaccess++;
  switch (sainseq->seqtype)
  {
    case GT_SAIN_PLAINSEQ:
      return (unsigned long) sainseq->seq.plainseq[position];
    case GT_SAIN_INTSEQ:
      return sainseq->seq.array[position];
    case GT_SAIN_ENCSEQ:
    return (unsigned long)
           gt_encseq_get_encoded_char_nospecial(sainseq->seq.encseq,
                                                position,
                                                GT_READMODE_FORWARD);
  }
  /*@ignore@*/
  return 0;
  /*@end@*/
}

static void gt_sain_endbuckets(GtSainseq *sainseq)
{
  unsigned long charidx;

  sainseq->bucketfillptr[0] = sainseq->bucketsize[0];
  for (charidx = 1UL; charidx < sainseq->numofchars; charidx++)
  {
    sainseq->bucketfillptr[charidx] = sainseq->bucketfillptr[charidx-1] +
                                      sainseq->bucketsize[charidx];
  }
}

static void gt_sain_startbuckets(GtSainseq *sainseq)
{
  unsigned long charidx;

  sainseq->bucketfillptr[0] = 0;
  for (charidx = 1UL; charidx < sainseq->numofchars; charidx++)
  {
    sainseq->bucketfillptr[charidx] = sainseq->bucketfillptr[charidx-1] +
                                      sainseq->bucketsize[charidx-1];
  }
}

typedef struct
{
#ifndef NDEBUG
  GtBitsequence *isStype;
#endif
  unsigned long countStype,
                countSstartype,
                totalSstarlength,
                longerthanmax,
                lendist[GT_SSTARLENGTH_MAX+1];
  GtSainseq *sainseq;
} GtSaininfo;

static GtSaininfo *gt_sain_info_new(GtSainseq *sainseq,unsigned long *suftab,
                                    GT_UNUSED unsigned long nonspecialentries)
{
  unsigned long position,
                idx,
                nextcc,
                nextSstartypepos,
                *fillptr = sainseq->bucketfillptr;
  bool nextisStype = true;
  GtSaininfo *saininfo;

  saininfo = gt_malloc(sizeof *saininfo);
  saininfo->sainseq = sainseq;
  saininfo->totalSstarlength = 0;
  saininfo->countStype = 1UL;
  saininfo->countSstartype = 0;
  saininfo->longerthanmax = 0;
#ifndef NDEBUG
  GT_INITBITTAB(saininfo->isStype,saininfo->sainseq->totallength+1);
  GT_SETIBIT(saininfo->isStype,saininfo->sainseq->totallength);
#endif
  nextSstartypepos = saininfo->sainseq->totallength;
  nextcc = GT_UNIQUEINT(saininfo->sainseq->totallength);
  for (idx = 0; idx <= (unsigned long) GT_SSTARLENGTH_MAX; idx++)
  {
    saininfo->lendist[idx] = 0;
  }
  gt_sain_endbuckets(saininfo->sainseq);
  gt_assert(fillptr[saininfo->sainseq->numofchars-1] == nonspecialentries);
  for (position = saininfo->sainseq->totallength-1; /* Nothing */; position--)
  {
    bool currentisStype;
    unsigned long currentcc;

    currentcc = gt_sain_seq_getchar(saininfo->sainseq,position);
    if (currentcc < nextcc || (currentcc == nextcc && nextisStype))
    {
      currentisStype = true;
      saininfo->countStype++;
#ifndef NDEBUG
      GT_SETIBIT(saininfo->isStype,position);
#endif
    } else
    {
      currentisStype = false;
    }
    /*printf("position %lu, char=%lu is %c-type\n",position,currentcc,
                                           currentisStype ? 'S' : 'L');*/
    if (!currentisStype && nextisStype)
    {
      unsigned long currentlen, putidx, cc;

      saininfo->countSstartype++;
      gt_assert(position + 1 < nextSstartypepos);
      currentlen = nextSstartypepos - position;
      saininfo->totalSstarlength += currentlen;
      if (currentlen <= (unsigned long) GT_SSTARLENGTH_MAX)
      {
        saininfo->lendist[currentlen]++;
      } else
      {
        saininfo->longerthanmax++;
      }
      nextSstartypepos = position + 1;
      cc = gt_sain_seq_getchar_nospecial(saininfo->sainseq,position+1);
      gt_assert(fillptr[cc] > 0);
      putidx = --fillptr[cc];
      gt_assert(putidx < nonspecialentries);
      suftab[putidx] = position;
#undef SAINSHOWSTATE
#ifdef SAINSHOWSTATE
      printf("Sstar.suftab[%lu]=%lu\n",putidx,position+1);
#endif
    }
    nextisStype = currentisStype;
    nextcc = currentcc;
    if (position == 0)
    {
      break;
    }
  }
  gt_assert(GT_MULT2(saininfo->countSstartype) <=
            saininfo->sainseq->totallength);
  return saininfo;
}

static void gt_sain_info_delete(GtSaininfo *saininfo)
{
  if (saininfo != NULL)
  {
#ifndef NDEBUG
    gt_free(saininfo->isStype);
#endif
    gt_free(saininfo);
  }
}

#ifdef CRITICAL
static bool gt_sain_info_isSstartype(const GtSaininfo *saininfo,
                                     unsigned long position)
{
  gt_assert(position <= saininfo->sainseq->totallength);
  return position == saininfo->sainseq->totallength ||
         (position > 0 &&
         GT_ISIBITSET(saininfo->isStype,position) &&
         !GT_ISIBITSET(saininfo->isStype,position-1)) ? true : false;
}

unsigned long gt_sain_countDcriticalsubstrings (const GtSaininfo *saininfo,
                                                unsigned long d)
{
  unsigned long int i = 0, j = 0;

  gt_assert(d >= 2UL);
  while (i < saininfo->sainseq->totallength)
  {
    unsigned long h;
    bool isLMS = false;

    for (h = 1UL; h <= d; h++)
    {
      if (gt_sain_info_isSstartype(saininfo,i+h))
      {
        isLMS = true;
        break;
      }
    }
    if (j == 0 && !isLMS)
    {
      i += d;
      continue;
    }
    i = isLMS ? i + h : i + d;
    gt_assert(i>0);
    /*printf("crititical %lu\n",i-1);*/
    j++;
  }
  return j;
}
#endif

static void gt_sain_info_show(const GtSaininfo *saininfo)
{
  unsigned long idx;

  printf("S-type: %lu (%.2f)\n",saininfo->countStype,
                (double) saininfo->countStype/saininfo->sainseq->totallength);
  printf("Sstar-type: %lu (%.2f)\n",saininfo->countSstartype,
              (double) saininfo->countSstartype/saininfo->sainseq->totallength);
  printf("Sstar-type.length: %lu (%.2f)\n",saininfo->totalSstarlength,
              (double) saininfo->totalSstarlength/saininfo->countSstartype);
  for (idx = 0; idx <= (unsigned long) GT_SSTARLENGTH_MAX; idx++)
  {
    if (saininfo->lendist[idx] > 0)
    {
      printf("%lu %lu (%.2f)\n",idx,saininfo->lendist[idx],
               (double) saininfo->lendist[idx]/saininfo->countSstartype);
    }
  }
  if (saininfo->longerthanmax)
  {
    printf(">%d %lu (%.2f)\n",GT_SSTARLENGTH_MAX,saininfo->longerthanmax,
                 (double) saininfo->longerthanmax/saininfo->countSstartype);
  }
#ifdef CRITICAL
  {
    unsigned long d;
    for (d=2UL; d<10UL; d++)
    {
      unsigned long critical = gt_sain_countDcriticalsubstrings (saininfo,d);
      printf("d=%lu,critical=%lu (%.2f)\n",d,critical,
                        (double) critical/saininfo->sainseq->totallength);
    }
  }
#endif
}

#ifndef NDEBUG
static bool gt_sain_info_isStype(const GtSaininfo *saininfo,
                                 unsigned long position)
{
  gt_assert(position <= saininfo->sainseq->totallength);
  return GT_ISIBITSET(saininfo->isStype,position) ? true : false;
}
#endif

static void gt_sain_induceLtypesuffixes1(const GtSaininfo *saininfo,
                                         long *suftab,
                                         unsigned long nonspecialentries)
{
  unsigned long idx,
                *fillptr = saininfo->sainseq->bucketfillptr;

  for (idx = 0; idx < nonspecialentries; idx++)
  {
    long position = suftab[idx];

    if (position > 0)
    {
      unsigned long nextcc;

      nextcc = gt_sain_seq_getchar(saininfo->sainseq,(unsigned long) position);
      if (nextcc < saininfo->sainseq->numofchars)
      {
        unsigned long cc, putidx;

        putidx = fillptr[nextcc]++;
        gt_assert(putidx < nonspecialentries && idx < putidx);
        position--;
        cc = gt_sain_seq_getchar(saininfo->sainseq,(unsigned long) position);
        suftab[putidx] = (cc < nextcc) ? ~position : position;
        suftab[idx] = 0;
#ifdef SAINSHOWSTATE
        printf("L-induce: suftab[%lu]=%ld\n",putidx,suftab[putidx]);
#endif
      }
    } else
    {
      if (position < 0)
      {
        suftab[idx] = ~position;
      }
    }
  }
}

static void gt_sain_induceLtypesuffixes2(const GtSaininfo *saininfo,
                                         long *suftab,
                                         unsigned long nonspecialentries)
{
  unsigned long idx, *fillptr = saininfo->sainseq->bucketfillptr;

  for (idx = 0; idx < nonspecialentries; idx++)
  {
    long position = suftab[idx];

    suftab[idx] = ~position;
    if (position > 0)
    {
      unsigned long nextcc;

      position--;
      nextcc = gt_sain_seq_getchar(saininfo->sainseq,(unsigned long) position);
      if (nextcc < saininfo->sainseq->numofchars)
      {
        unsigned long putidx;

        putidx = fillptr[nextcc]++;
        gt_assert(putidx < nonspecialentries && idx < putidx);
        if (position > 0)
        {
          unsigned long cc = gt_sain_seq_getchar(saininfo->sainseq,
                                                 (unsigned long) (position-1));
          suftab[putidx] = (cc < nextcc) ? ~position : position;
        } else
        {
          suftab[putidx] = position;
        }
#ifdef SAINSHOWSTATE
        printf("L-induce: suftab[%lu]=%ld\n",putidx,suftab[putidx]);
#endif
      }
    }
  }
}

static void gt_sain_singleSinduction1(const GtSaininfo *saininfo,
                                      long *suftab,
                                      long position,
                                      GT_UNUSED unsigned long nonspecialentries,
                                      unsigned long idx)
{
  unsigned long nextcc = gt_sain_seq_getchar(saininfo->sainseq,
                                             (unsigned long) position);

  if (nextcc < saininfo->sainseq->numofchars)
  {
    unsigned long cc,
                  putidx = --saininfo->sainseq->bucketfillptr[nextcc];

    gt_assert(putidx < nonspecialentries && position > 0 &&
              (idx == ULONG_MAX || putidx < idx));
    position--;
    cc = gt_sain_seq_getchar(saininfo->sainseq,(unsigned long) position);
    suftab[putidx] = (cc > nextcc) ? ~(position+1) : position;
#ifdef SAINSHOWSTATE
    printf("end S-induce: suftab[%lu]=%ld\n",putidx,suftab[putidx]);
#endif
  }
  if (idx != ULONG_MAX)
  {
    suftab[idx] = 0;
  }
}

static void gt_sain_induceStypes1fromspecialranges(
                                   const GtSaininfo *saininfo,
                                   const GtEncseq *encseq,
                                   long *suftab,
                                   GT_UNUSED unsigned long nonspecialentries)
{
  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;
    sri = gt_specialrangeiterator_new(encseq,false);
    while (gt_specialrangeiterator_next(sri,&range))
    {
      if (range.start > 1UL)
      {
        gt_sain_singleSinduction1(saininfo,
                                  suftab,
                                  (long) (range.start - 1),
                                  nonspecialentries,
                                  ULONG_MAX);
      }
    }
    gt_specialrangeiterator_delete(sri);
  }
}

static void gt_sain_induceStypesuffixes1(const GtSaininfo *saininfo,
                                         long *suftab,
                                         unsigned long nonspecialentries)
{
  unsigned long idx;

  gt_sain_singleSinduction1(saininfo,
                            suftab,
                            (long) (saininfo->sainseq->totallength-1),
                            nonspecialentries,
                            ULONG_MAX);
  if (saininfo->sainseq->seqtype == GT_SAIN_ENCSEQ)
  {
    gt_sain_induceStypes1fromspecialranges(saininfo,
                                           saininfo->sainseq->seq.encseq,
                                           suftab,
                                           nonspecialentries);
  }
  if (nonspecialentries == 0)
  {
    return;
  }
  for (idx = nonspecialentries - 1; /* Nothing */; idx--)
  {
    long position = suftab[idx];

    if (position > 0)
    {
      gt_sain_singleSinduction1(saininfo,
                                suftab,
                                position,
                                nonspecialentries,
                                idx);
    }
    if (idx == 0)
    {
      break;
    }
  }
}

static void gt_sain_moveSstar2front(const GtSaininfo *saininfo,
                                    long *suftab,
                                    unsigned long nonspecialentries)
{
  unsigned long idx, writeidx = 0;
  long position;

  for (idx = 0; (position = suftab[idx]) < 0; idx++)
  {
    suftab[idx] = ~position;
    gt_assert ((idx + 1) < nonspecialentries);
  }
  if (idx < saininfo->countSstartype)
  {
    for (writeidx = idx, idx++; /* Nothing */; idx++)
    {
      gt_assert (idx < nonspecialentries);
      if ((position = suftab[idx]) < 0)
      {
        gt_assert(writeidx < idx);
        suftab[writeidx++] = ~position;
        suftab[idx] = 0;
        if (writeidx == saininfo->countSstartype)
        {
          break;
        }
      }
    }
  } else
  {
    writeidx = idx;
  }
  gt_assert(writeidx == saininfo->countSstartype);
}

static void gt_sain_singleSinduction2(const GtSaininfo *saininfo,
                                      long *suftab,
                                      long position,
                                      GT_UNUSED unsigned long nonspecialentries,
                                      GT_UNUSED unsigned long idx)
{
  unsigned long nextcc;

  position--;
  nextcc = gt_sain_seq_getchar(saininfo->sainseq,(unsigned long) position);
  if (nextcc < saininfo->sainseq->numofchars)
  {
    unsigned long putidx = --saininfo->sainseq->bucketfillptr[nextcc];

    gt_assert(putidx < nonspecialentries &&
              (idx == ULONG_MAX || idx >= putidx));
    if (position == 0)
    {
      suftab[putidx] = ~0L;
    } else
    {
      unsigned long cc = gt_sain_seq_getchar(saininfo->sainseq,
                                             (unsigned long) (position-1));
      suftab[putidx] = (cc > nextcc) ? ~position : position;
    }
#ifdef SAINSHOWSTATE
    printf("end S-induce: suftab[%lu]=%ld\n",putidx,suftab[putidx]);
#endif
  }
}

static void gt_sain_induceStypes2fromspecialranges(
                                   GT_UNUSED const GtSaininfo *saininfo,
                                   const GtEncseq *encseq,
                                   long *suftab,
                                   GT_UNUSED unsigned long nonspecialentries)
{
  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;
    sri = gt_specialrangeiterator_new(encseq,false);
    while (gt_specialrangeiterator_next(sri,&range))
    {
      if (range.start > 0)
      {
        gt_sain_singleSinduction2(saininfo,
                                  suftab,
                                  (long) range.start,
                                  nonspecialentries,
                                  ULONG_MAX);
      }
    }
    gt_specialrangeiterator_delete(sri);
  }
}

static void gt_sain_induceStypesuffixes2(const GtSaininfo *saininfo,
                                         long *suftab,
                                         unsigned long nonspecialentries)
{
  unsigned long idx;

  gt_sain_singleSinduction2(saininfo,
                            suftab,
                            (long) saininfo->sainseq->totallength,
                            nonspecialentries,
                            ULONG_MAX);
  if (saininfo->sainseq->seqtype == GT_SAIN_ENCSEQ)
  {
    gt_sain_induceStypes2fromspecialranges(saininfo,
                                           saininfo->sainseq->seq.encseq,
                                           suftab,
                                           nonspecialentries);
  }
  if (nonspecialentries == 0)
  {
    return;
  }
  for (idx = nonspecialentries - 1; /* Nothing */; idx--)
  {
    long position = suftab[idx];

    if (position > 0)
    {
      gt_sain_singleSinduction2(saininfo,
                                suftab,
                                position,
                                nonspecialentries,
                                idx);
    } else
    {
      suftab[idx] = ~position;
    }
    if (idx == 0)
    {
      break;
    }
  }
}

static int gt_sain_compare_Sstarstrings(const GtSainseq *sainseq,
                                        unsigned long start1,
                                        unsigned long start2,
                                        unsigned long len)
{
  unsigned long end1 = start1 + len;

  gt_assert(start1 <= sainseq->totallength &&
            start2 <= sainseq->totallength &&
            start1 != start2);

  while (start1 < end1)
  {
    unsigned long cc1, cc2;

    if (start1 == sainseq->totallength)
    {
      gt_assert(start1 > start2);
      return 1;
    }
    if (start2 == sainseq->totallength)
    {
      gt_assert(start1 < start2);
      return -1;
    }
    cc1 = gt_sain_seq_getchar(sainseq,start1);
    cc2 = gt_sain_seq_getchar(sainseq,start2);
    if (cc1 < cc2)
    {
      return -1;
    }
    if (cc1 > cc2)
    {
      return 1;
    }
    start1++;
    start2++;
  }
  return 0;
}

/* The following function has been used previously, but is
   no used. As we maybe use it later, we keep it for now. */

#ifndef NDEBUG
int gt_sain_compare_substrings_new(const GtSaininfo *saininfo,
                                      bool withtype,
                                      unsigned long *lenptr,
                                      unsigned long start1,
                                      unsigned long start2)
{
  unsigned long len = 0;
  bool firstcmp = true, previousisS = true;

  gt_assert(start1 <= saininfo->sainseq->totallength &&
            start2 <= saininfo->sainseq->totallength &&
            start1 != start2);

  for (len=0; /* Nothing */; len++)
  {
    unsigned long cc1, cc2;

    if (start1 == saininfo->sainseq->totallength)
    {
      gt_assert(start1 > start2);
      return 1;
    }
    if (start2 == saininfo->sainseq->totallength)
    {
      gt_assert(start1 < start2);
      return -1;
    }
    cc1 = gt_sain_seq_getchar(saininfo->sainseq,start1);
    cc2 = gt_sain_seq_getchar(saininfo->sainseq,start2);
    if (cc1 < cc2)
    {
      return -1;
    }
    if (cc1 > cc2)
    {
      return 1;
    }
    if (withtype)
    {
      if (gt_sain_info_isStype(saininfo,start1))
      {
        if (gt_sain_info_isStype(saininfo,start2))
        {
          if (!firstcmp && !previousisS)
          {
            if (lenptr != NULL)
            {
              *lenptr = len+1;
            }
            return 0; /* previous is L => Sstar in both stop with equality */
          }
        } else
        {
          /* S > L */
          return 1;
        }
      } else
      {
        if (gt_sain_info_isStype(saininfo,start2))
        {
          /* L < S */
          return -1;
        } else
        {
          /* L == L */
          previousisS = false;
        }
      }
      firstcmp = false;
    }
    start1++;
    start2++;
  }
}
#endif

static int gt_sain_compare_suffixes(const GtSainseq *sainseq,
                                    unsigned long start1,
                                    unsigned long start2)
{
  gt_assert(start1 <= sainseq->totallength &&
            start2 <= sainseq->totallength &&
            start1 != start2);

  while (true)
  {
    unsigned long cc1, cc2;

    if (start1 == sainseq->totallength)
    {
      gt_assert(start1 > start2);
      return 1;
    }
    if (start2 == sainseq->totallength)
    {
      gt_assert(start1 < start2);
      return -1;
    }
    cc1 = gt_sain_seq_getchar(sainseq,start1);
    cc2 = gt_sain_seq_getchar(sainseq,start2);
    if (cc1 < cc2)
    {
      return -1;
    }
    if (cc1 > cc2)
    {
      return 1;
    }
    start1++;
    start2++;
  }
}

static void gt_sain_setundefined(unsigned long *suftab,
                                 unsigned long start, unsigned long end)
{
  unsigned long idx;

  for (idx = start; idx <= end; idx++)
  {
    suftab[idx] = 0;
  }
}

static void gt_sain_assignSstarlength(const GtSaininfo *saininfo,
                                      unsigned long *lentab)
{
  bool nextisStype = true;
  unsigned long position,
                nextSstartypepos = saininfo->sainseq->totallength,
                nextcc = GT_UNIQUEINT(saininfo->sainseq->totallength);

  for (position = saininfo->sainseq->totallength-1; /* Nothing */; position--)
  {
    bool currentisStype;
    unsigned long currentcc;

    currentcc = gt_sain_seq_getchar(saininfo->sainseq,position);
    if (currentcc < nextcc || (currentcc == nextcc && nextisStype))
    {
      currentisStype = true;
    } else
    {
      currentisStype = false;
    }
    if (!currentisStype && nextisStype)
    {
      unsigned long currentlen;

      gt_assert(position < nextSstartypepos);
      currentlen = nextSstartypepos - position;
      lentab[GT_DIV2(position+1)] = currentlen;
      nextSstartypepos = position + 1;
    }
    nextisStype = currentisStype;
    nextcc = currentcc;
    if (position == 0)
    {
      break;
    }
  }
}

static unsigned long gt_sain_assignSstarnames(const GtSaininfo *saininfo,
                                              unsigned long *suftab,
                                              GT_UNUSED unsigned long
                                              availableentries)
{
  unsigned long idx, previouspos, previouslen, currentname = 1UL;

  previouspos = suftab[0];
  previouslen = suftab[saininfo->countSstartype + GT_DIV2(previouspos)];
  suftab[saininfo->countSstartype + GT_DIV2(previouspos)] = currentname;
  for (idx = 1UL; idx < saininfo->countSstartype; idx++)
  {
    int cmp;
    unsigned long currentlen = 0, position = suftab[idx];

    currentlen = suftab[saininfo->countSstartype + GT_DIV2(position)];
    if (previouslen == currentlen)
    {
      cmp = gt_sain_compare_Sstarstrings(saininfo->sainseq,previouspos,position,
                                         currentlen);
      gt_assert(cmp != 1);
    } else
    {
      cmp = -1;
    }
    if (cmp == -1)
    {
      currentname++;
    }
    /* write the names in order of positions. As the positions of
       the Sstar suffixes differ by at least 2, the used address
       is unique */
    previouslen = currentlen;
    gt_assert(saininfo->countSstartype + GT_DIV2(position) < availableentries);
    suftab[saininfo->countSstartype + GT_DIV2(position)] = currentname;
    previouspos = position;
  }
  return currentname;
}

static void gt_sain_movenames2front(const GtSaininfo *saininfo,
                                    unsigned long *suftab,
                                    GT_UNUSED unsigned long availableentries)
{
  unsigned long ridx, widx;
  const unsigned long maxridx = saininfo->countSstartype +
                                GT_DIV2(saininfo->sainseq->totallength);

  for (ridx = widx = saininfo->countSstartype; ridx <= maxridx; ridx++)
  {
    if (suftab[ridx] != 0)
    {
      if (widx < ridx)
      {
        gt_assert(widx < availableentries);
        suftab[widx++] = suftab[ridx] - 1UL; /* As we have used names with
                                                offset 1 to distinguish them
                                                from the undefined values
                                                signified by 0 */
      } else
      {
        gt_assert(widx == ridx);
        suftab[widx]--;
        widx++;
      }
    }
  }
  gt_assert(widx == GT_MULT2(saininfo->countSstartype));
}

static void gt_sain_checkorder(const GtSainseq *sainseq,
                               const unsigned long *suftab,
                               unsigned long start,
                               unsigned long end)
{
  unsigned long idx;

  for (idx = start+1; idx <= end; idx++)
  {
    int cmp
      = gt_sain_compare_suffixes(sainseq,suftab[idx-1],suftab[idx]);

    if (cmp != -1)
    {
      fprintf(stderr,"%s: check interval [%lu,%lu] at idx=%lu: suffix "
                     "%lu >= %lu\n",__func__,start,end,idx,suftab[idx-1],
                     suftab[idx]);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

static void gt_sain_insertsortedSstarsuffixes(const GtSaininfo *saininfo,
                                              unsigned long *suftab,
                                              GT_UNUSED unsigned long
                                                        nonspecialsuffixes)
{
  unsigned long idx,
                *fillptr = saininfo->sainseq->bucketfillptr;

  if (saininfo->countSstartype == 0)
  {
    return;
  }
  for (idx = saininfo->countSstartype - 1; /* Nothing */; idx--)
  {
    unsigned long position = suftab[idx], putidx, cc;

    cc = gt_sain_seq_getchar_nospecial(saininfo->sainseq,position);
    gt_assert(fillptr[cc] > 0);
    putidx = --fillptr[cc];
    gt_assert(idx <= putidx);
    if (idx < putidx)
    {
      gt_assert(position > 0 && putidx < nonspecialsuffixes);
      suftab[putidx] = position;
      suftab[idx] = 0;
#ifdef SAINSHOWSTATE
      printf("insertsorted: suftab[%lu]=%lu\n",putidx,position);
      printf("insertsorted: suftab[%lu]=undef\n",idx);
#endif
    }
    if (idx == 0)
    {
      break;
    }
  }
}

static void gt_sain_expandorder2original(const GtSaininfo *saininfo,
                                         unsigned long *suftab)
{
  unsigned long idx = saininfo->countSstartype - 1, position,
                nextcc = GT_UNIQUEINT(saininfo->sainseq->totallength),
                *sstarsuffixes = suftab + saininfo->countSstartype;
  bool nextisStype = true;

  for (position = saininfo->sainseq->totallength-1; /* Nothing */; position--)
  {
    bool currentisStype;
    unsigned long currentcc;

    currentcc = gt_sain_seq_getchar(saininfo->sainseq,position);
    if (currentcc < nextcc || (currentcc == nextcc && nextisStype))
    {
      currentisStype = true;
    } else
    {
      currentisStype = false;
    }
    if (!currentisStype && nextisStype)
    {
      sstarsuffixes[idx--] = position+1;
    }
    nextisStype = currentisStype;
    nextcc = currentcc;
    if (position == 0)
    {
      break;
    }
  }
  for (idx = 0; idx < saininfo->countSstartype; idx++)
  {
    suftab[idx] = sstarsuffixes[suftab[idx]];
  }
}

static void gt_sain_filltailsuffixes(unsigned long *suftabtail,
                                     const GtEncseq *encseq)
{
  unsigned long specialcharacters = gt_encseq_specialcharacters(encseq),
                totallength = gt_encseq_total_length(encseq);

  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;
    unsigned long countspecial = 0;

    sri = gt_specialrangeiterator_new(encseq,true);
    while (gt_specialrangeiterator_next(sri,&range))
    {
      unsigned long idx;

      for (idx = range.start; idx < range.end; idx++)
      {
        gt_assert(countspecial < specialcharacters && idx < totallength);
        suftabtail[countspecial++] = idx;
      }
    }
    gt_assert(countspecial == specialcharacters);
    gt_specialrangeiterator_delete(sri);
  }
  suftabtail[specialcharacters] = totallength;
}

static void gt_sain_rec_sortsuffixes(unsigned int level,
                                     GtSainseq *sainseq,
                                     unsigned long *suftab,
                                     unsigned long nonspecialentries,
                                     unsigned long availableentries,
                                     unsigned long suftabentries,
                                     bool intermediatecheck,
                                     bool finalcheck,
                                     bool verbose,
                                     GtTimer *timer)
{
  GtSaininfo *saininfo;

  printf("level %u: sort sequence of length %lu over %lu symbols (%.2f)\n",
          level,sainseq->totallength,sainseq->numofchars,
          (double) sainseq->numofchars/sainseq->totallength);
  GT_SAIN_SHOWTIMER("insert Sstar suffixes");
  saininfo = gt_sain_info_new(sainseq,suftab,nonspecialentries);
  if (verbose)
  {
    gt_sain_info_show(saininfo);
  }
  if (saininfo->countSstartype > 0)
  {
    unsigned long startoccupied, numberofnames;

    gt_sain_startbuckets(saininfo->sainseq);
    GT_SAIN_SHOWTIMER("induce L suffixes");
    gt_sain_induceLtypesuffixes1(saininfo,(long *) suftab,nonspecialentries);
    gt_sain_endbuckets(saininfo->sainseq);
    GT_SAIN_SHOWTIMER("induce S suffixes");
    gt_sain_induceStypesuffixes1(saininfo,(long *) suftab,nonspecialentries);
    GT_SAIN_SHOWTIMER("moverStar2front");
    gt_sain_moveSstar2front(saininfo,(long *) suftab,nonspecialentries);
    gt_assert(availableentries > 0);
    gt_sain_setundefined(suftab,saininfo->countSstartype,availableentries-1);
    GT_SAIN_SHOWTIMER("assignSstarlength");
    gt_sain_assignSstarlength(saininfo,suftab + saininfo->countSstartype);
    GT_SAIN_SHOWTIMER("assignSstarnames");
    numberofnames = gt_sain_assignSstarnames(saininfo,suftab,availableentries);
    gt_assert(numberofnames <= saininfo->countSstartype);
    if (numberofnames < saininfo->countSstartype)
    {
    /* Now the name sequence is in the range from
       saininfo->countSstartype .. 2 * saininfo->countSstartype - 1 */
      unsigned long *subseq = suftab + saininfo->countSstartype;
      GtSainseq *sainseq_rec;

      GT_SAIN_SHOWTIMER("movenames2front");
      gt_sain_movenames2front(saininfo,suftab,availableentries);
      gt_sain_setundefined(suftab,0,saininfo->countSstartype-1);
      sainseq_rec = gt_sain_seq_new_from_array(subseq,
                                               saininfo->countSstartype,
                                               numberofnames,
                                               suftab,
                                               suftabentries);
      gt_sain_rec_sortsuffixes(level+1,sainseq_rec,suftab,
                               saininfo->countSstartype,
                               saininfo->countSstartype,
                               suftabentries,
                               intermediatecheck,
                               finalcheck,
                               verbose,
                               timer);
      startoccupied = gt_sain_seq_delete(sainseq_rec);
      gt_sain_setundefined(suftab,startoccupied,suftabentries-1);
      GT_SAIN_SHOWTIMER("expandorder2original");
      gt_sain_expandorder2original(saininfo,suftab);
    }
  }
  if (intermediatecheck && saininfo->countSstartype > 0)
  {
    gt_sain_checkorder(saininfo->sainseq,suftab,0,saininfo->countSstartype-1);
  }
  if (saininfo->sainseq->seqtype == GT_SAIN_INTSEQ)
  {
    gt_sain_determinebucketsize(saininfo->sainseq);
  }
  gt_sain_setundefined(suftab,saininfo->countSstartype,availableentries-1);
  gt_sain_endbuckets(saininfo->sainseq);
  gt_assert(saininfo->sainseq->bucketfillptr[saininfo->sainseq->numofchars-1]
            == nonspecialentries);
  GT_SAIN_SHOWTIMER("insert sorted Sstar suffixes");
  gt_sain_insertsortedSstarsuffixes (saininfo, suftab, nonspecialentries);
  gt_sain_startbuckets(saininfo->sainseq);
  GT_SAIN_SHOWTIMER("induce L suffixes");
  gt_sain_induceLtypesuffixes2(saininfo,(long *) suftab,nonspecialentries);
  gt_sain_endbuckets(saininfo->sainseq);
  GT_SAIN_SHOWTIMER("induce S suffixes");
  gt_sain_induceStypesuffixes2(saininfo,(long *) suftab,nonspecialentries);

  if (nonspecialentries > 0)
  {
    if (intermediatecheck)
    {
      gt_sain_checkorder(saininfo->sainseq,suftab,0,nonspecialentries-1);
    }
    if (saininfo->sainseq->seqtype == GT_SAIN_ENCSEQ && finalcheck)
    {
      GT_SAIN_SHOWTIMER("fill tail suffixes");
      gt_sain_filltailsuffixes(suftab + nonspecialentries,
                               saininfo->sainseq->seq.encseq);
      GT_SAIN_SHOWTIMER("check suffix order");
      gt_suftab_lightweightcheck(saininfo->sainseq->seq.encseq,
                                 GT_READMODE_FORWARD,
                                 saininfo->sainseq->totallength,
                                 suftab,
                                 NULL);
    }
  }
  gt_sain_info_delete(saininfo);
}

void gt_sain_encseq_sortsuffixesnew(const GtEncseq *encseq,
                                    bool intermediatecheck,
                                    bool finalcheck,bool verbose,
                                    GtTimer *timer)
{
  unsigned long nonspecialentries, suftabentries, totallength, *suftab;
  GtSainseq *sainseq;

  totallength = gt_encseq_total_length(encseq);
  nonspecialentries = totallength - gt_encseq_specialcharacters(encseq);
  suftabentries = totallength+1;
  suftab = gt_malloc(sizeof (*suftab) * suftabentries);
  gt_sain_setundefined(suftab,0,suftabentries-1);
  sainseq = gt_sain_seq_new_from_encseq(encseq);
  gt_sain_rec_sortsuffixes(0,sainseq,suftab,nonspecialentries,
                           suftabentries,
                           suftabentries,intermediatecheck,finalcheck,verbose,
                           timer);
  (void) gt_sain_seq_delete(sainseq);
  gt_free(suftab);
}

void gt_sain_plain_sortsuffixesnew(const GtUchar *plainseq,
                                   unsigned long len, bool intermediatecheck,
                                   bool verbose,
                                   GtTimer *timer)
{
  unsigned long suftabentries, *suftab;
  GtSainseq *sainseq;

  suftabentries = len+1;
  suftab = gt_malloc(sizeof (*suftab) * suftabentries);
  gt_sain_setundefined(suftab,0,suftabentries - 1);
  sainseq = gt_sain_seq_new_from_plainseq(plainseq,len);
  gt_sain_rec_sortsuffixes(0,sainseq,suftab,sainseq->totallength,
                           suftabentries,
                           suftabentries,intermediatecheck,false,verbose,
                           timer);
  printf("countcharaccess=%lu (%.2f)\n",gt_sain_countcharaccess,
          (double) gt_sain_countcharaccess/sainseq->totallength);
  (void) gt_sain_seq_delete(sainseq);
  gt_free(suftab);
}
