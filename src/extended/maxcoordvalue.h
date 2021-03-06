/*
  Copyright (c) 2015 Annika <annika.seidel@studium.uni-hamburg.de>
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
#ifndef MAXCOORDVALUE_H
#define MAXCOORDVALUE_H
#include "core/types_api.h"

typedef struct Gtmaxcoordvalue Gtmaxcoordvalue;

Gtmaxcoordvalue* gt_max_new(void);

void gt_max_delete(Gtmaxcoordvalue *max);

GtWord gt_max_get_value(const Gtmaxcoordvalue *max);

void gt_max_set_start(Gtmaxcoordvalue *max, GtUword starta, GtUword startb);

GtUwordPair gt_max_get_start(const Gtmaxcoordvalue *max);

void gt_max_set_end_with_pair(Gtmaxcoordvalue *max, GtUwordPair end);

GtUwordPair gt_max_get_end(const Gtmaxcoordvalue *max);

/*use this in linear space context*/
void gt_max_coord_update(Gtmaxcoordvalue *max, GtWord value,
                         GtUwordPair start,
                         GtUword enda, GtUword endb);
/*use this in square space context*/
void gt_max_coord_update_without_start (Gtmaxcoordvalue *max, GtWord value,
                                        GtUword enda, GtUword endb);
GtUword gt_max_get_row_length(const Gtmaxcoordvalue *max);

GtUword gt_max_get_col_length(const Gtmaxcoordvalue *max);

bool gt_max_get_length_safe(const Gtmaxcoordvalue *max);

void gt_max_reset(Gtmaxcoordvalue *max);
#endif
