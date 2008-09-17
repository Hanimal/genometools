/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include <math.h>
#include <string.h>
#include "core/bittab.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/range.h"
#include "core/unused_api.h"
#include "annotationsketch/canvas.h"
#include "annotationsketch/canvas_cairo_file.h"
#include "annotationsketch/canvas_rep.h"
#include "annotationsketch/graphics_cairo.h"
#include "annotationsketch/style.h"

#define MARGINS_DEFAULT           10
#define HEADER_SPACE              70

struct GtCanvasCairoFile {
  const GtCanvas parent_instance;
  GtGraphicsOutType type;
};

int gt_canvas_cairo_file_visit_diagram_pre(GtCanvas *canvas, GtDiagram *dia)
{
  double margins;

  assert(canvas && dia);

  if (gt_style_get_num(canvas->pvt->sty, "format", "margins", &margins, NULL))
    canvas->pvt->margins = margins;
  else
    canvas->pvt->margins = MARGINS_DEFAULT;

  if (!gt_style_get_bool(canvas->pvt->sty, "format", "show_track_captions",
                         &canvas->pvt->show_track_captions, NULL))
    canvas->pvt->show_track_captions = true;

  canvas->pvt->viewrange = gt_diagram_get_range(dia);
  if (canvas->pvt->g)
  {
    gt_graphics_delete(canvas->pvt->g);
    canvas->pvt->g = NULL;
  }
  canvas->pvt->g = gt_graphics_cairo_new(((GtCanvasCairoFile*) canvas)->type,
                                  canvas->pvt->width, 1);

  /* calculate scaling factor */
  canvas->pvt->factor = ((double) canvas->pvt->width
                     -(2*canvas->pvt->margins))
                    / gt_range_length(canvas->pvt->viewrange);
  return 0;
}

int gt_canvas_cairo_file_visit_diagram_post(GtCanvas *canvas, GtDiagram *dia)
{
  int had_err = 0;

  assert(canvas && dia);

  /* set initial image-specific values */
  canvas->pvt->y += HEADER_SPACE;
  canvas->pvt->height = gt_canvas_calculate_height(canvas, dia);
  if (canvas->pvt->ii)
    gt_image_info_set_height(canvas->pvt->ii, canvas->pvt->height);
  if (canvas->pvt->g)
  {
    gt_graphics_delete(canvas->pvt->g);
    canvas->pvt->g = NULL;
  }
  canvas->pvt->g = gt_graphics_cairo_new(((GtCanvasCairoFile*) canvas)->type,
                                 canvas->pvt->width, canvas->pvt->height);
  gt_graphics_set_margins(canvas->pvt->g, canvas->pvt->margins, 0);

  /* Add ruler/scale to the image */
  gt_canvas_draw_ruler(canvas);

  return had_err;
}

int gt_canvas_cairo_file_to_file(GtCanvasCairoFile *canvas,
                                 const char *filename, GtError *err)
{
  int had_err = 0;
  GtCanvas *c = (GtCanvas*) canvas;
  gt_error_check(err);
  assert(canvas && filename);
  /* write out result file */
  if (c->pvt->g)
    had_err = gt_graphics_save_to_file(c->pvt->g, filename, err);
  else
  {
    /* XXX: shouldn't this be an assertion? */
    gt_error_set(err, "No graphics has been created yet!");
    had_err = -1;
  }

  return had_err;
}

int gt_canvas_cairo_file_to_stream(GtCanvasCairoFile *canvas, GtStr *stream)
{
  int had_err = 0;
  GtCanvas *c = (GtCanvas*) canvas;
  assert(canvas && stream);

  /* write out result file */
  if (c->pvt->g)
    gt_graphics_save_to_stream(c->pvt->g, stream);

  return had_err;
}

const GtCanvasClass* gt_canvas_cairo_file_class(void)
{
  static const GtCanvasClass canvas_class =
    { sizeof (GtCanvasCairoFile),
      gt_canvas_cairo_file_visit_diagram_pre,
      gt_canvas_cairo_file_visit_diagram_post,
      NULL };
  return &canvas_class;
}

GtCanvas* gt_canvas_cairo_file_new(GtStyle *sty, GtGraphicsOutType type,
                                    unsigned long width, GtImageInfo *ii)
{
  GtCanvas *canvas;
  GtCanvasCairoFile *ccf;
  assert(sty && width > 0);
  canvas = gt_canvas_create(gt_canvas_cairo_file_class());
  canvas->pvt->sty = sty;
  canvas->pvt->ii = ii;
  canvas->pvt->width = width;
  canvas->pvt->bt = NULL;
  canvas->pvt->y = 0.5; /* 0.5 displacement to eliminate fuzzy horizontal lines */
  ccf = canvas_cairo_file_cast(canvas);
  ccf->type = type;
  return canvas;
}
