/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#include <assert.h>
#include "libgtcore/cstr.h"
#include "libgtcore/fa.h"
#include "libgtcore/io.h"
#include "libgtcore/ma.h"
#include "libgtcore/xansi.h"

struct IO {
  FILE *fp;
  char *path;
  unsigned long line_number;
  bool line_start;
};

IO* io_new(const char *path, const char *mode, Env *env)
{
  IO *io;
  assert(path && mode);
  assert(!strcmp(mode, "r")); /* XXX: only the read mode has been implemented */
  io = ma_malloc(sizeof (IO));
  io->fp = fa_xfopen(path, mode);
  io->path = cstr_dup(path, env);
  io->line_number = 1;
  io->line_start = true;
  return io;
}

int io_get_char(IO *io, char *c)
{
  int cc;
  assert(io && c);
  cc = xfgetc(io->fp);
  if (cc == '\n') {
    io->line_number++;
    io->line_start = true;
  }
  else
    io->line_start = false;
  if (cc == EOF)
    return -1; /* no character left */
  *c = cc;
  return 0;
}

void io_unget_char(IO *io, char c)
{
  assert(io);
  xungetc(c, io->fp);
}

bool io_line_start(const IO *io)
{
  assert(io);
  return io->line_start;
}

unsigned long io_get_line_number(const IO *io)
{
  assert(io);
  return io->line_number;
}

const char* io_get_filename(const IO *io)
{
  assert(io && io->path);
  return io->path;
}

void io_delete(IO *io, Env *env)
{
  if (!io) return;
  fa_xfclose(io->fp);
  ma_free(io->path);
  ma_free(io);
}
