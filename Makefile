#
# Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
# Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg 
# See LICENSE file or http://genometools.org/license.html for license details
#

CC:=gcc
CXX:=g++
INCLUDEOPT:= -I$(CURDIR)/src -I$(CURDIR)/obj \
             -I$(CURDIR)/src/external/lua-5.1.2/src \
             -I$(CURDIR)/src/external/expat-2.0.1/lib\
             -I$(CURDIR)/src/external/bzip2-1.0.4\
             -I$(CURDIR)/src/external/agg-2.4/include\
             -I$(CURDIR)/src/external/libpng-1.2.18\
             -I/usr/include/cairo\
             -I/usr/local/include/cairo

CFLAGS:=
LDFLAGS:=
CXXFLAGS:=
GT_CFLAGS:= -g -Wall -Werror -pipe -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 \
            $(INCLUDEOPT)
# expat needs -DHAVE_MEMMOVE
# lua needs -DLUA_USE_POSIX
# rnv needs -DUNISTD_H="<unistd.h>" -DEXPAT_H="<expat.h>" -DRNV_VERSION="\"1.7.8\""
EXT_FLAGS:= -DHAVE_MEMMOVE -DLUA_USE_POSIX -DUNISTD_H="<unistd.h>" \
            -DEXPAT_H="<expat.h>" -DRNV_VERSION="\"1.7.8\""
GT_CXXFLAGS:= -g -pipe $(INCLUDEOPT)
STEST_FLAGS:=
GT_LDFLAGS:=-L/usr/local/lib -L/usr/X11R6/lib
LDLIBS:=-lm -lz

# try to set RANLIB automatically
SYSTEM:=$(shell uname -s)
ifeq ($(SYSTEM),Darwin)
  RANLIB:=ranlib
  GT_CFLAGS+=-DHAVE_LLABS
endif
ifeq ($(SYSTEM),Linux)
  GT_CLAGS+=-D_ISOC99_SOURCE -D_XOPEN_SOURCE=600 -DHAVE_LLABS
endif

# the default GenomeTools libraries which are build
GTLIBS:=lib/libgtext.a\
        lib/libgtcore.a\
        lib/libgtmatch.a\
        lib/libbz2.a

# the core GenomeTools library (no other dependencies)
LIBGTCORE_SRC:=$(wildcard src/libgtcore/*.c)
LIBGTCORE_OBJ:=$(LIBGTCORE_SRC:%.c=obj/%.o)
LIBGTCORE_DEP:=$(LIBGTCORE_SRC:%.c=obj/%.d)

# the extended GenomeTools library (e.g., depends on Lua)
LIBGTEXT_C_SRC:=$(wildcard src/libgtext/*.c)
LIBGTEXT_C_OBJ:=$(LIBGTEXT_C_SRC:%.c=obj/%.o)
LIBGTEXT_C_DEP:=$(LIBGTEXT_C_SRC:%.c=obj/%.d)
LIBGTEXT_CXX_SRC:=$(wildcard src/libgtext/*.cxx)
LIBGTEXT_CXX_OBJ:=$(LIBGTEXT_CXX_SRC:%.cxx=obj/%.o)
LIBGTEXT_CXX_DEP:=$(LIBGTEXT_CXX_SRC:%.cxx=obj/%.d)

# the GenomeTools matching library
LIBGTMATCH_SRC:=$(wildcard src/libgtmatch/*.c)
LIBGTMATCH_OBJ:=$(LIBGTMATCH_SRC:%.c=obj/%.o)
LIBGTMATCH_DEP:=$(LIBGTMATCH_SRC:%.c=obj/%.d)

# the GenomeTools view library
LIBGTVIEW_C_SRC:=$(wildcard src/libgtview/*.c)
LIBGTVIEW_C_OBJ:=$(LIBGTVIEW_C_SRC:%.c=obj/%.o)
LIBGTVIEW_C_DEP:=$(LIBGTVIEW_C_SRC:%.c=obj/%.d)

TOOLS_SRC:=$(wildcard src/tools/*.c)
TOOLS_OBJ:=$(TOOLS_SRC:%.c=obj/%.o)
TOOLS_DEP:=$(TOOLS_SRC:%.c=obj/%.d)

LIBAGG_SRC:=$(wildcard src/external/agg-2.4/src/*.cpp src/external/agg-2.4/src/ctrl/*.cpp)
LIBAGG_OBJ:=$(LIBAGG_SRC:%.cpp=obj/%.o)
LIBAGG_DEP:=$(LIBAGG_SRC:%.cpp=obj/%.d)

EXPAT_DIR:=src/external/expat-2.0.1/lib
LIBEXPAT_SRC:=$(EXPAT_DIR)/xmlparse.c $(EXPAT_DIR)/xmlrole.c \
              $(EXPAT_DIR)/xmltok.c
LIBEXPAT_OBJ:=$(LIBEXPAT_SRC:%.c=obj/%.o)
LIBEXPAT_DEP:=$(LIBEXPAT_SRC:%.c=obj/%.d)

LUA_DIR:=src/external/lua-5.1.2/src
LIBLUA_SRC=$(LUA_DIR)/lapi.c $(LUA_DIR)/lcode.c $(LUA_DIR)/ldebug.c \
           $(LUA_DIR)/ldo.c $(LUA_DIR)/ldump.c $(LUA_DIR)/lfunc.c \
           $(LUA_DIR)/lgc.c $(LUA_DIR)/llex.c $(LUA_DIR)/lmem.c \
           $(LUA_DIR)/lobject.c $(LUA_DIR)/lopcodes.c $(LUA_DIR)/lparser.c \
           $(LUA_DIR)/lstate.c $(LUA_DIR)/lstring.c $(LUA_DIR)/ltable.c \
           $(LUA_DIR)/ltm.c $(LUA_DIR)/lundump.c $(LUA_DIR)/lvm.c \
           $(LUA_DIR)/lzio.c $(LUA_DIR)/lauxlib.c $(LUA_DIR)/lbaselib.c \
           $(LUA_DIR)/ldblib.c $(LUA_DIR)/liolib.c $(LUA_DIR)/lmathlib.c \
           $(LUA_DIR)/loslib.c $(LUA_DIR)/ltablib.c $(LUA_DIR)/lstrlib.c \
           $(LUA_DIR)/loadlib.c $(LUA_DIR)/linit.c
LIBLUA_OBJ:=$(LIBLUA_SRC:%.c=obj/%.o)
LIBLUA_DEP:=$(LIBLUA_SRC:%.c=obj/%.d)

PNG_DIR:=external/libpng-1.2.18
LIBPNG_SRC:=$(PNG_DIR)/png.o $(PNG_DIR)/pngset.o $(PNG_DIR)/pngget.o \
            $(PNG_DIR)/pngrutil.o $(PNG_DIR)/pngtrans.o $(PNG_DIR)/pngwutil.o \
            $(PNG_DIR)/pngread.o $(PNG_DIR)/pngrio.o $(PNG_DIR)/pngwio.o \
            $(PNG_DIR)/pngwrite.o $(PNG_DIR)/pngrtran.o $(PNG_DIR)/pngwtran.o \
            $(PNG_DIR)/pngmem.o $(PNG_DIR)/pngerror.o $(PNG_DIR)/pngpread.o
LIBPNG_OBJ:=$(LIBPNG_SRC:%.c=obj/%.o)
LIBPNG_DEP:=$(LIBPNG_SRC:%.c=obj/%.d)

RNV_DIR:=src/external/rnv-1.7.8
LIBRNV_SRC:=$(RNV_DIR)/rn.c $(RNV_DIR)/rnc.c $(RNV_DIR)/rnd.c $(RNV_DIR)/rnl.c \
            $(RNV_DIR)/rnv.c $(RNV_DIR)/rnx.c $(RNV_DIR)/drv.c \
            $(RNV_DIR)/ary.c $(RNV_DIR)/xsd.c $(RNV_DIR)/xsd_tm.c \
            $(RNV_DIR)/dxl.c $(RNV_DIR)/dsl.c $(RNV_DIR)/sc.c $(RNV_DIR)/u.c \
            $(RNV_DIR)/ht.c $(RNV_DIR)/er.c $(RNV_DIR)/xmlc.c $(RNV_DIR)/s.c \
            $(RNV_DIR)/m.c $(RNV_DIR)/rx.c
LIBRNV_OBJ:=$(LIBRNV_SRC:%.c=obj/%.o)
LIBRNV_DEP:=$(LIBRNV_SRC:%.c=obj/%.d)

BZ2_DIR:=src/external/bzip2-1.0.4
LIBBZ2_SRC:=$(BZ2_DIR)/blocksort.c $(BZ2_DIR)/huffman.c $(BZ2_DIR)/crctable.c \
            $(BZ2_DIR)/randtable.c $(BZ2_DIR)/compress.c \
            $(BZ2_DIR)/decompress.c $(BZ2_DIR)/bzlib.c
LIBBZ2_OBJ:=$(LIBBZ2_SRC:%.c=obj/%.o)
LIBBZ2_DEP:=$(LIBBZ2_SRC:%.c=obj/%.o)

SERVER=gordon@genometools.org
WWWBASEDIR=/var/www/servers

# process arguments
ifneq ($(opt),no)
  GT_CFLAGS += -Os
  GT_CXXFLAGS += -Os
endif

ifeq ($(assert),no)
  GT_CFLAGS += -DNDEBUG
  GT_CXXFLAGS += -DNDEBUG
endif

ifeq ($(static),yes)
  GT_LDFLAGS += -static
endif

ifeq ($(libgtview),yes)
  GTLIBS := $(GTLIBS) lib/libgtview.a
  GT_CFLAGS += -DLIBGTVIEW
  LDLIBS += -lcairo
  STEST_FLAGS += -libgtview
endif

# set prefix for install target
prefix ?= /usr/local

all: $(GTLIBS) bin/skproto bin/gt bin/rnv

lib/libexpat.a: $(LIBEXPAT_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBEXPAT_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libagg.a: $(LIBAGG_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBAGG_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libbz2.a: $(LIBBZ2_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBBZ2_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libgtcore.a: obj/gt_build.h obj/gt_cc.h obj/gt_cflags.h obj/gt_version.h \
                 $(LIBGTCORE_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBGTCORE_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libgtext.a: $(LIBGTEXT_C_OBJ) $(LIBGTEXT_CXX_OBJ) $(LIBLUA_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBGTEXT_C_OBJ) $(LIBGTEXT_CXX_OBJ) $(LIBLUA_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libgtmatch.a: $(LIBGTMATCH_OBJ)
	@echo "[link $@]"
	@ar ru $@ $(LIBGTMATCH_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libgtview.a: $(LIBGTVIEW_C_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBGTVIEW_C_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/libpng.a: $(LIBPNG_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBPNG_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/librnv.a: $(LIBRNV_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@ar ru $@ $(LIBRNV_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

bin/skproto: obj/src/skproto.o obj/src/tools/gt_skproto.o lib/libgtcore.a lib/libbz2.a
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CXX) $(LDFLAGS) $(GT_LDFLAGS) $^ $(LDLIBS) -o $@

bin/gt: obj/src/gt.o obj/src/gtr.o obj/src/gtlua.o $(TOOLS_OBJ) $(GTLIBS)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CXX) $(LDFLAGS) $(GT_LDFLAGS) $^ $(LDLIBS) -o $@

bin/rnv: obj/$(RNV_DIR)/xcl.o lib/librnv.a lib/libexpat.a
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CC) $(LDFLAGS) $(GT_LDFLAGS) $^ -o $@

obj/gt_build.h:
	@date +'#define GT_BUILT "%Y-%m-%d %H:%M:%S"' > $@

obj/gt_cc.h:
	@echo '#define GT_CC "'`$(CC) --version | head -n 1`\" > $@

obj/gt_cflags.h:
	@echo '#define GT_CFLAGS "$(CFLAGS) $(GT_CFLAGS)"' > $@

obj/gt_version.h: VERSION
	@echo '#define GT_VERSION "'`cat VERSION`\" > $@

Doxyfile: Doxyfile.in VERSION
	@echo '[rebuild $@]'
	@sed -e "s/@VERSION@/`cat VERSION`/g" <$< >$@

src/libgtcore/bitpackstringop8.c: src/libgtcore/bitpackstringop.template
	@echo '[rebuild $@]'
	@scripts/template2c.pl 8 $^

src/libgtcore/checkbitpackstring8.c: src/libgtcore/checkbitpackstring.template
	@echo '[rebuild $@]'
	@scripts/template2c.pl 8 $^

src/libgtcore/bitpackstringop16.c: src/libgtcore/bitpackstringop.template
	@echo '[rebuild $@]'
	@scripts/template2c.pl 16 $^

src/libgtcore/checkbitpackstring16.c: src/libgtcore/checkbitpackstring.template
	@echo '[rebuild $@]'
	@scripts/template2c.pl 16 $^

src/libgtcore/bitpackstringop32.c: src/libgtcore/bitpackstringop.template
	@echo '[rebuild $@]'
	@scripts/template2c.pl 32 $^

src/libgtcore/checkbitpackstring32.c: src/libgtcore/checkbitpackstring.template
	@echo '[rebuild $@]'
	@scripts/template2c.pl 32 $^

src/libgtcore/bitpackstringop64.c: src/libgtcore/bitpackstringop.template
	@echo '[rebuild $@]'
	@scripts/template2c.pl 64 $^

src/libgtcore/checkbitpackstring64.c: src/libgtcore/checkbitpackstring.template
	@echo '[rebuild $@]'
	@scripts/template2c.pl 64 $^

src/libgtcore/checkbitpackstring-int.c: src/libgtcore/checkbitpackstring.template
	@echo '[rebuild $@]'
	@scripts/template2c.pl '-int' $^

# we create the dependency files on the fly
obj/%.o: %.c
	@echo "[compile $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CC) -c $< -o $@  $(CFLAGS) $(GT_CFLAGS) $(EXT_FLAGS) \
        -MT $@ -MMD -MP -MF $(@:.o=.d)

obj/%.o: %.cxx
	@echo "[compile $@]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CXX) -c $< -o $@  $(CXXFLAGS) $(GT_CXXFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)

# read deps
-include obj/src/gt.d obj/src/gtlua.d obj/src/gtr.d obj/src/skproto.d \
         $(LIBGTCORE_DEP) $(LIBGTEXT_C_DEP) $(LIBGTEXT_CXX_DEP) \
         $(LIBGTMATCH_DEP) $(LIBGTVIEW_C_DEP) $(TOOLS_DEP) $(LIBAGG_DEP) \
         $(LIBEXPAT_DEP) $(LIBLUA_DEP) $(LIBPNG_DEP) $(LIBRNV_DEP) \
         $(LIBBBZ2_DEP)

.SUFFIXES:
.PHONY: dist srcdist release gt install splint test clean cleanup apidoc

dist: all
	tar cvzf gt-`cat VERSION`.tar.gz bin/gt

srcdist:
	git archive --format=tar --prefix=genometools-`cat VERSION`/ HEAD | \
        gzip -9 > genometools-`cat VERSION`.tar.gz

apidoc: Doxyfile
	doxygen

release: apidoc
	git tag "v`cat VERSION`"
	git archive --format=tar --prefix="genometools-`cat VERSION`"/ HEAD \
	>"genometools-`cat VERSION`.tar"
	mkdir -p "genometools-`cat VERSION`/doc"
	rsync -a doc/api "genometools-`cat VERSION`/doc/"
	tar -r -f "genometools-`cat VERSION`.tar" \
		"genometools-`cat VERSION`/doc/api"
	rm -rf "genometools-`cat VERSION`"
	gzip -9 "genometools-`cat VERSION`.tar"
	scp "genometools-`cat VERSION`.tar.gz" $(SERVER):$(WWWBASEDIR)/genometools.org/htdocs/pub
	rsync -rv doc/api/html/ $(SERVER):$(WWWBASEDIR)/genometools.org/htdocs/apidoc
	git push --tags

installwww:
# install genometools.org website
	rsync -rv www/genometools.org/ $(SERVER):$(WWWBASEDIR)/genometools.org

gt: bin/gt

install:
	test -d $(prefix)/bin || mkdir -p $(prefix)/bin
	cp bin/gt $(prefix)/bin
	cp -r gtdata $(prefix)/bin
	test -d $(prefix)/include/libgtcore || mkdir -p $(prefix)/include/libgtcore
	cp src/gtcore.h $(prefix)/include
	cp src/libgtcore/*.h $(prefix)/include/libgtcore
	test -d $(prefix)/include/libgtext || mkdir -p $(prefix)/include/libgtext
	cp src/gtext.h $(prefix)/include
	cp src/libgtext/*.h $(prefix)/include/libgtext
	cp src/gt.h $(prefix)/include
	test -d $(prefix)/lib || mkdir -p $(prefix)/lib
	cp lib/libgtcore.a $(prefix)/lib
ifdef RANLIB
	$(RANLIB) $(prefix)/lib/libgtcore.a
endif
	cp lib/libgtext.a $(prefix)/lib
ifdef RANLIB
	$(RANLIB) $(prefix)/lib/libgtext.a
endif

splint:
	splint -f $(CURDIR)/testdata/Splintoptions $(INCLUDEOPT) $(CURDIR)/src/*.c \
        $(CURDIR)/src/libgtcore/*.c \
        $(CURDIR)/src/libgtext/*.c \
        $(CURDIR)/src/tools/*.c

sgt:${addprefix obj/,${notdir ${subst .c,.splint,${wildcard ${CURDIR}/src/libgtmatch/*.c}}}}


obj/%.splint: ${CURDIR}/src/libgtmatch/%.c
	@echo "splint $<"
	@splint -Isrc -f $(CURDIR)/testdata/SKsplintoptions $<
	@touch $@

obj/%.prepro: ${CURDIR}/src/libgtmatch/%.c
	@echo "[generate $@]"
	${CC} -c $< -o $@ ${CFLAGS} ${GT_CFLAGS} -E -g3
	indent $@

test: all
	bin/gt -test
	cd testsuite && env -i ruby -I. testsuite.rb -testdata $(CURDIR)/testdata -bin $(CURDIR)/bin $(STEST_FLAGS)

clean:
	rm -rf obj
	rm -rf testsuite/stest_testsuite testsuite/stest_stest_tests

cleanup: clean
	rm -rf lib bin Doxyfile doc/api
