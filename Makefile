PATHTOCLASS = /usr/local/share/CLASS/class_public


CC = gcc
ARCHIVE = libhmpdf.a
SHARED = libhmpdf.so

all: $(SHARED) $(ARCHIVE)

CFLAGS = -std=gnu99 -fPIC -Wall -Wpedantic -Wno-variadic-macros -DHAVE_INLINE
OPTFLAGS = -O4 -ffast-math
OMPFLAGS = -fopenmp

INCLUDE = -I./include
INCLUDE += -I$(PATHTOCLASS)/include

LINKER = -L$(PATHTOCLASS)
LINKER += -lclass -lgsl -lgslcblas -lm -lfftw3

SRCDIR = ./src
OBJDIR = ./obj
OUTDIR = ./lib
SODIR  = ./lib

$(SHARED): $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.c))
	$(CC) -shared -o $(SODIR)/$@ $^ $(LINKER) $(OMPFLAGS)

$(ARCHIVE): $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.c))
	ar -r -o $(OUTDIR)/$@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) -c $(CFLAGS) $(INCLUDE) $(OPTFLAGS) $(OMPFLAGS) -o $@ $<

.PHONY: clean
clean:
	rm $(OBJDIR)/*.o
	rm $(OUTDIR)/$(SHARED)
	rm $(OUTDIR)/$(ARCHIVE)
