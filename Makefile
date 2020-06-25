PATHTOCLASS = /usr/local/share/class_public


CC = gcc
ARCHIVE = libhmpdf.a
SHARED = libhmpdf.so

CFLAGS = -std=gnu99 -fPIC -Wall -Wpedantic -Wno-variadic-macros -Wno-format -DHAVE_INLINE
OPTFLAGS = -O4 -ffast-math
OMPFLAGS = -fopenmp

INCLUDE = -I./include
INCLUDE += -I$(PATHTOCLASS)/include

LINKER = -L$(PATHTOCLASS)
LINKER += -lclass -lgsl -lgslcblas -lm -lfftw3

SRCDIR = ./src
OBJDIR = ./obj
OUTDIR = ./lib
SODIR  = .

.PHONY: directories

all: directories $(SHARED) $(ARCHIVE)

$(SHARED): $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.c))
	$(CC) -shared -o $(SODIR)/$@ $^ $(LINKER) $(OMPFLAGS)

$(ARCHIVE): $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.c))
	ar -r -o $(OUTDIR)/$@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) -c $(CFLAGS) $(INCLUDE) $(OPTFLAGS) $(OMPFLAGS) -o $@ $<

directories: $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)
	mkdir -p $(OUTDIR)

.PHONY: clean
clean:
	rm $(OBJDIR)/*.o
	rmdir $(OBJDIR)
	rm $(SODIR)/$(SHARED)
	rm $(OUTDIR)/$(ARCHIVE)
	rmdir $(OUTDIR)