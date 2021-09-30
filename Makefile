PATHTOCLASS = /usr/local/class_public-3.0.1
PATHTOFFTW = /usr/local/fftw-3.3.9

CC = gcc
ARCHIVE = libhmpdf.a
SHARED = libhmpdf.so

CFLAGS = --std=gnu99 -fPIC -Wall -Wextra -Wpedantic -Wno-variadic-macros -Winline -DHAVE_INLINE -DDEBUG
OPTFLAGS = -O4 -g3 -ffast-math
OMPFLAGS = -fopenmp

INCLUDE = -I./include
INCLUDE += -I$(PATHTOCLASS)/include \
           -I$(PATHTOCLASS)/external/HyRec2020 \
           -I$(PATHTOCLASS)/external/RecfastCLASS \
           -I$(PATHTOCLASS)/external/heating \
           -I$(PATHTOFFTW)/api

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

.PHONY: python
python:
	echo $(shell pwd) > hmpdf/PATHTOHMPDF.txt
	pip install . --user

.PHONY: clean
clean:
	rm $(OBJDIR)/*.o
	rmdir $(OBJDIR)
	rm $(SODIR)/$(SHARED)
	rm $(OUTDIR)/$(ARCHIVE)
	rmdir $(OUTDIR)
