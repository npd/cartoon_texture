# Copyright 2010 Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
#		 Jose-Luis Lisani <joseluis.lisani@uib.es>
#
# Copying and distribution of this file, with or without
# modification, are permitted in any medium without royalty provided
# the copyright notice and this notice are preserved.  This file is
# offered as-is, without any warranty.


##
# Set this line to compile with OpenMP multithreading.  Comment the
# line to disable OpenMP.
OPENMP=-fopenmp

# C source code
CSRC	= iio.c
# C++ source code
CXXSRC	= main.cpp cartoon.cpp

# all source code
SRC	= $(CSRC) $(CXXSRC)

# C objects
COBJ	= $(CSRC:.c=.o)
# C++ objects
CXXOBJ	= $(CXXSRC:.cpp=.o)
# all objects
OBJ	= $(COBJ) $(CXXOBJ)
# binary target
BIN	= cartoonTexture

default	: $(BIN)


# C optimization flags
COPT	= -O3 -ftree-vectorize -funroll-loops \
         -fomit-frame-pointer  -fno-tree-pre -falign-loops -ffast-math

# C++ optimization flags
CXXOPT	= $(COPT)

# C compilation flags
CFLAGS	= $(COPT) -Wall -Wextra -std=c99 \
	-Wno-write-strings
# C++ compilation flags
CXXFLAGS	= $(CXXOPT) -Wall -Wextra  \
	-Wno-write-strings -Wno-deprecated -ansi \
        -Weffc++   -pedantic  $(OPENMP)
# link flags
LDFLAGS	= -lpng -ljpeg -ltiff $(OPENMP)


# partial compilation of C source code
%.o: %.c %.h
	$(CC) -I/opt/local/include/ -I/usr/local/include/ -c -o $@  $< $(CFLAGS)
# partial compilation of C++ source code
%.o: %.cpp %.h
	$(CXX) -I/opt/local/include/ -I/usr/local/include/ -c -o $@  $< $(CXXFLAGS)

# link all the object code
$(BIN): $(OBJ) $(LIBDEPS)
	$(CXX)  -I/opt/local/include/ -I/usr/local/include/   -L/opt/local/lib/ -L/usr/local/lib/     -o $@ $(OBJ) $(LDFLAGS)

clean:
	rm -f *.o core



