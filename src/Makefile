# ***************************************************************************
# Timothy - Tissue Modelling Framework
# Copyright (C) 2014-2018 Maciej Cytowski
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# ***************************************************************************

CC=mpicc
CFLAGS=-O2 -g -fopenmp  -ffast-math -ftree-vectorize -finline-functions -Winline
LDFLAGS = -fopenmp
DEFS=-DHYPRE #-DDEBUG

ZOLTANLIB=-LPATH_TO_ZOLTAN_LIB -lzoltan
ZOLTANINC=-IPATH_TO_ZOLTAN_INCLUDE

HYPRELIB=-LPATH_TO_HYPRE_LIB -lHYPRE -lblas -llapack
HYPREINC=-IPATH_TO_HYPRE_INCLUDE

GMPLIB=-lgmp

CFLAGS += -I. $(ZOLTANINC) $(HYPREINC)
LDFLAGS += -lm $(ZOLTANLIB) $(GMPLIB) $(HYPRELIB)

OBJ = main.o io.o utils.o cells.o lb.o octree.o potential.o exchange.o initialisation.o grid.o patches.o fieldgradient.o environment.o step.o
DEPS = global.h io.h inline.h fieldgradient.h potential.h octree.h patches.h environment.h

all: timothy

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) $(DEFS) -o $@ -c $<

timothy: $(OBJ)
	$(CC) $(LDFLAGS) $(OBJ) -o $@

clean:
	rm -f $(OBJ); rm -f timothy

install:
	cp timothy ../bin
