# MAKEFILE FOR DYNUI PROGRAMS.

CFLAGS=-g -O3 -march=native -mtune=native -Wno-logical-op-parentheses

all:	circle LJ2D

dynui.o:	dynui.c dynui.h

circle:		dynui.o circle.o
	cc -o circle dynui.o circle.o \
	  -lm -lcsfml-graphics -lcsfml-window -lcsfml-system -lcsfml-audio

circle.o:	circle.c dynui.h

LJ2D:		dynui.o LJ2D.o
	cc -o LJ2D dynui.o LJ2D.o \
	  -lm -lcsfml-graphics -lcsfml-window -lcsfml-system -lcsfml-audio

LJ2D.o:		LJ2D.c dynui.h merge-sort.c quick-sort.c
