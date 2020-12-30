# MAKEFILE FOR DYNUI PROGRAMS.

all:	circle

circle:		dynui.o circle.o
	cc -o circle dynui.o circle.o \
	  -lm -lcsfml-graphics -lcsfml-window -lcsfml-system -lcsfml-audio

dynui.o:	dynui.c dynui.h
circle.o:	circle.c dynui.h
