# MAKEFILE FOR DYNUI PROGRAMS.

all:	dynui-circle

dynui-circle:	dynui.o circle.o
	cc -o dynui-circle dynui.o circle.o \
	  -lm -lcsfml-graphics -lcsfml-window -lcsfml-system -lcsfml-audio

dynui.o:	dynui.c dynui.h
circle.o:	circle.c dynui.h
