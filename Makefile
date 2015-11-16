CC = gcc
CFLAGS = -O3 -W -Wall -funroll-loops -ggdb
LIBS = -lm

BUILD = $(CC) $(CFLAGS) $(LIBS) -o $@ $?

all: mkvol

mkvol: mkvol.c
	$(BUILD)
