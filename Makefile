# Makefile for compiling NSGA-II source code
CC=gcc
LD=gcc
RM=rm -f
CFLAGS=-Wall -pedantic -g -std=gnu99
PY_CFLAGS:=$(shell python3-config --includes)
PY_LDFLAGS:=$(shell python3-config --embed --ldflags)
OBJS:=$(patsubst %.c,%.o,$(wildcard *.c))
MAIN=nsga2r
all:$(MAIN)
$(MAIN):$(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(MAIN) -lm $(PY_LDFLAGS)
%.o: %.c global.h rand.h
	$(CC) $(CFLAGS) $(PY_CFLAGS) -c $<
clean:
	$(RM) $(OBJS)
