# Makefile for compiling NSGA-II source code with C++
CXX=g++
LD=g++
RM=rm -f
CXXFLAGS=-Wall -Wextra -std=c++17 -g
OBJS:=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
MAIN=nsga2r

all: $(MAIN)

$(MAIN): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(MAIN) -lm

%.o: %.cpp global.h rand.h
	$(CXX) $(CXXFLAGS) -c $<

clean:
	$(RM) $(OBJS)
