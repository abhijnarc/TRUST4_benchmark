CXX ?= g++
CXXFLAGS ?= -O3 -Wall -Wextra -std=c++17 -fopenmp
LDFLAGS ?= -fopenmp

all: graph-trust4-qaos-clean

graph-trust4-qaos-clean: main.o clean_graph.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

main.o: main.cpp clean_graph.hpp
clean_graph.o: clean_graph.cpp clean_graph.hpp

clean:
	rm -f *.o graph-trust4-qaos-clean

.PHONY: all clean
