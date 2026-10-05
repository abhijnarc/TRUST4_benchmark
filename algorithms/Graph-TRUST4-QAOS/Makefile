CXX ?= g++
CXXFLAGS ?= -O3 -Wall -Wextra -std=c++17 -fopenmp
LDFLAGS ?= -fopenmp

all: graph-trust4-qaos

graph-trust4-qaos: main.o GraphAssemblerQAOS.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

main.o: main.cpp GraphAssemblerQAOS.hpp
GraphAssemblerQAOS.o: GraphAssemblerQAOS.cpp GraphAssemblerQAOS.hpp

self-test: graph-trust4-qaos
	./graph-trust4-qaos --self-test

clean:
	rm -f *.o graph-trust4-qaos
