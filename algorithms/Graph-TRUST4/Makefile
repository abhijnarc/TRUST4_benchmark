CXX = g++
CXXFLAGS = -O3 -Wall -std=c++11 -pthread

TARGETS = graph-trust4

SOURCES = main.cpp GraphAssembler.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(TARGETS)

graph-trust4: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(TARGETS)

.PHONY: all clean
