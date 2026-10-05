CXX := g++
CXXFLAGS := -std=c++11 -O3 -Wall -Wextra -pedantic -fopenmp
LDFLAGS := -fopenmp -lm
TARGET := graph-trust4-v3
SOURCES := main.cpp GraphAssemblerV3.cpp
OBJECTS := $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp GraphAssemblerV3.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all clean
