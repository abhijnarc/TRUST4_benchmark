CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall -fopenmp
LDFLAGS = -lm -fopenmp

TARGET = graph-trust4-v2-memoryfix
SOURCES = main.cpp GraphAssemblerV2.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "Build successful: $@ (OpenMP enabled)"

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all clean
