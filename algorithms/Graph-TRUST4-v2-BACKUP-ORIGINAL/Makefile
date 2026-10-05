CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall
LDFLAGS = -lm

TARGET = graph-trust4-v2
SOURCES = main.cpp GraphAssemblerV2.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "Build successful: $@"

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all clean
