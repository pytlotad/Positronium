CXX := g++
CXXFLAGS := -O2 $(shell root-config --cflags)
LDFLAGS := $(shell root-config --libs)

TARGET := positronium
SRC := positronium.cpp

.PHONY: all build run clean

all: run

build: $(TARGET)

run: $(TARGET)
	./$(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(TARGET)
