CXX := g++
CXXFLAGS := -O2 $(shell root-config --cflags)
LDFLAGS := $(shell root-config --libs)

TARGET := positronium
SRC := positronium.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(TARGET)
