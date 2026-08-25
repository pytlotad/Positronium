CXX := g++
CXXFLAGS := -O3 -march=native -Wall -Wextra -Wpedantic $(shell root-config --cflags)
LDFLAGS := $(shell root-config --libs)

TARGET := positronium
VALIDATION_TARGET := positronium_validation
SRC := positronium.cpp
HEADERS := $(wildcard modules/*.hpp)

.PHONY: all build validation validation-small validation-publication run clean

all: run

build: $(TARGET)

validation: $(VALIDATION_TARGET)

validation-small: $(VALIDATION_TARGET)
	./$(VALIDATION_TARGET) --statistics-profile small

validation-publication: $(VALIDATION_TARGET)
	./$(VALIDATION_TARGET) --statistics-profile publication

run: $(TARGET)
	./$(TARGET)

$(TARGET): $(SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC) $(LDFLAGS)

$(VALIDATION_TARGET): $(SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -Wno-unused-function -DPOSITRONIUM_ENABLE_FIELD_VALIDATION \
		-DPOSITRONIUM_VALIDATION_EXECUTABLE -o $@ $(SRC) $(LDFLAGS)

clean:
	rm -f $(TARGET) $(VALIDATION_TARGET)
