CXX := g++
CXXFLAGS := -O3 -march=native -Wall -Wextra -Wpedantic $(shell root-config --cflags)
LDFLAGS := $(shell root-config --libs)

TARGET := positronium
VALIDATION_TARGET := positronium_validation
SRC := positronium.cpp
HEADERS := bound_decay.hpp root_export.hpp statistics_archive.hpp \
	parameters/physical_constants.hpp objects/vector3.hpp objects/dipole_tensor.hpp objects/state.hpp \
	interactions/electrodynamics.hpp fields/maxwell_validation_backend.hpp \
	tests/maxwell_validation.hpp

.PHONY: all build validation run clean

all: run

build: $(TARGET)

validation: $(VALIDATION_TARGET)

run: $(TARGET)
	./$(TARGET)

$(TARGET): $(SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC) $(LDFLAGS)

$(VALIDATION_TARGET): $(SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -Wno-unused-function -DPOSITRONIUM_ENABLE_FIELD_VALIDATION \
		-DPOSITRONIUM_VALIDATION_EXECUTABLE -o $@ $(SRC) $(LDFLAGS)

clean:
	rm -f $(TARGET) $(VALIDATION_TARGET)
