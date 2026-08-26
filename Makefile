CXX := g++
SANITIZER_CXX ?= clang++
ROOT_CXXFLAGS := $(shell root-config --cflags)
ROOT_LDFLAGS := $(shell root-config --libs)
WARNINGS := -Wall -Wextra -Wpedantic
CXXFLAGS := -O3 -march=native $(WARNINGS) $(ROOT_CXXFLAGS)
LDFLAGS := $(ROOT_LDFLAGS)

# Portable, path-independent compiler profile for archived/publication builds.
# It deliberately contains neither -march=native nor -mtune=native.
REPRO_CXXFLAGS := -O3 $(WARNINGS) $(ROOT_CXXFLAGS) \
	-ffile-prefix-map=$(CURDIR)=. -fmacro-prefix-map=$(CURDIR)=. \
	-frandom-seed=positronium
SANITIZER_CXXFLAGS := -O1 -g $(WARNINGS) $(ROOT_CXXFLAGS) \
	-fno-omit-frame-pointer -ffile-prefix-map=$(CURDIR)=.

TARGET := positronium
VALIDATION_TARGET := positronium_validation
REPRO_TARGET := positronium_reproducible
REPRO_VALIDATION_TARGET := positronium_validation_reproducible
ASAN_VALIDATION_TARGET := positronium_validation_asan
UBSAN_VALIDATION_TARGET := positronium_validation_ubsan
SRC := positronium.cpp
HEADERS := $(wildcard modules/*.hpp)

.PHONY: all build validation validation-small validation-publication \
	reproducible reproducible-validation sanitizers sanitizers-check \
	references-check toolchain-info run clean

all: run

build: $(TARGET)

validation: $(VALIDATION_TARGET)

validation-small: $(VALIDATION_TARGET)
	./$(VALIDATION_TARGET) --statistics-profile small

validation-publication: $(VALIDATION_TARGET)
	./$(VALIDATION_TARGET) --statistics-profile publication

reproducible: $(REPRO_TARGET)

reproducible-validation: $(REPRO_VALIDATION_TARGET)

sanitizers: $(ASAN_VALIDATION_TARGET) $(UBSAN_VALIDATION_TARGET)

sanitizers-check: sanitizers
	ASAN_OPTIONS=detect_leaks=0 ./$(ASAN_VALIDATION_TARGET) --statistics-profile small
	UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 ./$(UBSAN_VALIDATION_TARGET) --statistics-profile small

toolchain-info:
	@$(CXX) --version | head -n 1
	@$(SANITIZER_CXX) --version | head -n 1
	@root-config --version
	@$(MAKE) --version | head -n 1

references-check:
	python3 tools/validate_references.py ScientificalReferences.txt

run: $(TARGET)
	./$(TARGET)

$(TARGET): $(SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC) $(LDFLAGS)

$(VALIDATION_TARGET): $(SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -Wno-unused-function -DPOSITRONIUM_ENABLE_FIELD_VALIDATION \
		-DPOSITRONIUM_VALIDATION_EXECUTABLE -o $@ $(SRC) $(LDFLAGS)

$(REPRO_TARGET): $(SRC) $(HEADERS)
	LC_ALL=C SOURCE_DATE_EPOCH=0 $(CXX) $(REPRO_CXXFLAGS) -o $@ $(SRC) $(LDFLAGS)

$(REPRO_VALIDATION_TARGET): $(SRC) $(HEADERS)
	LC_ALL=C SOURCE_DATE_EPOCH=0 $(CXX) $(REPRO_CXXFLAGS) -Wno-unused-function \
		-DPOSITRONIUM_ENABLE_FIELD_VALIDATION -DPOSITRONIUM_VALIDATION_EXECUTABLE \
		-o $@ $(SRC) $(LDFLAGS)

$(ASAN_VALIDATION_TARGET): $(SRC) $(HEADERS)
	$(SANITIZER_CXX) $(SANITIZER_CXXFLAGS) -fsanitize=address -Wno-unused-function \
		-DPOSITRONIUM_ENABLE_FIELD_VALIDATION -DPOSITRONIUM_VALIDATION_EXECUTABLE \
		-o $@ $(SRC) $(LDFLAGS) -fsanitize=address

$(UBSAN_VALIDATION_TARGET): $(SRC) $(HEADERS)
	$(SANITIZER_CXX) $(SANITIZER_CXXFLAGS) -fsanitize=undefined -Wno-unused-function \
		-DPOSITRONIUM_ENABLE_FIELD_VALIDATION -DPOSITRONIUM_VALIDATION_EXECUTABLE \
		-o $@ $(SRC) $(LDFLAGS) -fsanitize=undefined

clean:
	rm -f $(TARGET) $(VALIDATION_TARGET) $(REPRO_TARGET) \
		$(REPRO_VALIDATION_TARGET) $(ASAN_VALIDATION_TARGET) \
		$(UBSAN_VALIDATION_TARGET)
