CXX := g++
SANITIZER_CXX ?= clang++
ROOT_CXXFLAGS := $(shell root-config --cflags)
ROOT_LDFLAGS := $(shell root-config --libs)
WARNINGS := -Wall -Wextra -Wpedantic
PROJECT_CXX_STANDARD := -std=c++20
CXXFLAGS := -O3 -march=native $(WARNINGS) $(ROOT_CXXFLAGS) $(PROJECT_CXX_STANDARD)
LDFLAGS := $(ROOT_LDFLAGS)

# Portable, path-independent compiler profile for archived/publication builds.
# It deliberately contains neither -march=native nor -mtune=native.
REPRO_CXXFLAGS := -O3 $(WARNINGS) $(ROOT_CXXFLAGS) $(PROJECT_CXX_STANDARD) \
	-ffile-prefix-map=$(CURDIR)=. -fmacro-prefix-map=$(CURDIR)=. \
	-frandom-seed=positronium
SANITIZER_CXXFLAGS := -O1 -g $(WARNINGS) $(ROOT_CXXFLAGS) \
	$(PROJECT_CXX_STANDARD) \
	-fno-omit-frame-pointer -ffile-prefix-map=$(CURDIR)=.
PAIR ?= electron,positron
PAIR_SMOKE_ENERGY_EV ?= 20

TARGET := positronium
VALIDATION_TARGET := positronium_validation
REPRO_TARGET := positronium_reproducible
REPRO_VALIDATION_TARGET := positronium_validation_reproducible
ASAN_VALIDATION_TARGET := positronium_validation_asan
UBSAN_VALIDATION_TARGET := positronium_validation_ubsan
ASAN_PRODUCTION_TARGET := positronium_asan
UBSAN_PRODUCTION_TARGET := positronium_ubsan
LSAN_SUPPRESSIONS := $(CURDIR)/tools/lsan-root.supp
SRC := positronium.cpp
HEADERS := $(wildcard modules/*.hpp)

.PHONY: all build validation validation-small validation-publication \
	validation-pair production-pair-smoke \
	reproducible reproducible-validation sanitizers asan-check ubsan-check \
	sanitizers-check \
	references-check header-isolation-check toolchain-info run clean

all: run

build: $(TARGET)

validation: $(VALIDATION_TARGET)

validation-small: $(VALIDATION_TARGET)
	./$(VALIDATION_TARGET) --statistics-profile small

validation-publication: $(VALIDATION_TARGET)
	./$(VALIDATION_TARGET) --statistics-profile publication

# CI entry points deliberately remain separate.  A regression-suite failure
# must not be reported as a failed production trajectory (or vice versa).
validation-pair: $(VALIDATION_TARGET)
	./$(VALIDATION_TARGET) --statistics-profile small --pair "$(PAIR)"

production-pair-smoke: $(TARGET)
	./$(TARGET) --mode statistical --phenomenon 4 --runs 1 --seed 42 \
		--pair "$(PAIR)" --beam-energy-ev "$(PAIR_SMOKE_ENERGY_EV)" \
		--theta-min-deg 5 --angle-bins 10

reproducible: $(REPRO_TARGET)

reproducible-validation: $(REPRO_VALIDATION_TARGET)

sanitizers: $(ASAN_VALIDATION_TARGET) $(UBSAN_VALIDATION_TARGET) \
	$(ASAN_PRODUCTION_TARGET) $(UBSAN_PRODUCTION_TARGET)

asan-check: $(ASAN_VALIDATION_TARGET) $(ASAN_PRODUCTION_TARGET)
	ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
		LSAN_OPTIONS=suppressions=$(LSAN_SUPPRESSIONS):print_suppressions=1 \
		./$(ASAN_VALIDATION_TARGET) --statistics-profile small
	ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
		LSAN_OPTIONS=suppressions=$(LSAN_SUPPRESSIONS):print_suppressions=1 \
		./tools/run_sanitized_production_smoke.sh ./$(ASAN_PRODUCTION_TARGET)

ubsan-check: $(UBSAN_VALIDATION_TARGET) $(UBSAN_PRODUCTION_TARGET)
	UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 ./$(UBSAN_VALIDATION_TARGET) --statistics-profile small
	UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 \
		./tools/run_sanitized_production_smoke.sh ./$(UBSAN_PRODUCTION_TARGET)

sanitizers-check: asan-check ubsan-check

toolchain-info:
	@$(CXX) --version | head -n 1
	@$(SANITIZER_CXX) --version | head -n 1
	@root-config --version
	@$(MAKE) --version | head -n 1

references-check:
	python3 tools/validate_references.py ScientificalReferences.txt

# Every modules/*.hpp must compile on its own.  This is a syntax-only sweep, so
# it costs seconds; run it after touching any header's includes.
header-isolation-check:
	./tools/check_header_isolation.sh

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

$(ASAN_PRODUCTION_TARGET): $(SRC) $(HEADERS)
	$(SANITIZER_CXX) $(SANITIZER_CXXFLAGS) -fsanitize=address \
		-o $@ $(SRC) $(LDFLAGS) -fsanitize=address

$(UBSAN_PRODUCTION_TARGET): $(SRC) $(HEADERS)
	$(SANITIZER_CXX) $(SANITIZER_CXXFLAGS) -fsanitize=undefined \
		-o $@ $(SRC) $(LDFLAGS) -fsanitize=undefined

clean:
	rm -f $(TARGET) $(VALIDATION_TARGET) $(REPRO_TARGET) \
		$(REPRO_VALIDATION_TARGET) $(ASAN_VALIDATION_TARGET) \
		$(UBSAN_VALIDATION_TARGET) $(ASAN_PRODUCTION_TARGET) \
		$(UBSAN_PRODUCTION_TARGET)
