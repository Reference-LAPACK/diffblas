# Unified Makefile for building all differentiated BLAS functions
# Directory structure: src/ (sources), test/ (test programs), include/ (headers)

# Compilers and flags
FC = gfortran
CC = gcc
FFLAGS = -O2 -fPIC -ffree-line-length-none -Wuninitialized -Wmaybe-uninitialized -Iinclude
FFLAGS_F77 = -O2 -fPIC -ffixed-line-length-none -Wuninitialized -Wmaybe-uninitialized -Iinclude
CFLAGS = -O2 -fPIC

# Directory structure
SRC_DIR = src
TEST_DIR = test
INC_DIR = include
BUILD_DIR = build

# Create build directory
$(shell mkdir -p $(BUILD_DIR))

# Auto-detect functions from differentiated source files in src/
FUNCS_D := $(sort $(patsubst $(SRC_DIR)/%_d.f,%,$(wildcard $(SRC_DIR)/*_d.f)) \
                  $(patsubst $(SRC_DIR)/%_d.f90,%,$(wildcard $(SRC_DIR)/*_d.f90)))
FUNCS_B := $(sort $(patsubst $(SRC_DIR)/%_b.f,%,$(wildcard $(SRC_DIR)/*_b.f)) \
                  $(patsubst $(SRC_DIR)/%_b.f90,%,$(wildcard $(SRC_DIR)/*_b.f90)))
FUNCS_DV := $(sort $(patsubst $(SRC_DIR)/%_dv.f,%,$(wildcard $(SRC_DIR)/*_dv.f)) \
                   $(patsubst $(SRC_DIR)/%_dv.f90,%,$(wildcard $(SRC_DIR)/*_dv.f90)))
FUNCS_BV := $(sort $(patsubst $(SRC_DIR)/%_bv.f,%,$(wildcard $(SRC_DIR)/*_bv.f)) \
                   $(patsubst $(SRC_DIR)/%_bv.f90,%,$(wildcard $(SRC_DIR)/*_bv.f90)))

# All unique functions
FUNCS := $(sort $(FUNCS_D) $(FUNCS_B) $(FUNCS_DV) $(FUNCS_BV))

# Auto-detect test sources in test/, but ONLY for functions that have differentiated sources
# This filters out non-differentiable functions (like isamax, icamax, etc.)
TESTS_D_SRCS := $(wildcard $(TEST_DIR)/test_*.f90)
TESTS_D_ONLY := $(filter-out %_reverse.f90 %_vector_forward.f90 %_vector_reverse.f90,$(TESTS_D_SRCS))
# Filter to only include tests for functions in FUNCS_D
TESTS_D := $(foreach f,$(FUNCS_D),$(if $(wildcard $(TEST_DIR)/test_$(f).f90),$(BUILD_DIR)/test_$(f)))

TESTS_B_SRCS := $(wildcard $(TEST_DIR)/test_*_reverse.f90)
# Filter to only include tests for functions in FUNCS_B
TESTS_B := $(foreach f,$(FUNCS_B),$(if $(wildcard $(TEST_DIR)/test_$(f)_reverse.f90),$(BUILD_DIR)/test_$(f)_reverse))

TESTS_DV_SRCS := $(wildcard $(TEST_DIR)/test_*_vector_forward.f90)
# Filter to only include tests for functions in FUNCS_DV
TESTS_DV := $(foreach f,$(FUNCS_DV),$(if $(wildcard $(TEST_DIR)/test_$(f)_vector_forward.f90),$(BUILD_DIR)/test_$(f)_vector_forward))

TESTS_BV_SRCS := $(wildcard $(TEST_DIR)/test_*_vector_reverse.f90)
# Filter to only include tests for functions in FUNCS_BV
TESTS_BV := $(foreach f,$(FUNCS_BV),$(if $(wildcard $(TEST_DIR)/test_$(f)_vector_reverse.f90),$(BUILD_DIR)/test_$(f)_vector_reverse))

# Auto-detect original source files in src/ (for status display)
ORIG_SRCS_F := $(wildcard $(addprefix $(SRC_DIR)/,$(addsuffix .f,$(FUNCS))))
ORIG_SRCS_F90 := $(wildcard $(addprefix $(SRC_DIR)/,$(addsuffix .f90,$(FUNCS))))
ORIG_FUNCS := $(sort $(patsubst $(SRC_DIR)/%.f,%,$(ORIG_SRCS_F)) $(patsubst $(SRC_DIR)/%.f90,%,$(ORIG_SRCS_F90)))

# BLAS library for linking (lsame, xerbla, and original BLAS routines)
# Uses LAPACKDIR environment variable if set, otherwise assumes system BLAS
LAPACKDIR ?= $(shell echo $$LAPACKDIR)
ifneq ($(LAPACKDIR),)
BLAS_LIB ?= -L$(LAPACKDIR) -lrefblas
else
BLAS_LIB ?= -lrefblas
endif

# Optional: DIFFSIZES_access.o when using F77 ISIZE globals (run_tapenade_blas.py writes DIFFSIZES_access.f)
# Must be defined before any rule that uses it as a prerequisite, so "make forward" (etc.) builds it first.
ifneq ($(wildcard $(SRC_DIR)/DIFFSIZES_access.f),)
DIFFSIZES_ACCESS_OBJ := $(BUILD_DIR)/DIFFSIZES_access.o
else
DIFFSIZES_ACCESS_OBJ :=
endif

# Unified library targets (one library per mode containing all differentiated code)
# Note: Original BLAS functions come from $(BLAS_LIB) (librefblas in LAPACKDIR)
LIB_D := $(BUILD_DIR)/libdiffblas_d.a
LIB_B := $(BUILD_DIR)/libdiffblas_b.a
LIB_DV := $(BUILD_DIR)/libdiffblas_dv.a
LIB_BV := $(BUILD_DIR)/libdiffblas_bv.a

SHARED_D := $(BUILD_DIR)/libdiffblas_d.so
SHARED_B := $(BUILD_DIR)/libdiffblas_b.so
SHARED_DV := $(BUILD_DIR)/libdiffblas_dv.so
SHARED_BV := $(BUILD_DIR)/libdiffblas_bv.so

# Default target - build all modes
all: forward reverse vector-forward vector-reverse

# Mode-specific targets (unified libraries + tests)
# Original BLAS comes from $(BLAS_LIB) - no need to build liborigblas
forward: $(LIB_D) $(SHARED_D) $(TESTS_D)
	@echo "Forward mode build complete. Library: libdiffblas_d. Tests: $(words $(TESTS_D))"

reverse: $(LIB_B) $(SHARED_B) $(TESTS_B)
	@echo "Reverse mode build complete. Library: libdiffblas_b. Tests: $(words $(TESTS_B))"

vector-forward: $(LIB_DV) $(SHARED_DV) $(TESTS_DV)
	@echo "Vector forward mode build complete. Library: libdiffblas_dv. Tests: $(words $(TESTS_DV))"

vector-reverse: $(LIB_BV) $(SHARED_BV) $(TESTS_BV)
	@echo "Vector reverse mode build complete. Library: libdiffblas_bv. Tests: $(words $(TESTS_BV))"

# Libraries-only targets (no tests)
libs: libs-forward libs-reverse libs-vector-forward libs-vector-reverse

libs-forward: $(LIB_D) $(SHARED_D)
	@echo "Forward mode library complete: libdiffblas_d"

libs-reverse: $(LIB_B) $(SHARED_B)
	@echo "Reverse mode library complete: libdiffblas_b"

libs-vector-forward: $(LIB_DV) $(SHARED_DV)
	@echo "Vector forward mode library complete: libdiffblas_dv"

libs-vector-reverse: $(LIB_BV) $(SHARED_BV)
	@echo "Vector reverse mode library complete: libdiffblas_bv"

# =============================================================================
# Pattern rules for Fortran compilation (src/ -> build/)
# Use - prefix to continue on errors (some files may fail to compile)
# =============================================================================

# Forward mode objects
$(BUILD_DIR)/%_d.o: $(SRC_DIR)/%_d.f
	-$(FC) $(FFLAGS_F77) -c $< -o $@

$(BUILD_DIR)/%_d.o: $(SRC_DIR)/%_d.f90
	-$(FC) $(FFLAGS) -c $< -o $@

# Reverse mode objects
$(BUILD_DIR)/%_b.o: $(SRC_DIR)/%_b.f
	-$(FC) $(FFLAGS_F77) -c $< -o $@

$(BUILD_DIR)/%_b.o: $(SRC_DIR)/%_b.f90
	-$(FC) $(FFLAGS) -c $< -o $@

# Vector forward mode objects
$(BUILD_DIR)/%_dv.o: $(SRC_DIR)/%_dv.f $(BUILD_DIR)/DIFFSIZES.o
	-$(FC) $(FFLAGS_F77) -c $< -o $@

$(BUILD_DIR)/%_dv.o: $(SRC_DIR)/%_dv.f90 $(BUILD_DIR)/DIFFSIZES.o
	-$(FC) $(FFLAGS) -c $< -o $@

# Vector reverse mode objects
$(BUILD_DIR)/%_bv.o: $(SRC_DIR)/%_bv.f $(BUILD_DIR)/DIFFSIZES.o
	-$(FC) $(FFLAGS_F77) -c $< -o $@

$(BUILD_DIR)/%_bv.o: $(SRC_DIR)/%_bv.f90 $(BUILD_DIR)/DIFFSIZES.o
	-$(FC) $(FFLAGS) -c $< -o $@

# Original function objects
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f
	$(FC) $(FFLAGS_F77) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Dependency objects (xerbla_dep*.f, etc.)
$(BUILD_DIR)/%_dep.o: $(SRC_DIR)/%_dep.f
	$(FC) $(FFLAGS_F77) -c $< -o $@

$(BUILD_DIR)/%_dep1.o: $(SRC_DIR)/%_dep1.f
	$(FC) $(FFLAGS_F77) -c $< -o $@

$(BUILD_DIR)/%_dep2.o: $(SRC_DIR)/%_dep2.f
	$(FC) $(FFLAGS_F77) -c $< -o $@

# DIFFSIZES_access.f - global ISIZE storage and get/set/check (for _b, _bv when using F77 DIFFSIZES.inc)
$(BUILD_DIR)/DIFFSIZES_access.o: $(SRC_DIR)/DIFFSIZES_access.f
	$(FC) $(FFLAGS_F77) -c $< -o $@

# DIFFSIZES handling (supports both Fortran 90 module and Fortran 77 include)
# For F90: DIFFSIZES.f90 is compiled to produce DIFFSIZES.o and DIFFSIZES.mod
# For F77: DIFFSIZES.inc is included inline - no object file needed
# We use a marker file to track that F77 mode is being used
$(BUILD_DIR)/DIFFSIZES.o:
	@mkdir -p $(BUILD_DIR)
	@if [ -f $(INC_DIR)/DIFFSIZES.f90 ]; then \
		echo "Compiling DIFFSIZES.f90 module"; \
		$(FC) $(FFLAGS) -c $(INC_DIR)/DIFFSIZES.f90 -o $@; \
	elif [ -f $(INC_DIR)/DIFFSIZES.inc ]; then \
		echo "Using DIFFSIZES.inc (Fortran 77 include - no object needed)"; \
		echo "      subroutine diffsizes_stub_unused" > $(BUILD_DIR)/.diffsizes_stub.f; \
		echo "      end" >> $(BUILD_DIR)/.diffsizes_stub.f; \
		$(FC) $(FFLAGS_F77) -c $(BUILD_DIR)/.diffsizes_stub.f -o $@; \
		rm -f $(BUILD_DIR)/.diffsizes_stub.f; \
	else \
		echo "WARNING: No DIFFSIZES file found in $(INC_DIR)/"; \
		echo "      subroutine diffsizes_stub_unused" > $(BUILD_DIR)/.diffsizes_stub.f; \
		echo "      end" >> $(BUILD_DIR)/.diffsizes_stub.f; \
		$(FC) $(FFLAGS_F77) -c $(BUILD_DIR)/.diffsizes_stub.f -o $@; \
		rm -f $(BUILD_DIR)/.diffsizes_stub.f; \
	fi

# adStack for reverse mode (Tapenade runtime library)
# Try local copy first, then fall back to repo copy, then TAPENADEDIR
$(BUILD_DIR)/adStack.o: 
	@if [ -f $(SRC_DIR)/adStack.c ]; then \
		$(CC) $(CFLAGS) -c $(SRC_DIR)/adStack.c -o $@; \
	elif [ -f "$(CURDIR)/../TAPENADE/adStack.c" ]; then \
		$(CC) $(CFLAGS) -I"$(CURDIR)/../TAPENADE/include" -c "$(CURDIR)/../TAPENADE/adStack.c" -o $@; \
	elif [ -n "$$TAPENADEDIR" ] && [ -f "$$TAPENADEDIR/ADFirstAidKit/adStack.c" ]; then \
		$(CC) $(CFLAGS) -I$$TAPENADEDIR/ADFirstAidKit -c $$TAPENADEDIR/ADFirstAidKit/adStack.c -o $@; \
	else \
		echo "ERROR: adStack.c not found. Set TAPENADEDIR or copy to $(SRC_DIR)/"; \
		exit 1; \
	fi

# =============================================================================
# Compile all sources for each mode (phony targets, continue on errors)
# =============================================================================

# Source files for each mode
SRCS_D := $(wildcard $(SRC_DIR)/*_d.f) $(wildcard $(SRC_DIR)/*_d.f90)
SRCS_B := $(wildcard $(SRC_DIR)/*_b.f) $(wildcard $(SRC_DIR)/*_b.f90)
SRCS_DV := $(wildcard $(SRC_DIR)/*_dv.f) $(wildcard $(SRC_DIR)/*_dv.f90)
SRCS_BV := $(wildcard $(SRC_DIR)/*_bv.f) $(wildcard $(SRC_DIR)/*_bv.f90)

.PHONY: compile-d compile-b compile-dv compile-bv

compile-d:
	@echo "Compiling forward mode sources ($(words $(SRCS_D)) files)..."
	@failed=0; for src in $(SRCS_D); do \
		obj=$(BUILD_DIR)/$$(basename $${src%.*}).o; \
		if [ ! -f "$$obj" ] || [ "$$src" -nt "$$obj" ]; then \
			if echo "$$src" | grep -q '\.f90$$'; then \
				$(FC) $(FFLAGS) -c "$$src" -o "$$obj" 2>/dev/null || { echo "  FAILED: $$src"; failed=$$((failed+1)); }; \
			else \
				$(FC) $(FFLAGS_F77) -c "$$src" -o "$$obj" 2>/dev/null || { echo "  FAILED: $$src"; failed=$$((failed+1)); }; \
			fi; \
		fi; \
	done; echo "  Compiled $$(ls $(BUILD_DIR)/*_d.o 2>/dev/null | wc -w) objects ($$failed failed)"

compile-b: $(BUILD_DIR)/adStack.o
	@echo "Compiling reverse mode sources ($(words $(SRCS_B)) files)..."
	@failed=0; for src in $(SRCS_B); do \
		obj=$(BUILD_DIR)/$$(basename $${src%.*}).o; \
		if [ ! -f "$$obj" ] || [ "$$src" -nt "$$obj" ]; then \
			if echo "$$src" | grep -q '\.f90$$'; then \
				$(FC) $(FFLAGS) -c "$$src" -o "$$obj" 2>/dev/null || { echo "  FAILED: $$src"; failed=$$((failed+1)); }; \
			else \
				$(FC) $(FFLAGS_F77) -c "$$src" -o "$$obj" 2>/dev/null || { echo "  FAILED: $$src"; failed=$$((failed+1)); }; \
			fi; \
		fi; \
	done; echo "  Compiled $$(ls $(BUILD_DIR)/*_b.o 2>/dev/null | wc -w) objects ($$failed failed)"

compile-dv: $(BUILD_DIR)/DIFFSIZES.o
	@echo "Compiling vector forward mode sources ($(words $(SRCS_DV)) files)..."
	@failed=0; for src in $(SRCS_DV); do \
		obj=$(BUILD_DIR)/$$(basename $${src%.*}).o; \
		if [ ! -f "$$obj" ] || [ "$$src" -nt "$$obj" ]; then \
			if echo "$$src" | grep -q '\.f90$$'; then \
				$(FC) $(FFLAGS) -c "$$src" -o "$$obj" 2>/dev/null || { echo "  FAILED: $$src"; failed=$$((failed+1)); }; \
			else \
				$(FC) $(FFLAGS_F77) -c "$$src" -o "$$obj" 2>/dev/null || { echo "  FAILED: $$src"; failed=$$((failed+1)); }; \
			fi; \
		fi; \
	done; echo "  Compiled $$(ls $(BUILD_DIR)/*_dv.o 2>/dev/null | wc -w) objects ($$failed failed)"

compile-bv: $(BUILD_DIR)/adStack.o $(BUILD_DIR)/DIFFSIZES.o
	@echo "Compiling vector reverse mode sources ($(words $(SRCS_BV)) files)..."
	@failed=0; for src in $(SRCS_BV); do \
		obj=$(BUILD_DIR)/$$(basename $${src%.*}).o; \
		if [ ! -f "$$obj" ] || [ "$$src" -nt "$$obj" ]; then \
			if echo "$$src" | grep -q '\.f90$$'; then \
				$(FC) $(FFLAGS) -c "$$src" -o "$$obj" 2>/dev/null || { echo "  FAILED: $$src"; failed=$$((failed+1)); }; \
			else \
				$(FC) $(FFLAGS_F77) -c "$$src" -o "$$obj" 2>/dev/null || { echo "  FAILED: $$src"; failed=$$((failed+1)); }; \
			fi; \
		fi; \
	done; echo "  Compiled $$(ls $(BUILD_DIR)/*_bv.o 2>/dev/null | wc -w) objects ($$failed failed)"

# =============================================================================
# Unified libraries - built from whatever objects successfully compiled
# =============================================================================

# Single library for all forward mode differentiated code
$(BUILD_DIR)/libdiffblas_d.a: compile-d $(DIFFSIZES_ACCESS_OBJ)
	@ar rcs $@ $$(ls $(BUILD_DIR)/*_d.o 2>/dev/null) $(DIFFSIZES_ACCESS_OBJ)
	@echo "Created libdiffblas_d.a with $$(ls $(BUILD_DIR)/*_d.o 2>/dev/null | wc -w) objects"

$(BUILD_DIR)/libdiffblas_d.so: compile-d
	@$(FC) -shared -o $@ $$(ls $(BUILD_DIR)/*_d.o 2>/dev/null)

# Single library for all reverse mode differentiated code
$(BUILD_DIR)/libdiffblas_b.a: compile-b $(DIFFSIZES_ACCESS_OBJ)
	@ar rcs $@ $$(ls $(BUILD_DIR)/*_b.o 2>/dev/null) $(BUILD_DIR)/adStack.o $(DIFFSIZES_ACCESS_OBJ)
	@echo "Created libdiffblas_b.a with $$(ls $(BUILD_DIR)/*_b.o 2>/dev/null | wc -w) objects"

$(BUILD_DIR)/libdiffblas_b.so: compile-b $(DIFFSIZES_ACCESS_OBJ)
	@$(FC) -shared -o $@ $$(ls $(BUILD_DIR)/*_b.o 2>/dev/null) $(BUILD_DIR)/adStack.o $(DIFFSIZES_ACCESS_OBJ)

# Single library for all vector forward mode differentiated code
$(BUILD_DIR)/libdiffblas_dv.a: compile-dv $(DIFFSIZES_ACCESS_OBJ)
	@ar rcs $@ $$(ls $(BUILD_DIR)/*_dv.o 2>/dev/null) $(BUILD_DIR)/DIFFSIZES.o $(DIFFSIZES_ACCESS_OBJ)
	@echo "Created libdiffblas_dv.a with $$(ls $(BUILD_DIR)/*_dv.o 2>/dev/null | wc -w) objects"

$(BUILD_DIR)/libdiffblas_dv.so: compile-dv
	@$(FC) -shared -o $@ $$(ls $(BUILD_DIR)/*_dv.o 2>/dev/null) $(BUILD_DIR)/DIFFSIZES.o

# Single library for all vector reverse mode differentiated code
$(BUILD_DIR)/libdiffblas_bv.a: compile-bv $(DIFFSIZES_ACCESS_OBJ)
	@ar rcs $@ $$(ls $(BUILD_DIR)/*_bv.o 2>/dev/null) $(BUILD_DIR)/adStack.o $(BUILD_DIR)/DIFFSIZES.o $(DIFFSIZES_ACCESS_OBJ)
	@echo "Created libdiffblas_bv.a with $$(ls $(BUILD_DIR)/*_bv.o 2>/dev/null | wc -w) objects"

$(BUILD_DIR)/libdiffblas_bv.so: compile-bv $(DIFFSIZES_ACCESS_OBJ)
	@$(FC) -shared -o $@ $$(ls $(BUILD_DIR)/*_bv.o 2>/dev/null) $(BUILD_DIR)/adStack.o $(BUILD_DIR)/DIFFSIZES.o $(DIFFSIZES_ACCESS_OBJ)

# Note: Original BLAS functions come from $(BLAS_LIB) (librefblas in LAPACKDIR)
# No need to build a separate liborigblas

 $<

# =============================================================================
# Pattern rules for test executables (test/ -> build/)
# Link against unified libraries for efficiency
# =============================================================================

# Forward mode tests (link against libdiffblas_d + BLAS_LIB for original functions)
# libdiffblas_d.a already contains DIFFSIZES_access.o; do not link it again to avoid multiple definition
$(BUILD_DIR)/test_%: $(BUILD_DIR)/test_%.o $(BUILD_DIR)/libdiffblas_d.a
	-$(FC) $< $(BUILD_DIR)/libdiffblas_d.a $(BLAS_LIB) -o $@

$(BUILD_DIR)/test_%.o: $(TEST_DIR)/test_%.f90
	-$(FC) $(FFLAGS) -c $< -o $@

# Reverse mode tests (link against libdiffblas_b + BLAS_LIB for original functions)
$(BUILD_DIR)/test_%_reverse: $(BUILD_DIR)/test_%_reverse.o $(BUILD_DIR)/libdiffblas_b.a
	-$(FC) $< $(BUILD_DIR)/libdiffblas_b.a $(BLAS_LIB) -o $@

$(BUILD_DIR)/test_%_reverse.o: $(TEST_DIR)/test_%_reverse.f90
	-$(FC) $(FFLAGS) -c $< -o $@

# Vector forward mode tests (link against libdiffblas_dv + BLAS_LIB for original functions)
# libdiffblas_dv.a already contains DIFFSIZES_access.o; do not link it again to avoid multiple definition
$(BUILD_DIR)/test_%_vector_forward: $(BUILD_DIR)/test_%_vector_forward.o $(BUILD_DIR)/libdiffblas_dv.a
	-$(FC) $< $(BUILD_DIR)/libdiffblas_dv.a $(BLAS_LIB) -o $@

$(BUILD_DIR)/test_%_vector_forward.o: $(TEST_DIR)/test_%_vector_forward.f90 $(BUILD_DIR)/DIFFSIZES.o
	-$(FC) $(FFLAGS) -c $< -o $@

# Vector reverse mode tests (link against libdiffblas_bv + BLAS_LIB for original functions)
$(BUILD_DIR)/test_%_vector_reverse: $(BUILD_DIR)/test_%_vector_reverse.o $(BUILD_DIR)/libdiffblas_bv.a
	-$(FC) $< $(BUILD_DIR)/libdiffblas_bv.a $(BLAS_LIB) -o $@

$(BUILD_DIR)/test_%_vector_reverse.o: $(TEST_DIR)/test_%_vector_reverse.f90 $(BUILD_DIR)/DIFFSIZES.o
	-$(FC) $(FFLAGS) -c $< -o $@

# =============================================================================
# Utility targets
# =============================================================================

# Clean all build artifacts
clean:
	@echo "Cleaning build directory..."
	rm -rf $(BUILD_DIR)
	@echo "Clean complete."

# Rebuild everything
rebuild: clean all

# Run all tests
test: all
	@echo "Running all tests..."
	@passed=0; failed=0; \
	for t in $(TESTS_D) $(TESTS_B) $(TESTS_DV) $(TESTS_BV); do \
		if [ -f "$$t" ]; then \
			echo -n "Running $$(basename $$t)... "; \
			if $$t > /dev/null 2>&1; then \
				echo "PASS"; \
				passed=$$((passed + 1)); \
			else \
				echo "FAIL"; \
				failed=$$((failed + 1)); \
			fi; \
		fi; \
	done; \
	echo ""; \
	echo "Results: $$passed passed, $$failed failed"

# Show status
status:
	@echo "=== Directory Structure ==="
	@echo "  Sources:  $(SRC_DIR)/"
	@echo "  Tests:    $(TEST_DIR)/"
	@echo "  Python:   $(TEST_DIR)/python/"
	@echo "  Includes: $(INC_DIR)/"
	@echo "  Build:    $(BUILD_DIR)/"
	@echo ""
	@echo "=== Source Files Detected ==="
	@echo "Differentiated functions: $(words $(FUNCS))"
	@echo "  Forward mode (d):  $(words $(FUNCS_D)) sources"
	@echo "  Reverse mode (b):  $(words $(FUNCS_B)) sources"
	@echo "  Vector fwd (dv):   $(words $(FUNCS_DV)) sources"
	@echo "  Vector rev (bv):   $(words $(FUNCS_BV)) sources"
	@echo ""
	@echo "Original sources:    $(words $(ORIG_FUNCS)) files"
	@echo "BLAS library:        $(BLAS_LIB)"
	@echo ""
	@echo "Test sources found:"
	@echo "  Forward tests:     $(words $(TESTS_D))"
	@echo "  Reverse tests:     $(words $(TESTS_B))"
	@echo "  Vector fwd tests:  $(words $(TESTS_DV))"
	@echo "  Vector rev tests:  $(words $(TESTS_BV))"
	@echo -n "  Python interfaces: "; ls -1 $(TEST_DIR)/python/*.py 2>/dev/null | wc -l
	@echo ""
	@echo "=== Build Artifacts ==="
	@echo "Unified libraries (one per mode):"
	@echo -n "  libdiffblas_d:   "; [ -f $(BUILD_DIR)/libdiffblas_d.a ] && echo "built" || echo "not built"
	@echo -n "  libdiffblas_b:   "; [ -f $(BUILD_DIR)/libdiffblas_b.a ] && echo "built" || echo "not built"
	@echo -n "  libdiffblas_dv:  "; [ -f $(BUILD_DIR)/libdiffblas_dv.a ] && echo "built" || echo "not built"
	@echo -n "  libdiffblas_bv:  "; [ -f $(BUILD_DIR)/libdiffblas_bv.a ] && echo "built" || echo "not built"
	@echo "  Original BLAS:   from $(BLAS_LIB)"
	@echo -n "  Test executables:  "; ls -1 $(BUILD_DIR)/test_* 2>/dev/null | grep -v '\.o$$' | wc -l

# Help
help:
	@echo "Unified Makefile for Differentiated BLAS Functions"
	@echo ""
	@echo "Directory structure:"
	@echo "  src/         - Differentiated Fortran sources (*_d.f, *_b.f, etc.)"
	@echo "  test/        - Test programs (test_*.f90)"
	@echo "  test/python/ - Python interfaces (*_original.py, *_differentiated.py)"
	@echo "  include/     - Header files (DIFFSIZES.inc, DIFFSIZES.f90)"
	@echo "  build/       - Build artifacts (*.o, *.a, *.so, test executables)"
	@echo ""
	@echo "Build targets:"
	@echo "  all              - Build all modes with tests (default)"
	@echo "  forward          - Build forward mode libraries + tests"
	@echo "  reverse          - Build reverse mode libraries + tests"
	@echo "  vector-forward   - Build vector forward mode libraries + tests"
	@echo "  vector-reverse   - Build vector reverse mode libraries + tests"
	@echo ""
	@echo "Library-only targets (no tests):"
	@echo "  libs             - Build all libraries (no tests)"
	@echo "  libs-forward     - Build forward mode libraries only"
	@echo "  libs-reverse     - Build reverse mode libraries only"
	@echo "  libs-vector-forward - Build vector forward libraries only"
	@echo "  libs-vector-reverse - Build vector reverse libraries only"
	@echo ""
	@echo "Utility targets:"
	@echo "  clean            - Remove build/ directory"
	@echo "  rebuild          - Clean and rebuild everything"
	@echo "  test             - Run all test executables"
	@echo "  status           - Show detected files and build status"
	@echo "  help             - Show this help"
	@echo ""
	@echo "Detected $(words $(FUNCS)) functions"

.PHONY: all forward reverse vector-forward vector-reverse clean rebuild test status help
.PHONY: libs libs-forward libs-reverse libs-vector-forward libs-vector-reverse
