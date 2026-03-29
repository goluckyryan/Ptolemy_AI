# Makefile for Ptolemy DWBA Code (Fortran)
#
# Tested on Ubuntu 24.04 with gfortran 13/14
#
# ⚠️  WARNING: The 64-bit build compiles and runs, but produces
#     INCORRECT results for some reactions. Always use the 32-bit
#     build or the pre-built 32-bit binary for production work.
#     See README.md for details.
#
#   make             — 32-bit build (RECOMMENDED, requires gcc-multilib)
#   make ptolemy64   — 64-bit build (for debugging ONLY, requires -O0)
#   make test        — run test and verify DCS output
#   make clean

FC      = gfortran
CC      = gcc

# Source files
FSRC    = ptolemy-main.f source.f fortlib.f gfortran_stuff.f masstable.f \
          numbered_store.f linkule.f linkulesfitters.f keep.f keepsub.f \
          keptsub.f phiffer.f av18.f
CSRC    = dtime.c srread.c

FOBJ    = $(FSRC:.f=.o)
COBJ    = $(CSRC:.c=.o)

# ============================================================
# Default: 32-bit build (matches original development platform)
# ============================================================
all: ptolemy

ptolemy:
	@if $(CC) -m32 -xc -c /dev/null -o /dev/null 2>/dev/null; then \
		echo "=== Building 32-bit (recommended) ==="; \
	else \
		echo "ERROR: 32-bit support not found."; \
		echo "  sudo apt install gcc-multilib gfortran-multilib"; \
		exit 1; \
	fi
	@for f in $(FSRC); do \
		echo "  $$f"; \
		$(FC) -c -m32 -O2 -std=legacy -fallow-argument-mismatch \
			-fno-range-check -fallow-invalid-boz -w $$f || exit 1; \
	done
	$(CC) -c -m32 -w dtime.c
	$(CC) -c -m32 -w srread.c
	$(FC) -m32 -o ptolemy $(FOBJ) $(COBJ)
	@echo "=== Build complete: ./ptolemy (32-bit) ==="

# ============================================================
# 64-bit build: for DEBUGGING ONLY — results may be INCORRECT
# ============================================================
# Fortran flags for legacy code compatibility
# -O0 is REQUIRED on 64-bit (code breaks at -O1+)
FFLAGS64 = -O0 -std=legacy -fallow-argument-mismatch -fno-range-check \
           -fallow-invalid-boz -fd-lines-as-comments -w

ptolemy64: $(FSRC) $(CSRC)
	@echo "⚠️  WARNING: 64-bit build is for DEBUGGING ONLY."
	@echo "    Results may differ from the correct 32-bit output."
	@echo "    Use 'make' (32-bit) for production work."
	@echo ""
	@for f in $(FSRC); do \
		echo "  $$f"; \
		$(FC) -c $(FFLAGS64) $$f || exit 1; \
	done
	$(CC) -c -w dtime.c
	$(CC) -c -w srread.c
	$(FC) -o ptolemy $(FOBJ) $(COBJ)
	@echo "=== Build complete: ./ptolemy (64-bit, DEBUG ONLY) ==="

# Pattern rules (for ptolemy64 only — 32-bit build uses explicit loop)
%.o: %.f
	$(FC) -c $(FFLAGS64) $< -o $@

%.o: %.c
	$(CC) -c -w $< -o $@

# Test with standard DWBA input
test: ptolemy test_dwba.in
	@echo "=== Testing 16O(d,p)17O at 20 MeV ==="
	@./ptolemy < test_dwba.in > test_dwba.out 2>/dev/null; \
	if grep -q "TOTAL:" test_dwba.out; then \
		echo "✓ DCS computed successfully"; \
		grep "TOTAL:" test_dwba.out; \
	else \
		echo "✗ FAILED — no cross section output"; \
		echo "  Check test_dwba.out for errors"; \
		exit 1; \
	fi

clean:
	rm -f *.o ptolemy test_dwba.out

.PHONY: all clean test ptolemy64
