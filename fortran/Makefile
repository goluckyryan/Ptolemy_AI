# Makefile for Ptolemy DWBA Code (Fortran)
#
# Tested on Ubuntu 24.04 with gfortran 13/14
#
# CRITICAL: Use -O0 on 64-bit. The code relies on Fortran 77 behaviors
# (uninitialized variables, implicit SAVE, pointer-as-integer) that
# break with optimization. The 32-bit build tolerates -O2.
#
#   make             — auto-build (64-bit -O0)
#   make ptolemy32   — 32-bit build (requires gcc-multilib, can use -O2)
#   make test        — run test and verify DCS output
#   make clean

FC      = gfortran
CC      = gcc

# Fortran flags for legacy code compatibility
# -O0 is REQUIRED on 64-bit (code breaks at -O1+)
# -std=legacy allows old F77 constructs and dummy array size mismatches
FFLAGS  = -O0 -std=legacy -fallow-argument-mismatch -fno-range-check \
          -fallow-invalid-boz -fd-lines-as-comments -w

# Source files
FSRC    = ptolemy-main.f source.f fortlib.f gfortran_stuff.f masstable.f \
          numbered_store.f linkule.f linkulesfitters.f keep.f keepsub.f \
          keptsub.f phiffer.f av18.f
CSRC    = dtime.c srread.c

FOBJ    = $(FSRC:.f=.o)
COBJ    = $(CSRC:.c=.o)

# Default: 64-bit build with -O0
all: ptolemy

ptolemy: $(FOBJ) $(COBJ)
	@echo "Linking..."
	$(FC) -o $@ $(FOBJ) $(COBJ)
	@echo "=== Build complete: ./ptolemy ==="

# Pattern rules
%.o: %.f
	$(FC) -c $(FFLAGS) $< -o $@

%.o: %.c
	$(CC) -c -w $< -o $@

# 32-bit build (can use -O2, matches original 32-bit assumptions)
ptolemy32:
	@if $(CC) -m32 -xc -c /dev/null -o /dev/null 2>/dev/null; then \
		echo "=== Building 32-bit ==="; \
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

.PHONY: all clean test ptolemy32
