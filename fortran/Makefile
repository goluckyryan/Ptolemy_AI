# Makefile for Ptolemy DWBA Code (Fortran)
#
# Tested on Ubuntu 24.04 with gfortran 13/14
#
# The original code was written for 32-bit systems. On 64-bit systems,
# the internal memory allocator (GRAB/LOC) stores pointer offsets in
# INTEGER variables that may overflow. Two build modes are provided:
#
#   make ptolemy32   — 32-bit build (requires gcc-multilib), most reliable
#   make ptolemy64   — 64-bit build (may have pointer arithmetic issues)
#   make             — auto-selects 32-bit if multilib available, else 64-bit

FC      = gfortran
CC      = gcc

# Fortran flags for legacy code compatibility
FFLAGS  = -O2 -std=legacy -fallow-argument-mismatch -fno-range-check \
          -fallow-invalid-boz -fd-lines-as-comments -fno-automatic -w

# 32-bit flags (append to FFLAGS and CFLAGS)
M32     = -m32

# Source files
FSRC    = ptolemy-main.f source.f fortlib.f gfortran_stuff.f masstable.f \
          numbered_store.f linkule.f linkulesfitters.f keep.f keepsub.f \
          keptsub.f phiffer.f av18.f
CSRC_32 = dtime.c srread.c
CSRC_64 = dtime.c srread_64.c

FOBJ    = $(FSRC:.f=.o)
COBJ_32 = $(CSRC_32:.c=.o)
COBJ_64 = dtime.o srread_64.o

# Default: try 32-bit, fall back to 64-bit
.PHONY: all clean check32

all: check32

check32:
	@if $(CC) -m32 -xc -c /dev/null -o /dev/null 2>/dev/null; then \
		echo "=== Building 32-bit (recommended) ==="; \
		$(MAKE) ptolemy32; \
	else \
		echo "=== 32-bit libs not found, building 64-bit ==="; \
		echo "    Install gcc-multilib for 32-bit build:"; \
		echo "    sudo apt install gcc-multilib gfortran-multilib"; \
		$(MAKE) ptolemy64; \
	fi

# 32-bit build (recommended — matches original Ptolemy assumptions)
ptolemy32: FFLAGS += $(M32)
ptolemy32: $(FSRC) $(CSRC_32)
	@echo "Compiling Fortran sources (32-bit)..."
	@for f in $(FSRC); do \
		echo "  $$f"; \
		$(FC) -c $(FFLAGS) $(M32) $$f || exit 1; \
	done
	@echo "Compiling C sources (32-bit)..."
	$(CC) -c $(M32) -w dtime.c
	$(CC) -c $(M32) -w srread.c
	@echo "Linking..."
	$(FC) $(M32) -o ptolemy $(FOBJ) $(COBJ_32)
	@echo "=== Build complete: ./ptolemy (32-bit) ==="

# 64-bit build (uses srread_64.c with intptr_t fix)
ptolemy64: $(FSRC) $(CSRC_64)
	@echo "Compiling Fortran sources (64-bit)..."
	@for f in $(FSRC); do \
		echo "  $$f"; \
		$(FC) -c $(FFLAGS) $$f || exit 1; \
	done
	@echo "Compiling C sources (64-bit, with pointer fix)..."
	$(CC) -c -w dtime.c
	$(CC) -c -w srread_64.c
	@echo "Linking..."
	$(FC) -o ptolemy $(FOBJ) $(COBJ_64)
	@echo "=== Build complete: ./ptolemy (64-bit) ==="
	@echo "⚠️  64-bit build may have pointer arithmetic issues."
	@echo "   Install gcc-multilib for reliable 32-bit build."

# Test with a simple input
test: ptolemy
	@echo "Running test01.in..."
	@./ptolemy < test01.in > test01.out 2>/dev/null
	@if grep -q "TOTAL:" test01.out 2>/dev/null; then \
		echo "✓ Test passed — cross section computed"; \
		grep "TOTAL:" test01.out; \
	else \
		echo "✗ Test failed — no cross section output"; \
	fi

clean:
	rm -f *.o ptolemy test01.out

.PHONY: all clean test check32 ptolemy32 ptolemy64
