CXX      = g++
CXXFLAGS = -O2 -std=c++17 -Iinclude
LDFLAGS  = -lm

TARGET   = ptolemy_cpp

SRCS = src/main.cpp \
       src/elastic/elastic.cpp \
       src/dwba/a12.cpp \
       src/dwba/av18_potential.cpp \
       src/dwba/bound.cpp \
       src/dwba/coulin.cpp \
       src/dwba/dwba.cpp \
       src/dwba/grdset_ineldc_faithful.cpp \
       src/dwba/ineldc_collective.cpp \
       src/dwba/math_utils.cpp \
       src/dwba/potential_eval.cpp \
       src/dwba/rcwfn.cpp \
       src/dwba/setup.cpp \
       src/dwba/spline.cpp \
       src/dwba/stubs.cpp \
       src/dwba/thiele_cf.cpp \
       src/dwba/wavelj.cpp \
       src/dwba/xsectn.cpp \
       src/input/InputGenerator.cpp \
       src/input/Isotope.cpp \
       src/input/Potentials.cpp \
       src/input/PtolemyParser.cpp

OBJS = $(SRCS:.cpp=.o)

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(TARGET) $(OBJS)

test: $(TARGET)
	@echo "=== 206Hg(d,p) transfer ==="
	./$(TARGET) < test_inputs/test_hg206dp_gs.in | tail -5
	@echo "=== Done ==="
