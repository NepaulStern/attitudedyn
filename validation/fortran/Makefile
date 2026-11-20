# One GSS likelihood evaluation of the original program, then the C++ comparison.
# GQOPT and the IMSL-style routines are linked only into this oracle. They are
# not part of the R package.

FC = gfortran
FFLAGS = -O1 -std=legacy -ffixed-form -ffixed-line-length-none -fopenmp
CXX = g++
CXXFLAGS = -std=c++17 -O2 -I../../src/core/include

FORTRAN_SRC = src/opt.f src/nmsimp.f src/lnsr.f src/matp.f src/eigrs.f src/bdummy.f

.PHONY: all clean run

all: build/alecia_oracle build/compare_gss

build/alecia_oracle.f build/alecia.p: src/alecia_parallel.f src/alecia.p.gss src/gss_log.txt apply_patch.py
	python3 apply_patch.py

build/alecia_oracle: build/alecia_oracle.f oracle_extra.f shared_math.cpp $(FORTRAN_SRC)
	$(FC) $(FFLAGS) -c build/alecia_oracle.f -o build/alecia_oracle.o
	$(FC) $(FFLAGS) -c oracle_extra.f -o build/oracle_extra.o
	$(CXX) $(CXXFLAGS) -c shared_math.cpp -o build/shared_math.o
	$(FC) $(FFLAGS) -c src/opt.f -o build/opt.o
	$(FC) $(FFLAGS) -c src/nmsimp.f -o build/nmsimp.o
	$(FC) $(FFLAGS) -c src/lnsr.f -o build/lnsr.o
	$(FC) $(FFLAGS) -c src/matp.f -o build/matp.o
	$(FC) $(FFLAGS) -c src/eigrs.f -o build/eigrs.o
	$(FC) $(FFLAGS) -c src/bdummy.f -o build/bdummy.o
	$(FC) $(FFLAGS) -o build/alecia_oracle build/alecia_oracle.o build/oracle_extra.o \
	  build/shared_math.o build/opt.o build/nmsimp.o build/lnsr.o build/matp.o \
	  build/eigrs.o build/bdummy.o -lstdc++

build/compare_gss: ../compare_gss.cpp ../../src/core/src/layout.cpp ../../src/core/src/likelihood.cpp \
	../../src/core/src/draws.cpp ../../src/core/src/postestimation.cpp
	$(CXX) $(CXXFLAGS) -o build/compare_gss ../compare_gss.cpp ../../src/core/src/layout.cpp \
	  ../../src/core/src/likelihood.cpp ../../src/core/src/draws.cpp ../../src/core/src/postestimation.cpp -pthread

clean:
	rm -rf build
