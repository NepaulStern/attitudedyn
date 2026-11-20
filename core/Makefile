CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Iinclude -Isrc

SRC = src/layout.cpp src/draws.cpp src/likelihood.cpp src/postestimation.cpp

.PHONY: test clean

test: tests/test_core
	./tests/test_core

tests/test_core: tests/test_core.cpp $(SRC) include/attitudedyn/api.hpp include/attitudedyn/normal.hpp
	$(CXX) $(CXXFLAGS) -o $@ tests/test_core.cpp $(SRC)

clean:
	rm -f tests/test_core
