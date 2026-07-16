CXX ?= g++
CXXFLAGS ?= -O3 -march=native -flto
CPPFLAGS ?=
WARNINGS := -Wall -Wextra -Wpedantic -Wconversion -Wshadow
LDLIBS := -lz -pthread
TARGETS := v4mer v4mer-jf

.PHONY: all clean

all: $(TARGETS)

v4mer: v4mer.cpp v4mer_core.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 $(WARNINGS) $< -o $@ $(LDLIBS)

v4mer-jf: v4mer_jf.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 $(WARNINGS) $< -o $@

clean:
	rm -f $(TARGETS)
