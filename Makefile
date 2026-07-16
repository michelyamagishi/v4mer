CXX ?= g++
CXXFLAGS ?= -O3 -march=native -flto
CPPFLAGS ?=
WARNINGS := -Wall -Wextra -Wpedantic -Wconversion -Wshadow
LDLIBS := -lz -pthread
TARGET := v4mer

.PHONY: all clean

all: $(TARGET)

$(TARGET): v4mer.cpp v4mer_core.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++17 $(WARNINGS) $< -o $@ $(LDLIBS)

clean:
	rm -f $(TARGET)
