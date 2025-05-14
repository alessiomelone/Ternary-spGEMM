# Compiler and flags
CXX = g++
CXXFLAGS = -O2 -march=native -mtune=native -fstrict-aliasing -DNDEBUG -DPMU
LDFLAGS =

# Source files and output
SRCS = main.cpp perf.cpp
TARGET = SparseGEMM.out

# Default rule: build the target
all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)

# Clean rule: remove the target and any other build artifacts
clean:
	rm -f $(TARGET)

# Phony targets
.PHONY: all clean 