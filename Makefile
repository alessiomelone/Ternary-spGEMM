# Compiler and flags
CXX = g++
CXXFLAGS = -O2 -march=native -mtune=native -fstrict-aliasing -DNDEBUG -DPMU
LDFLAGS =

# Source files and output
SRCS = cpp_impl/main.cpp \
       cpp_impl/perf.cpp

TARGET = sparseGEMM

# Default rule: build the target
all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)

# Clean rule: remove the target and any other build artifacts
clean:
	rm -f $(TARGET)
	rm -f *.o

# Phony targets
.PHONY: all clean 