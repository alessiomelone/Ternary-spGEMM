# Compiler and flags
CXX = g++
# Base CXXFLAGS, these will always be applied
CXXFLAGS_BASE = -std=c++17 -O2 -march=native -mtune=native -fstrict-aliasing -DNDEBUG -DPMU -Wno-vla-cxx-extension

# Initialize CXXFLAGS with the base flags
CXXFLAGS = $(CXXFLAGS_BASE)
CXXFLAGS += -Icodegen

# Add CALIBRATE by default unless NO_CALIBRATE is specified
ifndef NO_CALIBRATE
	CXXFLAGS += -DCALIBRATE
	BUILD_INFO = (Calibrated Build)
else
	BUILD_INFO = (Non-Calibrated Build)
endif

# Instrumentation flag if requested
ifdef INSTRUMENT
	CXXFLAGS += -DINSTRUMENTATION_RUN
	BUILD_INFO := $(BUILD_INFO) (Instrumented)
endif

LDFLAGS =

# Source files and outputs
SRCS = cpp_impl/main.cpp \
       cpp_impl/perf.cpp

TARGET = sparseGEMM.out

# Default rule
all: $(TARGET)

$(TARGET): $(SRCS)
	@echo "--- Building $(TARGET) $(BUILD_INFO) ---"
	@echo "Using CXXFLAGS: $(CXXFLAGS)"
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)
	@echo "--- Finished building $(TARGET) ---"

# Clean rule
clean:
	@echo "--- Cleaning $(TARGET) and *.o ---"
	rm -f $(TARGET)
	rm -f *.o

.PHONY: all clean

# -----------------------------------------------------------------------------
# Codegen-only mode: minimal flags, include generated.h
CODEGEN_TARGET   = sparseGEMM_codegen.out
CODEGEN_FLAGS    = -std=c++17 -g -Icodegen

.PHONY: codegen
codegen: codegen/generated.h $(SRCS)
	@echo "--- Building $(CODEGEN_TARGET) (Codegen mode) ---"
	$(CXX) $(CODEGEN_FLAGS) -O0 $(SRCS) -o $(CODEGEN_TARGET)
	@echo "--- Finished building $(CODEGEN_TARGET) ---"
