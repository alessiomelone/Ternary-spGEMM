# Compiler and flags
CXX = g++
# Base CXXFLAGS, these will always be applied
CXXFLAGS_BASE = -std=c++17 -O2 -march=native -mtune=native -fstrict-aliasing -DNDEBUG -DPMU -Wno-vla-cxx-extension

# Initialize CXXFLAGS with the base flags
CXXFLAGS = $(CXXFLAGS_BASE)

# Add CALIBRATE by default unless NO_CALIBRATE is specified
ifndef NO_CALIBRATE
    CXXFLAGS += -DCALIBRATE
    BUILD_INFO = (Calibrated Build)
else
    BUILD_INFO = (Non-Calibrated Build)
endif

# This is where the conditional logic is added:
# Check if the 'INSTRUMENT' variable is defined when make is run
# Example: make INSTRUMENT=1 or make INSTRUMENT=true
ifdef INSTRUMENT
    # You could add more specific checks here if needed, e.g.:
    # ifeq ($(INSTRUMENT), true)
    # CXXFLAGS += -DINSTRUMENT_RUN
    # endif
    # For simplicity, if INSTRUMENT is defined to any non-empty value, we add the flag.
	CXXFLAGS += -DINSTRUMENTATION_RUN
	BUILD_INFO := $(BUILD_INFO) (Instrumented)
endif

LDFLAGS =

# Source files and output
SRCS = cpp_impl/main.cpp \
       cpp_impl/perf.cpp

TARGET = sparseGEMM.out

# Default rule: build the target
all: $(TARGET)

$(TARGET): $(SRCS)
	@echo "--- Building $(TARGET) $(BUILD_INFO) ---"
	@echo "Using CXXFLAGS: $(CXXFLAGS)"
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)
	@echo "--- Finished building $(TARGET) ---"

# Clean rule: remove the target and any other build artifacts
clean:
	@echo "--- Cleaning $(TARGET) and *.o ---"
	rm -f $(TARGET)
	rm -f *.o

# Phony targets
.PHONY: all clean