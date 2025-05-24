CXX = g++
CXXFLAGS = -std=c++11 -Wall

all: test_blocked_tcsc

test_blocked_tcsc: test_blocked_tcsc.cpp BlockedTCSC.h
	$(CXX) $(CXXFLAGS) -o test_blocked_tcsc test_blocked_tcsc.cpp

clean:
	rm -f test_blocked_tcsc 