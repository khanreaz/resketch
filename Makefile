CXX = g++
CXXFLAGS = -fpermissive -pedantic -Wall -Wextra -Wno-error=unused-variable -O2 -std=c++11
OBJECTS = resketch.o

all: resketch
	rm -f $(OBJECTS)

resketch: resketch.o
	$(CXX) -o resketch resketch.o $(CXXFLAGS) -lm

resketch.o: resketch.c ./lib/least-squares-cpp/include/lsqcpp.h
	$(CXX) -c $(CXXFLAGS) -Ilib/least-squares-cpp/include/ -Ilib/eigen/ resketch.c

.PHONY: all clean
clean:
	rm -f $(OBJECTS) resketch
