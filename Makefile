CC = g++
CFLAGS = -fpermissive -pedantic -Wall -Wextra -Wno-error=unused-variable -O2
OBJECTS = resketch.o

all: resketch
	rm -f $(OBJECTS)

resketch: resketch.o
	$(CC) -o resketch resketch.o $(CFLAGS) -lm

resketch.o: resketch.c ./lib/least-squares-cpp/include/lsqcpp.h
	$(CC) -c $(CFLAGS) -Ilib/least-squares-cpp/include/ -Ilib/eigen/ resketch.c

.PHONY: all clean
clean:
	rm -f $(OBJECTS) resketch