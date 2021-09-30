CC=cc
MYCFLAGS=
CFLAGS=-Wall -Wextra -std=c99 -g -march=native -flto $(MYCFLAGS)
MYLDFLAGS= -Ofast
LDFLAGS=-lm -lgomp $(MYLDFLAGS)

TARGETS= main
SOURCES=$(shell echo *.c *.cpp)
OBJECTS=

all: $(TARGETS)

main: main.o
	$(CC) $(CFLAGS) -o tidal $^ $(LDFLAGS)

clean:
	rm -rf *.o tidal

.depend: $(SOURCES)
	$(CC) -MM $^ > $@

-include .depend

.PHONY: clean all