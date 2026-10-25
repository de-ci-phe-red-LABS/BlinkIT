# Build the two C evaluators into build/.
CC      ?= gcc
CFLAGS  ?= -O3 -march=native -fopenmp
BIN     := build/blink_cube build/blink_attack

all: $(BIN)

build/%: c/%.c
	@mkdir -p build
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -rf build

.PHONY: all clean
