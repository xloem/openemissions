LDLIBS=-lpigpiod_if2 -lrt
CFLAGS=-Wall -Werror -pthread -std=c11 -pedantic

input: input.c impl*.c
