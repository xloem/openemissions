CC=$(CXX)
LD=$(CXX)
CXXFLAGS=-std=c++14 -pthread $$(sdl2-config --cflags) -ggdb -Werror -Wall -Wpedantic
LDLIBS=-lpthread -lrtlsdr -litpp $$(sdl2-config --libs) -ggdb -Werror -Wall -Wpedantic

Main: Main.o RtlSdr.o Source.o SDL.o

clean:
	-rm *.o Main
