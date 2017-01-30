CC=$(CXX)
LD=$(CXX)
CXXFLAGS=-std=c++14 $$(sdl2-config --cflags) -ggdb -Werror -Wall -Wpedantic
LDLIBS=-lrtlsdr -litpp $$(sdl2-config --libs) -ggdb -Werror -Wall -Wpedantic

Main: Main.o RtlSdr.o Source.o SDL.o

clean:
	-rm *.o Main
