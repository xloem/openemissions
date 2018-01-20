CC=$(CXX)
LD=$(CXX)
CXXFLAGS=-std=c++14 -pthread $$(pkg-config --cflags itpp) $$(sdl2-config --cflags) $$(pkg-config --cflags tinyxml2) $$(pkg-config --cflags flac++) -ggdb -Werror -Wall -Wpedantic -O0
LDLIBS=-lpthread $$(pkg-config --libs itpp) -lrtlsdr $$(sdl2-config --libs) -lSDL2_ttf $$(pkg-config --libs tinyxml2) $$(pkg-config --libs flac++) -ggdb -Werror -Wall -Wpedantic -O0

Main: Main.o Source.o GUI.o RtlSdr.o SDL.o Oscilloscope.o Waterfall.o PeriodFinder.o PeriodViewer.o

Test: Test.o GUI.o SDL.o

run: Main
	./Main

runtest: Test
	./Test

clean:
	-rm *.o Main
