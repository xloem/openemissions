CC=$(CXX)
LD=$(CXX)
CPPFLAGS=-DEIGEN_NO_DEBUG
CXXFLAGS=-std=c++14 -pthread $$(pkg-config --cflags eigen3) $$(sdl2-config --cflags) $$(pkg-config --cflags tinyxml2) $$(pkg-config --cflags flac++) -ggdb -Werror -Wall -Wpedantic -O0
#CXXFLAGS=-std=c++14 -pthread $$(pkg-config --cflags eigen3) $$(sdl2-config --cflags) $$(pkg-config --cflags tinyxml2) $$(pkg-config --cflags flac++) -ggdb         -Wall -Wpedantic -O3
LDLIBS=-lpthread $$(pkg-config --libs fftw3) -lrtlsdr $$(sdl2-config --libs) -lSDL2_ttf $$(pkg-config --libs tinyxml2) $$(pkg-config --libs flac++) -ggdb -Werror -Wall -Wpedantic -O0

Main: Main.o Source.o GUI.o RtlSdr.o SDL.o Oscilloscope.o Waterfall.o PeriodFinder.o PeriodViewer.o Math.o

Test: Test.o GUI.o SDL.o

run: Main
	./Main

runtest: Test
	./Test

clean:
	-rm *.o Main
