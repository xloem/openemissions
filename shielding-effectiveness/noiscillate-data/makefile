all: test check pcmsummary

%: %.cpp *.hpp
	g++ $@.cpp -o $@ $$(pkg-config --cflags eigen3) $$(pkg-config --libs fftw3) -std=c++14 -ggdb
