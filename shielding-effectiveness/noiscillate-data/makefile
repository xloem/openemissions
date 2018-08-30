test: test.cpp *.hpp
	clang++ test.cpp -o test $$(pkg-config --cflags eigen3) $$(pkg-config --libs fftw3) -std=c++14 -ggdb
