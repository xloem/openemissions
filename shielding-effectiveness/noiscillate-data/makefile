all: test check pcmsummary

%: %.cpp *.hpp
	clang++ -pedantic $@.cpp -o $@ $$(pkg-config --cflags eigen3) $$(pkg-config --libs fftw3) $$(pkg-config --libs SoapySDR) -std=c++14 -ggdb -O3 #-DEIGEN_DONT_VECTORIZE
