all: 1-1-prof-env test check pcmsummary

%: %.cpp *.hpp
	g++ -pedantic $*.cpp -o $@ $$(pkg-config --cflags eigen3) $$(pkg-config --libs fftw3) $$(pkg-config --libs SoapySDR) $$(root-config --cflags --libs) --std=c++14 -ggdb -O3

%.dbg: %.cpp *.hpp
	g++ -pedantic $*.cpp -o $@ $$(pkg-config --cflags eigen3) $$(pkg-config --libs fftw3) $$(pkg-config --libs SoapySDR) $$(root-config --cflags --libs) --std=c++14 -ggdb -O0 -DEIGEN_DONT_VECTORIZE
