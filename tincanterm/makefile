LDLIBS=-lpigpiod_if2 -lpigpio -lrt
CPPFLAGS=-I src/tincanlib
CXXFLAGS=-Wall -Werror -pthread -pedantic

%.cpp: %.h src/tincanlib/config.h

tincanterm: src/tincanlib/*.cpp src/tincanterm/*.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDLIBS) $^ -o $@

tincandiag: src/tincanlib/*.cpp src/tincandiag/*.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDLIBS) $^ -o $@
