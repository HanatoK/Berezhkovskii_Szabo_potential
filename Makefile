CXX=g++
SWIG=swig

all: main library

main: main.cpp cv.cpp cv.h
	$(CXX) main.cpp cv.cpp -Wall -std=c++20 -O2 -o main

library: cv.h
	$(SWIG) -c++ -Wall -tcl cv.i
	$(CXX) -shared cv_wrap.cxx cv.cpp -std=c++20 -fPIC -o cv.so -O2

clean:
	rm -f cv_wrap.cxx main cv.so
