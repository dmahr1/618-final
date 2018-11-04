default: contour

contour: CXXFLAGS = -I. -Wall -fopenmp -Wno-unknown-pragmas

contour: contours.o
	g++ -std=c++11 -o contours contours.o

contours.o: contours.cpp
	g++ -std=c++11 -c contours.cpp
