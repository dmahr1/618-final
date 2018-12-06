default: contour

contour: CXX = g++ -m64 -std=c++11
contour: CXXFLAGS = -I. -O3 -Wall -fopenmp -Wno-unknown-pragmas -g

contour: contours.o
	$(CXX) $(CXXFLAGS) -o contours contours.o -lgdal

contours.o: contours.cpp contours.h
	$(CXX) $(CXXFLAGS) -c contours.cpp

clean:
	rm contours contours.o
