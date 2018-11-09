default: contour

contour: CXX = g++ -m64 -std=c++11
contour: CXXFLAGS = -I. -O3 -Wall -fopenmp -Wno-unknown-pragmas

contour: contours.o
	$(CXX) $(CXXFLAGS) -o contours contours.o

contours.o: contours.cpp
	$(CXX) $(CXXFLAGS) -c contours.cpp

clean:
	rm contours contours.o
