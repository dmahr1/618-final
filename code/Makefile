default: contour

contour: GDAL_INCLUDE = /afs/andrew.cmu.edu/usr16/dmahr/public/gdal/include
contour: GDAL_LIB = /afs/andrew.cmu.edu/usr16/dmahr/public/gdal/lib
contour: CXX = g++ -m64 -std=c++11
contour: CXXFLAGS = -I$(GDAL_INCLUDE) -O3 -Wall -fopenmp -Wno-unknown-pragmas -g

contour: contours.o
	$(CXX) $(CXXFLAGS) -o contours contours.o -lgdal -L$(GDAL_LIB)

contours.o: contours.cpp contours.h
	$(CXX) $(CXXFLAGS) -c contours.cpp

clean:
	rm contours contours.o
