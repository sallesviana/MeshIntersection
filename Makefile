all: 3dIntersection

3dIntersection: src/3dIntersection.cpp src/3d_objects.cpp src/floodFillScanline.cpp src/rationals.h  src/common2.h  src/nested3DGrid.cpp src/tritri_isectline.c
	g++ src/3dIntersection.cpp -lgmp -lgmpxx -std=c++11 -O3  -fopenmp  -o 3dIntersection


