all: meshIntersection

meshIntersection: src/PinMesh.cpp src/meshIntersection.cpp src/3d_objects.cpp src/floodFillScanline.cpp src/rationals.h  src/common2.h  src/nested3DGrid.cpp src/tritri_isectline.c
	g++ src/meshIntersection.cpp -lgmp -lgmpxx -std=c++11 -O3  -fopenmp  -o meshIntersection


