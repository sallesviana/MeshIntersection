flags = -lgmp -lgmpxx -std=c++11 -O3  -fopenmp -ltcmalloc 

all: meshIntersection

meshIntersection: pinMesh.o floodFillScanline.o meshIntersection.o nested3DGrid.o 3dGeometry.o utils.o triangleRetesselation.o triangleClassification.o
		g++ pinMesh.o floodFillScanline.o meshIntersection.o nested3DGrid.o 3dGeometry.o triangleRetesselation.o triangleClassification.o utils.o $(flags) -o  meshIntersection


meshIntersection.o: src/meshIntersection.cpp
	g++ src/meshIntersection.cpp  $(flags) -c  -o meshIntersection.o

3dGeometry.o: src/3dGeometry.cpp src/3dGeometry.h src/3dGeometryGeometricalPart.cpp src/tritri_isectline.c  src/boundaryPolygon.cpp
	g++ src/3dGeometry.cpp  $(flags) -c  -o 3dGeometry.o


triangleRetesselation.o: src/triangleRetesselation.h src/triangleRetesselation.cpp
	g++ src/triangleRetesselation.cpp  $(flags) -c  -o triangleRetesselation.o

triangleClassification.o: src/triangleClassification.h src/triangleClassification.cpp
		g++ src/triangleClassification.cpp $(flags) -c -o triangleClassification.o


nested3DGrid.o: src/nested3DGrid.h src/nested3DGrid.cpp
		g++ src/nested3DGrid.cpp $(flags) -c -o nested3DGrid.o

utils.o: src/utils.h src/utils.cpp
		g++ src/utils.cpp $(flags) -c -o utils.o

pinMesh.o: src/PinMesh.h src/PinMesh.cpp src/3dGeometry.h src/nested3DGrid.h 
		g++ src/PinMesh.cpp $(flags) -c -o pinMesh.o

floodFillScanline.o: src/floodFillScanline.h src/floodFillScanline.cpp
		g++ src/floodFillScanline.cpp $(flags) -c -o floodFillScanline.o