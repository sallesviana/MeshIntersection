flags = -lgmp -lgmpxx -std=c++11 -O3  -fopenmp -ltcmalloc

all: meshIntersection

meshIntersection: meshIntersection.o triangleClassification.o nested3DGrid.o 3d_objects.o utils.o
		g++ meshIntersection.o triangleClassification.o nested3DGrid.o 3d_objects.o utils.o $(flags) -o  meshIntersection


meshIntersection.o: src/meshIntersection.cpp
	g++ src/meshIntersection.cpp  $(flags) -c  -o meshIntersection.o

3d_objects.o: src/3d_objects.cpp src/3d_objects.h
	g++ src/3d_objects.cpp  $(flags) -c  -o 3d_objects.o

triangleClassification.o: src/triangleClassification.h src/triangleClassification.cpp
		g++ src/triangleClassification.cpp $(flags) -c -o triangleClassification.o


nested3DGrid.o: src/nested3DGrid.h src/nested3DGrid.cpp
		g++ src/nested3DGrid.cpp $(flags) -c -o nested3DGrid.o

utils.o: src/utils.h src/utils.cpp
		g++ src/utils.cpp $(flags) -c -o utils.o