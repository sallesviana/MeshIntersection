#ifndef TRIANGLE_CLASSIFICATION_H
#define TRIANGLE_CLASSIFICATION_H

#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <cstdlib>
#include <unordered_set>
#include <unordered_map>
#include <time.h>
#include "rationals.h"
#include "utils.h"
#include "3d_objects.h"
#include "nested3DGrid.h"
#include <omp.h>



using namespace std;

void classifyTrianglesAndGenerateOutput(const Nested3DGridWrapper *uniformGrid, 
																				const unordered_set<const Triangle *> trianglesThatIntersect[2],
																				vector<BoundaryPolygon> polygonsFromRetesselation[2],
																				const vector<Point> vertices[3],
																				const vector<Triangle> triangles[2],
																				ostream &outputStream);

#endif