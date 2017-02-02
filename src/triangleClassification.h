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
#include "nested3DGrid.h"
#include "3dGeometry.h"
#include <omp.h>



using namespace std;



void classifyTrianglesAndGenerateOutput(const Nested3DGridWrapper *uniformGrid, 
																				MeshIntersectionGeometry &geometry, 
                                        const unordered_set<const InputTriangle *> trianglesThatIntersect[2],
                                        vector<BoundaryPolygon> polygonsFromRetesselation[2],                                                                             
                                        ostream &outputStream);

#endif