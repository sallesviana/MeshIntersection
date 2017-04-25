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


//returns the time to classify (without the time to write the output)
double classifyTrianglesAndGenerateOutput(const Nested3DGridWrapper *uniformGrid, 
																				MeshIntersectionGeometry &geometry, 
                                        const unordered_set<const InputTriangle *> trianglesThatIntersect[2],
                                        vector< pair<const InputTriangle *,vector<BoundaryPolygon>> > polygonsFromRetesselationOfEachTriangle[2],                                                                             
                                        ostream &outputStream);

#endif