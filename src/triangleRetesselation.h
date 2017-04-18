#include <iostream>
#include <set>
#include <sstream>
#include "rationals.h"
#include "utils.h"
#include "nested3DGrid.h"
#include "3dGeometry.h"

using namespace std;

#define COLLECT_STATISTICS

void retesselateIntersectingTriangles(MeshIntersectionGeometry & meshIntersectionGeometry, 
																			const vector< pair<VertexFromIntersection, 
																			VertexFromIntersection> >  &edgesFromIntersection, 
																			const vector< pair<InputTriangle *,InputTriangle *> > &intersectingTrianglesThatGeneratedEdges,
																			vector< pair<const InputTriangle *,vector<BoundaryPolygon>> > polygonsFromRetesselationOfEachTriangle[2]);
