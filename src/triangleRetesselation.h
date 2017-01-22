#include <iostream>
#include <set>
#include "rationals.h"
#include "utils.h"
#include "nested3DGrid.h"
#include "3dGeometry.h"

using namespace std;

void retesselateIntersectingTriangles(MeshIntersectionGeometry & meshIntersectionGeometry, 
																			const vector< pair<VertexFromIntersection, 
																			VertexFromIntersection> >  &edgesFromIntersection, 
																			const vector< pair<InputTriangle *,InputTriangle *> > &intersectingTrianglesThatGeneratedEdges,
																			vector<BoundaryPolygon> polygonsFromRetesselation[2]);
