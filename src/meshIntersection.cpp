/*
Copyright 2016 Salles V. G. Magalhaes, W. R. Franklin, Marcus Andrade

This file is part of PinMesh.

PinMesh is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PinMesh is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PinMesh.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <cstdlib>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <queue>
#include <sstream>
#include <time.h>
#include "3d_objects.h"
#include "rationals.h"
#include "utils.h"
#include "nested3DGrid.h"
#include "triangleClassification.h"
#include "common2.h"
#include <omp.h>



using namespace std;


//===============================================================
// Constants, flags...
//===============================================================

//#define DEBUGGING_MODE

//This may slow down the algorithm
#define COLLECT_STATISTICS
//#define COLLECT_STATISTICS_PRINT_TRIANGLES_INTERSECTIONS

//#define SANITY_CHECKS

//===============================================================





double timeReadData,timeCreateGrid,timeDetectIntersections,timeRetesselate,timeClassifyTriangles;




#include "tritri_isectline.c"


//Each vector represents the vertices of a layer
//The first layear is mesh 0, the second is mesh1 and the third contains vertices generated from the intersection of the two meshes...
vector<Point> vertices[3]; 

//Each vector represents a set of objects in the same layer
//The objects are represented by a set of triangles (defining their boundaries)
vector<Triangle> triangles[2]; 
//we don't need bounding-boxes for the output triangles...
vector<TriangleNoBB> trianglesFromRetesselation[2]; //triangles formed by retesselating triangles from each mesh...

vector<BoundaryPolygon> polygonsFromRetesselation[2]; //polygons generated after retesselation of triangles from each mesh..


Point bBox[2][2] ;//bounding box of each mesh (each boundingbox has two vertices)
Point boundingBoxTwoMeshesTogetter[2]; //bounding box considering both meshes togetter

Nested3DGridWrapper uniformGrid;

timespec t0BeginProgram, t0AfterDatasetRead;

//if vertexId negative, the point is shared...
Point * getPointFromVertexId(int vertexId, int meshId) {
 if(vertexId>=0) assert(vertexId < vertices[meshId].size()); //TODO: remover assertions...
 else assert(-vertexId -1 < vertices[2].size());

 return (vertexId>=0)?&vertices[meshId][vertexId]:&vertices[2][-vertexId -1];
}

bool isVertexSharedBetweenMeshes(int vertexId) {
  return vertexId<0;
}

//=======================================================================================================================




#include "debuggingUtils.cpp"




// TODO: perform ONLY necessary computation..
//TODO: example we do not need z to check if the z component of the normal is 0...
//Given a triangle, returns a plane (X=0,Y=0 or Z=0) such that the triangle is not perpendicular to that plane..
//A triangle is perpendicular to a plane iff its normal has the corresponding coordinate equal to 0
//a = p2-p1
//b = p3-p1
//normal: (a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x)
//tempVars should have at least 8 elements...
int getPlaneTriangleIsNotPerpendicular(const Point &p1, const Point &p2,const Point &p3, VertCoord tempVars[]) {
  //a = (tempVars[0],tempVars[1],tempVars[2])
  tempVars[0]=p2[0];
  tempVars[0]-=p1[0];

  tempVars[1]=p2[1];
  tempVars[1]-=p1[1];

  
  //b = (tempVars[3],tempVars[4],tempVars[5])
  tempVars[3]=p3[0]; 
  tempVars[3]-=p1[0];

  tempVars[4]=p3[1];
  tempVars[4]-=p1[1];

  

  //Let's check if the z-component of the normal is 0
  //if it is --> the triangle is perpendicular to z=0...
  tempVars[6] = tempVars[0]; //a.x
  tempVars[6] *= tempVars[4]; //b.y
  tempVars[7] = tempVars[1]; //a.y
  tempVars[7] *=tempVars[3]; //b.x
  tempVars[6] -= tempVars[7]; //a.x*b.y - a.y*b.x
  if(sgn(tempVars[6])!=0) return PLANE_Z0; //the triangle is NOT perpendicular to z=0...


  tempVars[2]=p2[2];
  tempVars[2]-=p1[2];

  tempVars[5]=p3[2];
  tempVars[5]-=p1[2];

  tempVars[6] = tempVars[2]; //a.z
  tempVars[6] *= tempVars[3]; //b.x
  tempVars[7] = tempVars[0]; //a.x
  tempVars[7] *=tempVars[5]; //b.z
  tempVars[6] -= tempVars[7]; //a.z*b.x - a.x*b.z
  if(sgn(tempVars[6])!=0) return PLANE_Y0; //the triangle is NOT perpendicular to y=0...
  
  return PLANE_X0; //we do not need to check... (unless the triangle is just a point...)

}




void extractPairsTrianglesInGridCell(const Nested3DGrid *grid,int i,int j, int k, int gridSize,vector<pair<Triangle *,Triangle *> > &pairsTrianglesToProcess) {
	int numTrianglesMesh0 = grid->numTrianglesInGridCell(0, gridSize,i,j,k);//cell.triangles[0].size();
  int numTrianglesMesh1 = grid->numTrianglesInGridCell(1, gridSize,i,j,k);



  Triangle **ptrTriMesh0 = grid->getPointerStartListTriangles(0,gridSize,i,j,k);
  Triangle **ptrTriMesh1Temp = grid->getPointerStartListTriangles(1,gridSize,i,j,k);

  for(int tA = 0;tA < numTrianglesMesh0; tA++) {
  	Triangle **ptrTriMesh1 = ptrTriMesh1Temp;
    for(int tB=0;tB<numTrianglesMesh1;tB++) {
      pairsTrianglesToProcess.push_back(pair<Triangle *,Triangle *>(*ptrTriMesh0, *ptrTriMesh1));
      ptrTriMesh1++;
    }
    ptrTriMesh0++;
  }

}


unsigned long long sumPairsTrianglesInGridCell(const Nested3DGrid *grid,int i,int j, int k, int gridSize,vector<pair<Triangle *,Triangle *> > &pairsTrianglesToProcess) {
	int numTrianglesMesh0 = grid->numTrianglesInGridCell(0, gridSize,i,j,k);//cell.triangles[0].size();
  int numTrianglesMesh1 = grid->numTrianglesInGridCell(1, gridSize,i,j,k);

  return ((unsigned long long)numTrianglesMesh0)*numTrianglesMesh1;

}

unsigned long long computeNumPairsTrianglesToProcessBeforeUnique(const Nested3DGridWrapper *uniformGrid,vector<pair<Triangle *,Triangle *> > &pairsTrianglesToProcess) {
	int gridSizeLevel1 =  uniformGrid->gridSizeLevel1;
	int gridSizeLevel2 =  uniformGrid->gridSizeLevel2;

	unsigned long long ans = 0;
	for(int i=0;i<gridSizeLevel1;i++) 
    for(int j=0;j<gridSizeLevel1;j++) 
      for(int k=0;k<gridSizeLevel1;k++) {
        if (uniformGrid->grid.hasSecondLevel(i,j,k) ) {
        	Nested3DGrid *secondLevelGrid = uniformGrid->grid.getChildGrid(i,j,k); 
			    for(int iLevel2=0;iLevel2<gridSizeLevel2;iLevel2++) 
			    	for(int jLevel2=0;jLevel2<gridSizeLevel2;jLevel2++) 
			      	for(int kLevel2=0;kLevel2<gridSizeLevel2;kLevel2++) {
			      		  ans += sumPairsTrianglesInGridCell(secondLevelGrid,iLevel2,jLevel2,kLevel2,gridSizeLevel2,pairsTrianglesToProcess);   
			      	}    
        } else {
        	ans += sumPairsTrianglesInGridCell(&(uniformGrid->grid),i,j,k,gridSizeLevel1,pairsTrianglesToProcess);         
        }
      }
  return ans;
}

//the vector will be filled with pointers to the triangles in the same uniform grid cells
//pairs will appear maximum once in the vector
//the first element in the pair is a triangle from mesh 0 and the second one is from mesh 1
void getPairsTrianglesInSameUnifGridCells(const Nested3DGridWrapper *uniformGrid,vector<pair<Triangle *,Triangle *> > &pairsTrianglesToProcess) {
	//timespec t0,t1;
	//clock_gettime(CLOCK_REALTIME, &t0);

	//TODO: remove this for performance purposess.... just for debugging
	//cerr << "Num pairs to process according to this configuration: " << computeNumPairsTrianglesToProcessBeforeUnique(uniformGrid,pairsTrianglesToProcess) << "\n";
	//exit(0);

	pairsTrianglesToProcess.reserve(min(uniformGrid->trianglesInGrid[0]->size(),uniformGrid->trianglesInGrid[1]->size()));

  int gridSizeLevel1 =  uniformGrid->gridSizeLevel1;
	int gridSizeLevel2 =  uniformGrid->gridSizeLevel2;



  for(int i=0;i<gridSizeLevel1;i++) 
    for(int j=0;j<gridSizeLevel1;j++) 
      for(int k=0;k<gridSizeLevel1;k++) {
        if (uniformGrid->grid.hasSecondLevel(i,j,k) ) {
        	Nested3DGrid *secondLevelGrid = uniformGrid->grid.getChildGrid(i,j,k); 
			    for(int iLevel2=0;iLevel2<gridSizeLevel2;iLevel2++) 
			    	for(int jLevel2=0;jLevel2<gridSizeLevel2;jLevel2++) 
			      	for(int kLevel2=0;kLevel2<gridSizeLevel2;kLevel2++) {
			      		  extractPairsTrianglesInGridCell(secondLevelGrid,iLevel2,jLevel2,kLevel2,gridSizeLevel2,pairsTrianglesToProcess);   
			      	}    
        } else {
        	extractPairsTrianglesInGridCell(&(uniformGrid->grid),i,j,k,gridSizeLevel1,pairsTrianglesToProcess);         
        }
      }

//  clock_gettime(CLOCK_REALTIME, &t1);
  //cerr << "T before sort: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 

  cerr << "Pairs before unique: " << pairsTrianglesToProcess.size() << "\n";
  sort(pairsTrianglesToProcess.begin(),pairsTrianglesToProcess.end());
  vector<pair<Triangle *,Triangle *> >::iterator it = std::unique (pairsTrianglesToProcess.begin(), pairsTrianglesToProcess.end());
  pairsTrianglesToProcess.resize( std::distance(pairsTrianglesToProcess.begin(),it) ); // 10 20 30 20 10 

  cerr << "Pairs after unique: " << pairsTrianglesToProcess.size() << "\n";

 // clock_gettime(CLOCK_REALTIME, &t1);
  //cerr << "T after sort: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n";    
}



/*
void saveEdges(const string &path) {
  ofstream fout(path.c_str());

  int numVert = 3*edges.size();
  int numEdges = 3*edges.size();
  int numFaces = edges.size();
  fout << numVert << " " << numEdges << " " << numFaces << "\n";
  for(pair< array<VertCoord,3>,array<VertCoord,3> > e:edges) {
    fout << e.first[0].get_d() << " " << e.first[1].get_d() << " " << e.first[2].get_d() << "\n";
    fout << e.second[0].get_d() << " " << e.second[1].get_d() << " " << e.second[2].get_d() << "\n";
    fout << e.second[0].get_d() << " " << e.second[1].get_d() << " " << e.second[2].get_d() << "\n";
  }
  fout << endl;

  int start = 1;
  for(int i=0;i<numFaces;i++) {
    fout << start << " " << start+1 << "\n";
    fout << start+1 << " " << start+2 << "\n";
    fout << start+2 << " " << start << "\n";
    start+= 3;
  }
  fout << endl;
  start = 1;
  for(int i=0;i<numFaces;i++) {
    fout << start << " " << start+1 << " " << start+2 <<  "\n";
    start+= 3;
  }
}
*/



//--------------------------------------------------------------
//Source: http://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
//tempVars should have at least 3 elements...
//
int orientation(const Point &p, const Point &q, const Point &r, int whatPlaneProjectTo,  VertCoord tempVars[])
{
    /*VertCoord val = (q[1] - p[1]) * (r[0] - q[0]) -
              (q[0] - p[0]) * (r[1] - q[1]);
 
    if (val == 0) return 0;  // colinear
 
    return (val > 0)? 1: 2; // clock or counterclock wise*/

    //we will project the points to the plane whatPlaneProjectTriangleTo
    //what coordinates should we check during computations?
    //if the points are projected to z=0 --> we have to use x and y
    //if the points are projected to y=0 --> we have to use x and z
    //...............................x=0 --> we have to use y and z
    int coord1=0,coord2=1;
    if(whatPlaneProjectTo==PLANE_Y0) {
      coord1=2;
      coord2=0;
    } else if(whatPlaneProjectTo==PLANE_X0) {
      coord1=1;
      coord2=2;
    }

    tempVars[0] = q[coord2];
    tempVars[0] -= p[coord2];
    tempVars[1] = r[coord1];
    tempVars[1] -= q[coord1];
    tempVars[0]*=tempVars[1];//tempVars[0] = (q[1] - p[1]) * (r[0] - q[0])

    tempVars[2] = q[coord1];
    tempVars[2] -=p[coord1];
    tempVars[1] = r[coord2];
    tempVars[1] -= q[coord2];
    tempVars[1] *= tempVars[2]; //tempVars[1] =  (q[0] - p[0]) * (r[1] - q[1])

   /* tempVars[0] -= tempVars[1]; 
    int sign = sgn(tempVars[0]);
    if(sign==0) return 0;
    return sign>0?1:2; */
    if(tempVars[0]==tempVars[1]) return 0;
    if(tempVars[0]>tempVars[1]) return 1;
    return 2;
}

//TODO : reduce amount of tests...
//TODO: maybe use address of vertices to reduce amount of tests...
//TODO: avoid memory copy... (remove max, min, etc);
// Given three colinear points p, q, r, the function checks if
// point q lies on line segment 'pr'
//All segments are co-planar (obvious) and are projected to plane whatPlaneProjectTriangleTo
bool onSegment(const Point & p, const Point & q, const Point & r, int whatPlaneProjectTo)
{
		//if (&q==&r || &q==&p) return false; //the endpoints coincide...
    //we will project the points to the plane whatPlaneProjectTriangleTo
    //what coordinates should we check during computations?
    //if the points are projected to z=0 --> we have to use x and y
    //if the points are projected to y=0 --> we have to use x and z
    //...............................x=0 --> we have to use y and z
    int coord1=0,coord2=1;
    if(whatPlaneProjectTo==PLANE_Y0) {
      coord1=2;
      coord2=0;
    } else if(whatPlaneProjectTo==PLANE_X0) {
      coord1=1;
      coord2=2;
    }

    if (q[coord1] < max(p[coord1], r[coord1]) && q[coord1] > min(p[coord1], r[coord1]) &&
        q[coord2] <= max(p[coord2], r[coord2]) && q[coord2] >= min(p[coord2], r[coord2]))
       return true;

    if (q[0] <= max(p[coord1], r[coord1]) && q[coord1] >= min(p[coord1], r[coord1]) &&
        q[coord2] < max(p[coord2], r[coord2]) && q[coord2] > min(p[coord2], r[coord2])) 
    	return true;
 
    return false;
}
 
// TODO: special cases!!!
//do p1-q1 intersect p2-q2 ?
//the two segments are co-planar
//during computation, we project them to whatPlaneProjectTriangleTo plane... 
bool doIntersect(const Point &p1, const Point &q1, const Point &p2, const Point &q2, int whatPlaneProjectTriangleTo, VertCoord tempVars[])
{
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(p1, q1, p2,whatPlaneProjectTriangleTo,tempVars);

    // Special Cases
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1,whatPlaneProjectTriangleTo)) return true;


    int o2 = orientation(p1, q1, q2,whatPlaneProjectTriangleTo,tempVars);

    // p1, q1 and p2 are colinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1,whatPlaneProjectTriangleTo)) return true;

    //the two segments are collinear, but they do not intersect 
    if(o1==0 && o2==0) return false;


    int o3 = orientation(p2, q2, p1,whatPlaneProjectTriangleTo,tempVars);

    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2,whatPlaneProjectTriangleTo)) return true;

    int o4 = orientation(p2, q2, q1,whatPlaneProjectTriangleTo,tempVars); 
    
 
     // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2,whatPlaneProjectTriangleTo)) return true;

    // General case
    if (o1!=0 && o2!=0 && o3!=0  && o4!=0)
      if (o1 != o2 && o3 != o4) {
    	 return true;
      }

 
    return false; // Doesn't fall in any of the above cases
}



//returns true iff the candidate edge (represented using the ids of the vertices in meshIdToProcess -- negative vertices are in the "common layer")
//intersects an edge from the set edgesToTest (these edges are in the same layer).
//we do not consider intersections in endpoints..
//TODO: consider special cases: vertical triangle, parallel edges, etc...
long long ctEdgeIntersect  = 0; //how many of the intersection tests performed are true?
long long ctEdgeDoNotIntersect = 0;
bool intersects(const pair<int,int> &candidateEdge,const set<pair<int,int> > &edgesToTest,int meshIdToProcess, int whatPlaneProjectTriangleTo, VertCoord tempVars[]) {
	const Point *v1CandidateEdge = (candidateEdge.first>=0)?&vertices[meshIdToProcess][candidateEdge.first]:&vertices[2][-candidateEdge.first -1];
  const Point *v2CandidateEdge = (candidateEdge.second>=0)?&vertices[meshIdToProcess][candidateEdge.second]:&vertices[2][-candidateEdge.second -1];
  
  //cerr << "Testing intersection: " << candidateEdge.first << " " << candidateEdge.second << endl;

  for(const pair<int,int> &edgeToTest:edgesToTest) {


  	const Point *v1EdgeToTest = (edgeToTest.first>=0)?&vertices[meshIdToProcess][edgeToTest.first]:&vertices[2][-edgeToTest.first -1];
  	const Point *v2EdgeToTest = (edgeToTest.second>=0)?&vertices[meshIdToProcess][edgeToTest.second]:&vertices[2][-edgeToTest.second -1];
  	
  	//check if edges (v1EdgeToTest,v2EdgeToTest) intersects (v1CandidateEdge,v2CandidateEdge) (it is ok if they intersect at common vertices..)
  	//if (v1EdgeToTest==v1CandidateEdge || v1EdgeToTest==v2CandidateEdge || v2EdgeToTest==v1CandidateEdge || v2EdgeToTest==v2CandidateEdge )
  	//	continue; //they intersect at common vertices...

  	bool inter = doIntersect(*v1CandidateEdge,*v2CandidateEdge,*v1EdgeToTest,*v2EdgeToTest,whatPlaneProjectTriangleTo,tempVars);
  	
    if(inter) ctEdgeIntersect++;
    else ctEdgeDoNotIntersect++;

  	if(inter) 
  		return true;

  }
  return false;
}



//TODO: avoid memory allocation!!!
bool isInTriangleProj(const Point &p, const Point &p0,const Point &p1,const Point &p2, int whatPlaneProjectTriangleTo) {
  if ( p==p0 || p==p1 || p==p2) return false; //is the point directly above a vertex of the triangle?

  //what coordinates should we check during computations?
  //if the points are projected to z=0 --> we have to use x and y
  //if the points are projected to y=0 --> we have to use x and z
  //...............................x=0 --> we have to use y and z
  int coord1=0,coord2=1;
  if(whatPlaneProjectTriangleTo==PLANE_Y0) {
    coord1=2;
    coord2=0;
  } else if(whatPlaneProjectTriangleTo==PLANE_X0) {
    coord1=1;
    coord2=2;
  }

  
  VertCoord denominator = ((p1[coord2] - p2[coord2])*(p0[coord1] - p2[coord1]) + (p2[coord1] - p1[coord1])*(p0[coord2] - p2[coord2]));
  if (denominator==0) { //TODO: check this.... degenerate triangles or vertical triangles (a segment never intersects a vertical triangle...)
    return false;
  }
  VertCoord a = ((p1[coord2] - p2[coord2])*(p[coord1] - p2[coord1]) + (p2[coord1] - p1[coord1])*(p[coord2] - p2[coord2])) / denominator;
  if ( a<=0 || a >=1) return false;
  
  VertCoord b = ((p2[coord2] - p0[coord2])*(p[coord1] - p2[coord1]) + (p0[coord1] - p2[coord1])*(p[coord2] - p2[coord2])) / denominator;

  if (b<=0 || b>=1) return false;
  VertCoord c = 1 - a - b;
  
  return 0 < c && c < 1; 
}


//returns true iff one of the vertices (represented by their id) is completely inside the triangle a,b,c (represented by ids of vertices)
//the positives ids represent vertices from the meshId mesh
//all vertices are coplanar...thus, we just check if a vertex is in the projection of a,b,c to z=0 (TODO: consider vertical triangle..try to project to x=0 and also y=0 depending on the triangle...)
//we check if the point is completely inside the triangle (i.e., not in the border or outside..)
bool isThereAVertexInsideTriangle(const vector<int> &vertices,const int p0, const int p1, const int p2, const int meshId, int whatPlaneProjectTriangleTo) {
	for(const int v:vertices) {
		if(isInTriangleProj(*getPointFromVertexId(v,meshId),*getPointFromVertexId(p0,meshId),*getPointFromVertexId(p1,meshId),*getPointFromVertexId(p2,meshId), whatPlaneProjectTriangleTo)) {
			return true;
		}
	}
	return false;
}

//seedEdge is an edge containint originalTriangleRetesselated.p[0] and that is collinear with (originalTriangleRetesselated.p[0],originalTriangleRetesselated.p[1])
//it will be use to determine the orientation of the retesselated triangles...
//meshIdToProcess is the id of the mesh of the triangle being split..
void createNewTrianglesFromRetesselationAndOrient(const set<pair<int,int> > &edgesUsedInThisTriangle,const Triangle &originalTriangleRetesselated,vector<TriangleNoBB> &newTrianglesFromRetesselation,  pair<int,int> seedEdge,const int meshIdToProcess, const int  whatPlaneProjectTriangleTo) {
	map<int,int>  verticesStartingFrom0; //TODO: optimize this part...
	vector<int> vertexFrom0ToVertexId;

	vector<pair<int,int> > edgesInTriangleCountingFrom0;
	for(const pair<int,int> &p:edgesUsedInThisTriangle) {
		assert(p.first!=p.second);
		if(verticesStartingFrom0.count(p.first)==0) {
			verticesStartingFrom0[p.first] = vertexFrom0ToVertexId.size();
			vertexFrom0ToVertexId.push_back(p.first);
		}
		if(verticesStartingFrom0.count(p.second)==0) {
			verticesStartingFrom0[p.second] = vertexFrom0ToVertexId.size();
			vertexFrom0ToVertexId.push_back(p.second);
		}

		edgesInTriangleCountingFrom0.push_back(make_pair(verticesStartingFrom0[p.first],verticesStartingFrom0[p.second]));
	}
	int numVertices = vertexFrom0ToVertexId.size();

	int idOfTV0 = verticesStartingFrom0[originalTriangleRetesselated.p[0]]; //what is the label (in our graph) of the two vertices (v[0] and v[1]) from the original triangle 

	bool seedEdgeSameOrientationOriginalTriangle = seedEdge.first==originalTriangleRetesselated.p[0];

	seedEdge.first = verticesStartingFrom0[seedEdge.first];
	seedEdge.second = verticesStartingFrom0[seedEdge.second];

	bool adjMatrix[numVertices][numVertices];
	for(int i=0;i<numVertices;i++)
		for(int j=0;j<numVertices;j++)
			adjMatrix[i][j] = false;

	//cerr << "Number of edges to triangulate: " << edgesInTriangleCountingFrom0.size() << endl;
	//cerr << "Number of vertices: " << numVertices << endl;
	for(const pair<int,int> &p:edgesInTriangleCountingFrom0) {
		//cerr << p.first << "  " << p.second << endl;
		adjMatrix[p.first][p.second] = adjMatrix[p.second][p.first] = true;
		//cerr << p.first << " " << p.second << endl;
		assert(p.first!=p.second);
	}

	set< pair<int,pair<int,int>> > trianglesProcessed;
	queue<pair<pair<int,int>,pair<int,int> > > edgesToProcess;

	if(seedEdgeSameOrientationOriginalTriangle)
		edgesToProcess.push( make_pair(seedEdge,make_pair(originalTriangleRetesselated.above, originalTriangleRetesselated.below)));
	else 
		edgesToProcess.push( make_pair(seedEdge,make_pair(originalTriangleRetesselated.below,originalTriangleRetesselated.above)));

	set<pair<int,int> > processedEdges;
	while(!edgesToProcess.empty()) {
		int i = edgesToProcess.front().first.first;
		int j = edgesToProcess.front().first.second;
		int above = edgesToProcess.front().second.first;
		int below = edgesToProcess.front().second.second;
		edgesToProcess.pop();


		pair<int,int> edge(min(i,j),max(i,j)); //let's use this to uniquely identify each edge so that we do not process the same edge twice...
		if(processedEdges.count(edge)!=0) {
			continue; //the triangles containting that edge were already processed...
		}
		processedEdges.insert(edge);

		
		//we have an edge (i,j) with the same orientation as originalTriangleRetesselated
		for(int k=0;k<numVertices;k++) {
			if(adjMatrix[i][k] && adjMatrix[j][k]) {

				pair<int,pair<int,int>> triangle;
				if (i<j && i<k) {
					triangle.first = i;

					if(j<k) {
						triangle.second.first = j;
						triangle.second.second = k;
					} else {
						triangle.second.first = k;
						triangle.second.second = j;
					}
				} else 	if (j<i && j<k) {
					triangle.first = j;

					if(i<k) {
						triangle.second.first = i;
						triangle.second.second = k;
					} else {
						triangle.second.first = k;
						triangle.second.second = i;
					}
				} else {
					triangle.first = k;

					if(i<j) {
						triangle.second.first = i;
						triangle.second.second = j;
					} else {
						triangle.second.first = j;
						triangle.second.second = i;
					}
				}
			
				if(trianglesProcessed.count(triangle)!=0) {
					continue;
				}

				//if a vertex from the retesselation is inside this candidate triangle, we will not create this triangle..
				//test all "numVertices" vertices in vertexFrom0ToVertexId
				if(isThereAVertexInsideTriangle(vertexFrom0ToVertexId,vertexFrom0ToVertexId[i],vertexFrom0ToVertexId[j],vertexFrom0ToVertexId[k],meshIdToProcess,whatPlaneProjectTriangleTo)) {
					//cerr << "Ignoring triangle: " << vertexFrom0ToVertexId[i] << " " << vertexFrom0ToVertexId[j] << " " << vertexFrom0ToVertexId[k] << endl;
					continue;
				}

				//cerr << "Adding: " << vertexFrom0ToVertexId[i] << " " << vertexFrom0ToVertexId[j] << " " << vertexFrom0ToVertexId[k] << endl;
					//we have a triangle i,j,k...
				assert(vertexFrom0ToVertexId[i] != vertexFrom0ToVertexId[j]);
				assert(vertexFrom0ToVertexId[i] != vertexFrom0ToVertexId[k]);
				assert(vertexFrom0ToVertexId[j] != vertexFrom0ToVertexId[k]);
				newTrianglesFromRetesselation.push_back(TriangleNoBB(vertexFrom0ToVertexId[i],vertexFrom0ToVertexId[j],vertexFrom0ToVertexId[k],above,below)); //TODO: fix orientation...
				trianglesProcessed.insert(triangle);
				edgesToProcess.push(make_pair(make_pair(j,k),make_pair(below,above)));
				edgesToProcess.push(make_pair(make_pair(k,i),make_pair(below,above)));
			}
		}
	}

	/*
	for(const pair<int,int> &p:edgesInTriangleCountingFrom0) {
			int i=p.first;
			int j=p.second;
			for(int k=0;k<numVertices;k++) {
				if(adjMatrix[i][k] && adjMatrix[j][k]) {
					//we have a triangle i,j,k...
					assert(vertexFrom0ToVertexId[i] != vertexFrom0ToVertexId[j]);
					assert(vertexFrom0ToVertexId[i] != vertexFrom0ToVertexId[k]);
					assert(vertexFrom0ToVertexId[j] != vertexFrom0ToVertexId[k]);
					newTrianglesFromRetesselation.push_back(TriangleNoBB(vertexFrom0ToVertexId[i],vertexFrom0ToVertexId[j],vertexFrom0ToVertexId[k],-1,0)); //TODO: fix orientation...
				}
			}
		}
	*/
}


int numVerticesInEdges = 0;
//stores some counters (statistics)
struct StatisticsAboutRetesseation {
  StatisticsAboutRetesseation() {
    ctTrianglesRetesselate = 0;
    ctEdgesActuallyInsertedInRetesselation =0;
    ctEdgesInsertedBecauseOfConstraint=0;
    ctVerticesInsertedTrianglesToRetesselate=0;
    ctEdgesTestedToInsertInRetesselation=0;
    simple = difficult = 0;
    ctConvexPolygonsInTriangleRetesselations  = ctConcavePolygonsInTriangleRetesselations = 0;
    numConnectedPolygonSubdivisionOfTriangles = numDisconnectedPolygonSubdivisionOfTriangles = 0;
  }

  void printStatsSoFar() {
    cerr << numVerticesInEdges << endl;
    cerr << "Num triangles to retesselate                                 : "  << ctTrianglesRetesselate << "\n";
    cerr << "Total number of edges we tried to insert during retesselation: "  << ctEdgesTestedToInsertInRetesselation << "\n";  
    cerr << "Total edges tried insert and actually inserted               : " << ctEdgesActuallyInsertedInRetesselation  <<"\n";
    cerr << "Total edges inserted because are constraint edges            : " << ctEdgesInsertedBecauseOfConstraint <<"\n";
    cerr << "Total number of edges in retesselated triangles (add prev. 2): " << ctEdgesInsertedBecauseOfConstraint+ctEdgesActuallyInsertedInRetesselation <<"\n";
    cerr << "Total number of vertices in retesselated triangles           : "  << ctVerticesInsertedTrianglesToRetesselate << "\n";
    cerr << "Simple case triangles                                        : " << simple <<"\n";
    cerr << "Difficult case triangles                                     : " << difficult << "\n";
    cerr << "Number of convex polygons in triangle retesselations         : " << ctConvexPolygonsInTriangleRetesselations << "\n";
    cerr << "Number of concave polygons in triangle retesselations        : " << ctConcavePolygonsInTriangleRetesselations << "\n";
    cerr << "Number of connected polygon subdivision in triangles         : " << numConnectedPolygonSubdivisionOfTriangles << "\n";
    cerr << "Number of disconnected polygon subdivision in triangles      : " << numDisconnectedPolygonSubdivisionOfTriangles << "\n";
  }
  long long simple, difficult;
  long long ctTrianglesRetesselate;
  long long ctEdgesActuallyInsertedInRetesselation;
  long long ctEdgesInsertedBecauseOfConstraint;
  long long ctVerticesInsertedTrianglesToRetesselate;
  long long ctEdgesTestedToInsertInRetesselation;
  long long ctConvexPolygonsInTriangleRetesselations;
  long long ctConcavePolygonsInTriangleRetesselations;
  long long numDisconnectedPolygonSubdivisionOfTriangles;
  long long numConnectedPolygonSubdivisionOfTriangles;
};


//temporary variables... we avoid reallocating them every time this function is calles
struct TempVarsRetesselateTriangleFunction {
  vector<int> polygons;
  vector<array<int,3> > wedgesTemp;
  vector<pair<int,int> > raggedArraySortedEdges;
  vector<bool> usedWedgesTemp;
  VertCoord tempVertCoords[8];
};

//given a triangle t, this function will split t at the intersection edges (intersection with other triangles)
//and retesselate t, creating more triangles.
void retesselateTriangle(const vector<pair<int,int> > &edgesUsingVertexId, const vector<int> &edgesFromIntersection,const Triangle &t, const int meshWhereTriangleIs, vector<TriangleNoBB> &newTrianglesGeneratedFromRetesselation,VertCoord tempVars[], StatisticsAboutRetesseation &statistics) {
  //A triangle (except if it is degenerate to a single point) cannot be perpendicular to all 3 planes (x=0,y=0,z=0)
  //Thus, we chose one of these planes to project the triangle during computations...
  const int whatPlaneProjectTriangleTo = getPlaneTriangleIsNotPerpendicular(vertices[meshWhereTriangleIs][t.p[0]], vertices[meshWhereTriangleIs][t.p[1]],vertices[meshWhereTriangleIs][t.p[2]], tempVars);

  //we need to "retesselate" t and orient all the new triangles properly
  //this set will store the edges we used in the retesselated triangle
  //we need to choose what edges to create and, then, use these new edges to reconstruct the triangulation..
  //This set will have the edges that will form the retriangulation of ts..
  set<pair<int,int> > edgesUsedInThisTriangle; //TODO: maybe use unordered_set (see overhead difference...)


  int ct =0;
  //The edges from intersection will, necessarelly, be in the triangulation...
  for(int edgeId:edgesFromIntersection) {
    //cerr << edgesUsingVertexId[edgeId].first << " " << edgesUsingVertexId[edgeId].second << endl;
    assert(edgesUsingVertexId[edgeId].first != edgesUsingVertexId[edgeId].second);
    edgesUsedInThisTriangle.insert(edgesUsingVertexId[edgeId]);   

    #ifdef COLLECT_STATISTICS
      #pragma omp atomic
      statistics.ctEdgesInsertedBecauseOfConstraint++;   
    #endif
  }


  //Now, we need to add more edges to fill the parts of ts that still do not form triangle
  vector<int> verticesToTesselate;
  //what vertices will we have in the new triangulation of ts? (the original vertices + the vertices of the edges from intersections)
  for(const auto &p:edgesUsedInThisTriangle) {
    verticesToTesselate.push_back(p.first);
    verticesToTesselate.push_back(p.second);
  }
  verticesToTesselate.push_back(t.p[0]);
  verticesToTesselate.push_back(t.p[1]);
  verticesToTesselate.push_back(t.p[2]);
  sort(verticesToTesselate.begin(),verticesToTesselate.end());
  auto newEnd = unique(verticesToTesselate.begin(),verticesToTesselate.end());
  verticesToTesselate.resize(newEnd-verticesToTesselate.begin());
  int numVTriangle = verticesToTesselate.size();

  #ifdef COLLECT_STATISTICS
    #pragma omp atomic
    statistics.ctVerticesInsertedTrianglesToRetesselate += numVTriangle;
  #endif

  

  const Point &triangleVertex0 = vertices[meshWhereTriangleIs][t.p[0]];
  const Point &triangleVertex1 = vertices[meshWhereTriangleIs][t.p[1]];
  const Point &triangleVertex2 = vertices[meshWhereTriangleIs][t.p[2]];
  //cerr << triangleVertex0[0].get_d( ) << " " << triangleVertex0[1].get_d() << " " << triangleVertex0[2].get_d() << endl;;
  //cerr << triangleVertex1[0].get_d( ) << " " << triangleVertex1[1].get_d() << " " << triangleVertex1[2].get_d() << endl;;
  //cerr << triangleVertex2[0].get_d( ) << " " << triangleVertex2[1].get_d() << " " << triangleVertex2[2].get_d() << endl;;
  //cerr << "Project to: " << whatPlaneProjectTriangleTo << endl;


  //verticesIncidentEachEdgeOriginalTriangle[0] is for vertex t.p[0]-t.p[1]
  //verticesIncidentEachEdgeOriginalTriangle[1] is for vertex t.p[1]-t.p[2]
  //verticesIncidentEachEdgeOriginalTriangle[2] is for vertex t.p[2]-t.p[0]
  list<int> verticesIncidentEachEdgeOriginalTriangle[3];
  for(int i=0;i<numVTriangle;i++) {
    int v = verticesToTesselate[i];
    if(v==t.p[0] || v==t.p[1] || v==t.p[2]) continue;
    //for each vertex v that is not one of the original vertices of the triangle, let's see if v is in one of the edges of t...

    const Point &vertex = *getPointFromVertexId(v, meshWhereTriangleIs);
    //cerr << "V: " << vertex[0].get_d( ) << " " << vertex[1].get_d() << " " << vertex[2].get_d() << endl;;
    int o1 = orientation(triangleVertex0, triangleVertex1,  vertex,whatPlaneProjectTriangleTo,tempVars);
    
    
    //cerr << o1 << o2 << o3 << endl;
    if(o1==0) {
      verticesIncidentEachEdgeOriginalTriangle[0].push_back(i);
      continue;
    }
    int o2 = orientation(triangleVertex1, triangleVertex2,  vertex,whatPlaneProjectTriangleTo,tempVars);
    if(o2==0) {
      //cerr << "12" << endl;
      verticesIncidentEachEdgeOriginalTriangle[1].push_back(i);
      continue;
    }
    int o3 = orientation(triangleVertex2, triangleVertex0,  vertex,whatPlaneProjectTriangleTo,tempVars);
    if(o3==0) {
     // cerr << 20 << endl;
      verticesIncidentEachEdgeOriginalTriangle[2].push_back(i);
      continue;
    }
  }
  //cerr << "In edges: " << numVerticesInEdges << "\n";
 // cerr << "Vertices: " << numVTriangle << "\n";
  //cerr << "Constraint edges: " << edgesFromIntersection.size() << "\n\n";

  pair<int,int> seedEdge(-1,-1); //we need a seed output edge that will be oriented in the same way t is oriented..
                           //we will use this seed to reorient all the triangles resulting from the retesselation of t

  bool allSimple = true;
  //no vertex intersect the edge t.p[0]-t.p[1] --> it will be in the triangulation!!!
  int verticesIncidentEdge;
  verticesIncidentEdge = verticesIncidentEachEdgeOriginalTriangle[0].size();
  if(verticesIncidentEdge==0) {
    int v1 = t.p[0];
    int v2 = t.p[1];
            
    if(v2<v1) swap(v1,v2); //the first vertex is always smaller than the second in the edges... (to simplify tests)
    pair<int,int> candidateEdge(v1,v2);

    if (v1 == t.p[0] && v2 == t.p[1]) {
      seedEdge.first = t.p[0];
      seedEdge.second = t.p[1];
    } else if  (v1 == t.p[1] && v2 == t.p[0]) {
      seedEdge.first = t.p[1];
      seedEdge.second = t.p[0];
    }

    #ifdef COLLECT_STATISTICS
      #pragma omp atomic
      statistics.ctEdgesActuallyInsertedInRetesselation++;
    #endif

    assert(candidateEdge.first != candidateEdge.second);
    edgesUsedInThisTriangle.insert(candidateEdge);
  } else {
    if(verticesIncidentEdge==1) {
      int v1 = t.p[0];
      int v2 = t.p[1]; 
      int v =  verticesToTesselate[*(verticesIncidentEachEdgeOriginalTriangle[0].begin())];          
      
      pair<int,int> candidateEdge;    

      if(v1<v) {candidateEdge.first = v1; candidateEdge.second = v;}
      else {candidateEdge.first = v; candidateEdge.second = v1;}

      #ifdef COLLECT_STATISTICS
        #pragma omp atomic
        statistics.ctEdgesActuallyInsertedInRetesselation++;
      #endif

      assert(candidateEdge.first != candidateEdge.second);
      edgesUsedInThisTriangle.insert(candidateEdge);

      if(v2<v) {candidateEdge.first = v2; candidateEdge.second = v;}
      else {candidateEdge.first = v; candidateEdge.second = v2;}

      #ifdef COLLECT_STATISTICS
        #pragma omp atomic
        statistics.ctEdgesActuallyInsertedInRetesselation++;
      #endif
      assert(candidateEdge.first != candidateEdge.second);
      edgesUsedInThisTriangle.insert(candidateEdge);
    } else {
      allSimple= false;
    }
  }

  verticesIncidentEdge = verticesIncidentEachEdgeOriginalTriangle[1].size();
  if(verticesIncidentEdge==0) {
    int v1 = t.p[1];
    int v2 = t.p[2];
            
    if(v2<v1) swap(v1,v2); //the first vertex is always smaller than the second in the edges... (to simplify tests)
    pair<int,int> candidateEdge(v1,v2);    

    #ifdef COLLECT_STATISTICS
        #pragma omp atomic
        statistics.ctEdgesActuallyInsertedInRetesselation++;
    #endif

    assert(candidateEdge.first != candidateEdge.second);
    edgesUsedInThisTriangle.insert(candidateEdge);
  } else {
    if(verticesIncidentEdge==1) {
      int v1 = t.p[1];
      int v2 = t.p[2]; 
      int v =  verticesToTesselate[*(verticesIncidentEachEdgeOriginalTriangle[1].begin())];          
      
      pair<int,int> candidateEdge;    

      if(v1<v) {candidateEdge.first = v1; candidateEdge.second = v;}
      else {candidateEdge.first = v; candidateEdge.second = v1;}

      #ifdef COLLECT_STATISTICS
        #pragma omp atomic
        statistics.ctEdgesActuallyInsertedInRetesselation++;
      #endif
      assert(candidateEdge.first != candidateEdge.second);
      edgesUsedInThisTriangle.insert(candidateEdge);

      if(v2<v) {candidateEdge.first = v2; candidateEdge.second = v;}
      else {candidateEdge.first = v; candidateEdge.second = v2;}

      #ifdef COLLECT_STATISTICS
        #pragma omp atomic
        statistics.ctEdgesActuallyInsertedInRetesselation++;
      #endif
      assert(candidateEdge.first != candidateEdge.second);
      edgesUsedInThisTriangle.insert(candidateEdge);
    } else {
      allSimple= false;
    }
  }

  verticesIncidentEdge = verticesIncidentEachEdgeOriginalTriangle[2].size();
  if(verticesIncidentEdge==0) {
    int v1 = t.p[2];
    int v2 = t.p[0];
            
    if(v2<v1) swap(v1,v2); //the first vertex is always smaller than the second in the edges... (to simplify tests)
    pair<int,int> candidateEdge(v1,v2);

    #ifdef COLLECT_STATISTICS
      #pragma omp atomic
      statistics.ctEdgesActuallyInsertedInRetesselation++;
    #endif
    assert(candidateEdge.first != candidateEdge.second);
    edgesUsedInThisTriangle.insert(candidateEdge);
  } else {
    if(verticesIncidentEdge==1) {
      int v1 = t.p[2];
      int v2 = t.p[0]; 
      int v =  verticesToTesselate[*(verticesIncidentEachEdgeOriginalTriangle[2].begin())];          
      
      pair<int,int> candidateEdge;    

      if(v1<v) {candidateEdge.first = v1; candidateEdge.second = v;}
      else {candidateEdge.first = v; candidateEdge.second = v1;}
      #ifdef COLLECT_STATISTICS
        #pragma omp atomic
        statistics.ctEdgesActuallyInsertedInRetesselation++;
      #endif
      assert(candidateEdge.first != candidateEdge.second);
      edgesUsedInThisTriangle.insert(candidateEdge);

      if(v2<v) {candidateEdge.first = v2; candidateEdge.second = v;}
      else {candidateEdge.first = v; candidateEdge.second = v2;}

      #ifdef COLLECT_STATISTICS
        #pragma omp atomic
        statistics.ctEdgesActuallyInsertedInRetesselation++;
      #endif
      assert(candidateEdge.first != candidateEdge.second);
      edgesUsedInThisTriangle.insert(candidateEdge);
    } else {
      allSimple= false;
    }
  }
  
  #ifdef COLLECT_STATISTICS        
      if(allSimple) {
        #pragma omp atomic
        statistics.simple++;
      }
      else {
        #pragma omp atomic
        statistics.difficult++;
      }
  #endif

  if(allSimple) {
    //vertices should connect: internal - boundary , two boundary iff from intersection...

    set<pair<int,int> > internalEdgesThisTriangle;
    for(int edgeId:edgesFromIntersection) {   
      internalEdgesThisTriangle.insert(edgesUsingVertexId[edgeId]);       
    }
    /*set<int> verticesOriginalTriangle;
    for(int i=0;i<3;i++) {
      verticesOriginalTriangle.insert(t.p[i]);
    }*/
    

    for(int i=0;i<numVTriangle;i++) {
      int v1 = verticesToTesselate[i];
      bool v1Original = (v1==t.p[0])||(v1==t.p[1])||(v1==t.p[2]); 
     // bool v1Original = verticesOriginalTriangle.count(v1)!=0; 
      for(int j=i+1;j<numVTriangle;j++) {
        
        int v2 = verticesToTesselate[j];
         
        bool v2Original = (v2==t.p[0])||(v2==t.p[1])||(v2==t.p[2]); 
        if(v1Original && v2Original) continue; //v1 is a boundary vertex...

        pair<int,int> candidateEdge;
        if(v1<v2) {candidateEdge.first = v1; candidateEdge.second = v2;}
        else {candidateEdge.first = v2; candidateEdge.second = v1;}

        if(edgesUsedInThisTriangle.count(candidateEdge)!=0) continue; 
        assert(candidateEdge.first != candidateEdge.second);

        #ifdef COLLECT_STATISTICS   
          #pragma omp atomic
          statistics.ctEdgesTestedToInsertInRetesselation++;
        #endif
        
        //we do not need to test the candidate edges for intersection with the edges in the boundary of the triangle
        if(!intersects(candidateEdge,internalEdgesThisTriangle,meshWhereTriangleIs,whatPlaneProjectTriangleTo,tempVars)) {        

          #ifdef COLLECT_STATISTICS   
            #pragma omp atomic
            statistics.ctEdgesActuallyInsertedInRetesselation++;
          #endif
          edgesUsedInThisTriangle.insert(candidateEdge); //if it does not intersect, we can safely add this edge to the triangulation...
          internalEdgesThisTriangle.insert(candidateEdge); 
          //cerr << "Dont intersect\n";
        }

      }
    }

  } else { 
    //for each pair of vertices (forming an edge), let's try to add this edge e to the triangulation
    //

    /*set<pair<int,int> > internalEdgesThisTriangle;

    set<int> verticesOriginalTriangle;
    for(int i=0;i<3;i++) {
      verticesOriginalTriangle.insert(t.p[i]);
    } */ 

    for(int i=0;i<numVTriangle;i++) {
      int v1 = verticesToTesselate[i];

     // bool v1Original = verticesOriginalTriangle.count(v1)!=0; 
      for(int j=i+1;j<numVTriangle;j++) {
        
        int v2 = verticesToTesselate[j];


              
        assert(v1<v2); //the vector was previously sorted! also, all elements should be unique!
        //let's try to insert edge (v1,v2)...

        //bool v2Original = verticesOriginalTriangle.count(v2)!=0;


        pair<int,int> candidateEdge(v1,v2);
        if(edgesUsedInThisTriangle.count(candidateEdge)!=0) continue; //the edge was already used...


       /* if(allSimple && v1Original && v2Original) {
          //cerr << "continuing..." << endl;
          continue;
        }*/
        #ifdef COLLECT_STATISTICS   
          #pragma omp atomic
          statistics.ctEdgesTestedToInsertInRetesselation++;
        #endif

        //if edge e=(v1,v2) does not intersect other edges already inserted in this triangle,
        //we will add e to the new triangulation...
        //any triangulation is fine as soon as the edges from the intersection of the meshes are there...
        if(!intersects(candidateEdge,edgesUsedInThisTriangle,meshWhereTriangleIs,whatPlaneProjectTriangleTo,tempVars)) {        

          #ifdef COLLECT_STATISTICS   
            #pragma omp atomic
            statistics.ctEdgesActuallyInsertedInRetesselation++;
          #endif
          edgesUsedInThisTriangle.insert(candidateEdge); //if it does not intersect, we can safely add this edge to the triangulation...
        
          //cerr << "Dont intersect\n";
        } else {
          //cerr << "Intersects\n";
        }
      }
    }
  }


  if(seedEdge.first == -1) { //the edge t.p[0] - t.p[1] is not in the output because other edges intersects it...
    for(const pair<int,int> &e:edgesUsedInThisTriangle) 
      if(e.first == t.p[0] ) { //let's try to find an edge (t.p[0],x) (or (x,t.p[0])) collinear with (t.p[0], t.p[1])
        //is e.second on line (t.p[0]-t.p[1]) ??
        //TODO: this is computed with the projected triangle... we should consider vertical triangles...
        if (orientation(vertices[meshWhereTriangleIs][t.p[0]], vertices[meshWhereTriangleIs][t.p[1]], *getPointFromVertexId(e.second, meshWhereTriangleIs),whatPlaneProjectTriangleTo, tempVars)==0) {
          seedEdge.first = t.p[0];
          seedEdge.second = e.second;
          //cerr << "Aqui" << endl;
          break;
        }
      } else if (e.second == t.p[0]) {
        if (orientation(vertices[meshWhereTriangleIs][t.p[0]], vertices[meshWhereTriangleIs][t.p[1]], *getPointFromVertexId(e.first, meshWhereTriangleIs),whatPlaneProjectTriangleTo,tempVars)==0) {
          seedEdge.first = e.first;
          seedEdge.second = t.p[0];
          //cerr << "Ali" << endl;
          break;
        }
      }         
  }



  createNewTrianglesFromRetesselationAndOrient(edgesUsedInThisTriangle,t,newTrianglesGeneratedFromRetesselation,seedEdge,meshWhereTriangleIs, whatPlaneProjectTriangleTo);
          

  //cerr << meshWhereTriangleIs << " " << i << " " << myNewTrianglesFromRetesselation.size() << endl;
 // for(auto &elem:edgesUsedInThisTriangle)
   // myNewTriEdgesFromEachMap.push_back(elem); ;//newTriEdgesFromEachMap[meshWhereTriangleIs].push_back(elem);   
}



//returns the index i of the vector v, where v[i][0] = key.first and v[i][1] = key.second 
int binarySearch(vector<array<int,3> > &v,const pair<int,int> &key) {
  int lo = 0, hi = v.size()-1, mid;
  while (lo<=hi) {
    mid = lo + (hi-lo)/2;
    if(v[mid][0]==key.first && v[mid][1]==key.second) 
      return mid;
    else if ( (v[mid][0] < key.first) || ((v[mid][0] == key.first) && (v[mid][1] < key.second)) ) //is the element v[mid] smaller than the key?
      lo = mid+1;
    else 
      hi = mid-1;
  }
  return -1;
}


//tempVars should have at least size 4
void sortTriangleEdgesByAngle(int meshIdToProcess, int whatPlaneProjectTo, vector<pair<int,int> > &edgesToSort,VertCoord tempVars[]) {
  int xCoord = 0;
  int yCoord = 1;
  if(whatPlaneProjectTo==PLANE_Y0) {
    xCoord= 2;
    yCoord= 0;
  } else if(whatPlaneProjectTo==PLANE_X0) {
    xCoord= 1;
    yCoord= 2;    
  } 

  auto edgesComparisonFunctionProject = [&](const pair<int,int> &e1, const pair<int,int> &e2) { 
      if(e1.first!=e2.first) return e1.first<e2.first; //sort first by the first vertex..
      //if the first vertex is equal, we need to sort basing on the angle...
      //let's check if the cross product between e1 and e2 is positive... cross = v1.x * v2.y - v2.x * v1.y;
      const Point &e10 = *getPointFromVertexId(e1.first, meshIdToProcess); //if the vertex is shared (that is, it does not belong only to meshIdToProcess) its id will be negative (the function 
      const Point &e11 = *getPointFromVertexId(e1.second, meshIdToProcess); // getPointFromVertexId treats this automatically).

      const Point &e20 = *getPointFromVertexId(e2.first, meshIdToProcess); //if the vertex is shared (that is, it does not belong only to meshIdToProcess) its id will be negative (the function 
      const Point &e21 = *getPointFromVertexId(e2.second, meshIdToProcess); // getPointFromVertexId treats this automatically).

      /*
      const VertCoord v1x = (e11[xCoord]-e10[xCoord]);
      const VertCoord v1y = (e11[yCoord]-e10[yCoord]);
      const VertCoord v2x = (e21[xCoord]-e20[xCoord]);
      const VertCoord v2y = (e21[yCoord]-e20[yCoord]);


      if(sgn(v1y)>=0 && sgn(v2y)<0) return true; //check if the two vectors are in different sides of the x axis...
      if(sgn(v1y)<0 && sgn(v2y)>=0) return false;
      if(sgn(v1y) ==0 && sgn(v2y)==0) { //both are on the x axis... 
        if(sgn(v1x)>=0) return true; //they both cannot have the same sign simultaneously (otherwise they would coincide..)
        else return false;
      }


      //is the cross product positive?
      const VertCoord component1 = v1x*v2y; //v1.x * v2.y 
      const VertCoord component2 = v2x*v1y; //v2.x*v1.y

      //the cross product is = component1-component2
      //is it positive?
      return (component1 > component2);*/
      VertCoord &v1x = tempVars[0];
      VertCoord &v1y = tempVars[1];
      VertCoord &v2x = tempVars[2];
      VertCoord &v2y = tempVars[3];
      v1x = e11[xCoord];
      v1x -= e10[xCoord];

      v1y = e11[yCoord];
      v1y -= e10[yCoord];

      v2x = e21[xCoord];
      v2x -= e20[xCoord];

      v2y = e21[yCoord];
      v2y -= e20[yCoord];

      const int sgnV1y = sgn(v1y);
      const int sgnV2y = sgn(v2y);

      if(sgnV1y>=0 && sgnV2y<0) return true; //check if the two vectors are in different sides of the x axis...
      if(sgnV1y<0 && sgnV2y>=0) return false;
      if(sgnV1y==0 && sgnV2y==0) { //both are on the x axis... 
        if(sgn(v1x)>=0) return true; //they both cannot have the same sign simultaneously (otherwise they would coincide..)
        else return false;
      }

      //is the cross product positive?
      v1x *= v2y; //const VertCoord component1 = v1x*v2y; //v1.x * v2.y 
      v2x *= v1y; //const VertCoord component2 = v2x*v1y; //v2.x*v1.y

      //the cross product is = component1-component2
      //is it positive?
      return v1x > v2x;//return (component1 > component2);
  };

  sort(edgesToSort.begin(),edgesToSort.end(), edgesComparisonFunctionProject );
}








#define ERROR_CODE -2123123123
int binarySearchFirstOccurrency(const vector<pair<int,int> > &v,int key) {
  int mid;
  int lo = 0;
  int hi = v.size()-1;
  while (lo < hi) {
      mid = lo + (hi-lo)/2;   // note: division truncates
      if (v[mid].first>=key)
         hi = mid;
      else
         lo = mid+1;
  }
          
  if (v[lo].first!=key)
      return ERROR_CODE;                // p(x) is true for all x in S!
      
  return lo;         // lo is the least x for which p(x) is true
}


bool checkIfPolygonsAreConnected(const vector<pair<int,int> > &raggedArraySortedEdges) {
  queue<int> verticesToProcess;
  set<int> processedVertices;
  verticesToProcess.push(raggedArraySortedEdges[0].first);
  processedVertices.insert(raggedArraySortedEdges[0].first);

  
  

  //visit all vertices reachable from v...
  while(!verticesToProcess.empty()) {
    int v = verticesToProcess.front();
    verticesToProcess.pop();

    int posStartVRaggedArray =  binarySearchFirstOccurrency(raggedArraySortedEdges,v);
    //cerr << "v and pos v " << v << " " << posStartVRaggedArray << endl;
    for(int i=posStartVRaggedArray;i<raggedArraySortedEdges.size()&&raggedArraySortedEdges[i].first==v;i++) {
      int adj = raggedArraySortedEdges[i].second; //adj is one of the vertices adjacent to v
      if(processedVertices.count(adj)==0) {
        processedVertices.insert(adj);
        verticesToProcess.push(adj);
      }
    }
  }

  //cerr << "Checking connectivity... " << endl;
  for(const pair<int,int> &p:raggedArraySortedEdges) {
    if(processedVertices.count(p.first)==0) return false; //at least one vertex was not visited...
    //cerr << p.first << " " << p.second << endl;

  } 
 /* cerr << "Visited: " << endl;
  for(int i:processedVertices) cerr << i << endl;  
  */

  return true;
}

//given a triangle t, this function will split t at the intersection edges (intersection with other triangles)
//and retesselate t, creating more triangles.
//each edge is represented by a pair of ids of the vertices connected by that edge
//t is retesselated by:
//- using an algorithm that sorts the edges creating wedges to extract the polygons (whose union is t) that need to be retesselated
//- it is faster to retesselate these polygons individually than to retesselate the whole t (this is true mainly when there are several edges from intersection in t)

//Output: polygons
//-- Polygons will be a list of vertices ids representing the polygons
//-- Example: v1,v2,v3,v1,v9,v12,v10,v9 represents two polygons: v1,v2,v3,v1 and v9,v12,v10.v9

//raggedArraySortedEdges,usedWedgesTemp and wedgesTemp are temporary vectors (we pass them to the function to avoid too much memory reallocation)
//The first numEdgesFromIntersection edges in edgesUsingVertexId are from the intersections
//The last ones are edges from the boundary of the original rectangle being retesselated

//raggedArraySortedEdges will be filled with the edges in the triangle (if we have an edge (a,b) , the array will have both (a,b) and (b,a)). Also,
//at the end of the function the edges in raggedArraySortedEdges will be sorted by the first vertex and, if there is a tie, by the slope of the edge (w.r.t. the projection to whatPlaneProjectTo)
//at the end of this function raggedArraySortedEdges will contain the edges (two directed edges for each edge) sorted
//by the first vertex and by the angle of the edge 
void sortEdgesAndExtractPolygonsFromEdgeListUsingWedges(const vector<pair<int,int> > &edgesUsingVertexId,const int numEdgesFromIntersection, const int meshIdToProcess, const int whatPlaneProjectTo,  TempVarsRetesselateTriangleFunction &tempVars ) {
  //Algorithm: http://ac.els-cdn.com/016786559390104L/1-s2.0-016786559390104L-main.pdf?_tid=7586b72e-a059-11e6-afdd-00000aacb35d&acdnat=1478021856_7213cc56dd2148587a06891102f5d518
   
  
  vector<array<int,3> > &wedgesTemp = tempVars.wedgesTemp;
  vector<pair<int,int> > &raggedArraySortedEdges = tempVars.raggedArraySortedEdges;
  vector<bool> &usedWedgesTemp = tempVars.usedWedgesTemp;
  vector<int> &polygons = tempVars.polygons;

  //First, we need to find all the wedges:

  //Duplicate each undirected edge such that we will have two directed edges (u,v) and (v,u)
  //sort the edges basing first on the first vertex and second on the angle the edge makes with the horizontal line (supposing the triangle is projected to planeProjectTriangleTo)
  //Process the groups in the sorted list:
  //-- A Group is a set of edges with the same first vertex
  //-- Within each group, combine each pair of consecutive edges (including the last with the first) to form a wedge: (a,b) (a,c) --> (c,a,b)
  const int numDirectedEdges = (edgesUsingVertexId.size())*2 ;
  raggedArraySortedEdges.resize(numDirectedEdges);
  const int numOrigEdges = edgesUsingVertexId.size();

  //cerr << "Edge list: " << endl;
  //for(auto p:edgesUsingVertexId) cerr << p.first << " " << p.second << endl; cerr << endl;

  for(int i=0;i<numOrigEdges;i++) { //duplicate the edges
    raggedArraySortedEdges[2*i] = edgesUsingVertexId[i]; //add (u,v)
    raggedArraySortedEdges[2*i+1].first = edgesUsingVertexId[i].second; //add (v,u)
    raggedArraySortedEdges[2*i+1].second = edgesUsingVertexId[i].first;
  }
  /*for(int i=numEdgesFromIntersection;i<numDirectedEdges;i++) { //duplicate the edges
    edgesTemp[2*i] = edgesUsingVertexId[i]; //add (u,v)
  }*/

  int x = 0;
  //now we need to sort the edges basing on the first vertex and, then, on the angle with horizon (considering whatPlaneProjectTo)
  
  sortTriangleEdgesByAngle(meshIdToProcess, whatPlaneProjectTo, raggedArraySortedEdges,tempVars.tempVertCoords);  
  

  wedgesTemp.resize(0);
  int numWedgesFound = 0;
  //Now let's process the groups of edges (set of edges with the same first vertex..) to extract the wedges  
  for(int firstElementGroup =0;firstElementGroup<numDirectedEdges;) {
    for(int currElemnt = firstElementGroup;currElemnt<numDirectedEdges;currElemnt++) {
      //combine the currentElement with the next one to form a wedge...
      int nextElement = currElemnt+1;
      //(a,b)         (a,c) --> (c,a,b)
      //currElement    nextElement

      wedgesTemp.resize(numWedgesFound+1);
      numWedgesFound++;

      array<int,3> &wedgeToAdd = wedgesTemp.back();
      //did we reach the end of a group?
      if(nextElement >= numDirectedEdges || raggedArraySortedEdges[nextElement].first != raggedArraySortedEdges[currElemnt].first) {
        //if yes, we will create a wedge from currElement to the first one! (to complete the wedges from this group..)
        //the next element is actually the first element in this group...
        wedgeToAdd[0] = raggedArraySortedEdges[firstElementGroup].second; //c
        wedgeToAdd[1] = raggedArraySortedEdges[firstElementGroup].first; //a
        wedgeToAdd[2] = raggedArraySortedEdges[currElemnt].second;  //b

        firstElementGroup = currElemnt+1; //let's process the next group!
        break;
      } else {
        //else, we will create a wedge from the current to the next...
        wedgeToAdd[0] = raggedArraySortedEdges[nextElement].second; //c
        wedgeToAdd[1] = raggedArraySortedEdges[nextElement].first; //a
        wedgeToAdd[2] = raggedArraySortedEdges[currElemnt].second;  //b
      }
    }
  }

  //Now, wedgesTemp contain all the wedges we need to extract the polygons!

  //Once we have the wedges, we need to group the wedges to form the polygons:
  //Sort the list of wedges by the first and, then, second parameters
  //For each unused wedge w=(a,b,c):
  //- Mark w as used
  //- Start a new region with w
  //- Binary search the list to find w2 = (x1,x2,x3)
  //- Assert the binary search does not fail (it shouldn't!)
  //- If x2==a and x3==b (i.e., we returned to w): stop... we finished constructing a polygon...
  //- Go to the first step...

  //sort wedges:
  sort(wedgesTemp.begin(),wedgesTemp.end());
  const int numWedges = wedgesTemp.size();
  usedWedgesTemp.resize(numWedges);
  for(int i=0;i<numWedges;i++) usedWedgesTemp[i] = false;

  //Let's use a ragged array for performance...
  polygons.resize(0);  



  for(int wedgeStart=0;wedgeStart<numWedges;wedgeStart++) {
    if(usedWedgesTemp[wedgeStart]) continue;
    usedWedgesTemp[wedgeStart] = true; //mark w as used...
    polygons.push_back(wedgesTemp[wedgeStart][0]);//first vertex of the polygon...

    pair<int,int> keyForBinarySearch(wedgesTemp[wedgeStart][1],wedgesTemp[wedgeStart][2]);
    while(true) {
      int next = binarySearch(wedgesTemp,keyForBinarySearch);// binary search the list to find w2 = (x1,x2,x3) where (x1,x2) is the keyForBinarySearch
      //assert...
      assert(next!=-1); //this shouldn't happen..

      usedWedgesTemp[next] = true;
      polygons.push_back(wedgesTemp[next][0]);
      keyForBinarySearch.first = wedgesTemp[next][1];
      keyForBinarySearch.second = wedgesTemp[next][2];

      if(wedgesTemp[next][1] == wedgesTemp[wedgeStart][0]) { //have I completed a polygon?
        polygons.push_back(wedgesTemp[wedgeStart][0]);
        break; //let's find a (possible) next vertex to start a new polygon...
      }
    }
  }

  
}

//sort the vertices from meshIdToProcess (or created from intersection) basing on the distance to vertex orig
//tempVars should have at least 4 slots..
void sortVerticesBasingOnDistance(vector<int> &vertices,int origId,int meshIdToProcess,TempVarsRetesselateTriangleFunction &tempVars) {
  const Point &orig = *getPointFromVertexId(origId, meshIdToProcess);

  VertCoord *tempVertCoords = tempVars.tempVertCoords;

  sort(vertices.begin(),vertices.end(),
    [&](const int &v1Id, const int &v2Id) { //is v1 closer to orig than v2?
      const Point &v1 = *getPointFromVertexId(v1Id, meshIdToProcess);
      const Point &v2 = *getPointFromVertexId(v2Id, meshIdToProcess);

      //tempVars[2] will be the distance^2 between v1 and orig..
      tempVertCoords[2] =  0;
      for(int i=0;i<3;i++) {
        tempVertCoords[0] = v1[i];
        tempVertCoords[0] -= orig[i];

        tempVertCoords[1]= tempVertCoords[0];
        tempVertCoords[1]*=tempVertCoords[0];

        tempVertCoords[2]+= tempVertCoords[1];
      }

      //tempVars[3] will be the distance^2 between v2 and orig..
      tempVertCoords[3] =  0;
      for(int i=0;i<3;i++) {
        tempVertCoords[0] = v2[i];
        tempVertCoords[0] -= orig[i];

        tempVertCoords[1] = tempVertCoords[0];
        tempVertCoords[1]*= tempVertCoords[0];

        tempVertCoords[3]+= tempVertCoords[1];
      }
      return tempVertCoords[2]<tempVertCoords[3];

    });

}



//returns the index i of the array v, where v[i][0] = key.first and v[i][1] = key.second 
int binarySearch(array<int,3> v[],int size,int first,int second) {
  int lo = 0, hi = size-1, mid;
  while (lo<=hi) {
    mid = lo + (hi-lo)/2;
    if(v[mid][0]==first && v[mid][1]==second) 
      return mid;
    else if ( (v[mid][0] < first) || ((v[mid][0] == first) && (v[mid][1] < second)) ) //is the element v[mid] smaller than the key?
      lo = mid+1;
    else 
      hi = mid-1;
  }
  return -1;
};

//we suppose the polygons are connected
//after this function all polygons will be oriented in clockwise orientation
void orientPolygonsClockwise(vector<BoundaryPolygon>::iterator itBeginPolygonsToOrient, 
                             vector<BoundaryPolygon>::iterator itEndPolygonsToOrient,pair<int,int> seedEdge,
                             bool isSeedEdgeClockwisedOriented,bool isOrientationOfRetesselatedEqualToTriangle, int meshIdToProcess) {

  const int numPolygons = itEndPolygonsToOrient-itBeginPolygonsToOrient;
  cerr << "Num Polygons: " << numPolygons << endl;
  bool alreadyFound[numPolygons];
  memset(alreadyFound,0,sizeof(alreadyFound));

  queue< pair<int,bool> > polygonsToProcess;


  //find the polygon containing the seed edge and add to queue with correct orientation label
  int polygonContainingSeedEdge = -1;
  for(int id=0;id<numPolygons;id++) {
    BoundaryPolygon &polygon = *(itBeginPolygonsToOrient+id);
    vector<VertexId> &verticesPolygon = polygon.vertexSequence;
    int numVerticesPolygon = verticesPolygon.size()-1; //the last vertex is equal to the first...
    for(int i=0;i<numVerticesPolygon;i++) {
      if(isOrientationOfRetesselatedEqualToTriangle) {
        if(verticesPolygon[i]==seedEdge.first && verticesPolygon[i+1]==seedEdge.second) {
          polygonContainingSeedEdge = id;  
          break;    
        }
      } else
        if(verticesPolygon[i]==seedEdge.second && verticesPolygon[i+1]==seedEdge.first) {
          polygonContainingSeedEdge = id;  
          break;    
        }
    }
  }
  assert(polygonContainingSeedEdge!=-1);
  cerr << "Polygon containing seed: " << polygonContainingSeedEdge << endl;

  if(isOrientationOfRetesselatedEqualToTriangle)
    polygonsToProcess.push(make_pair(polygonContainingSeedEdge,isSeedEdgeClockwisedOriented));
  else
    polygonsToProcess.push(make_pair(polygonContainingSeedEdge,!isSeedEdgeClockwisedOriented));

  cerr << "Seed clockwise: " << isSeedEdgeClockwisedOriented << " retesselated equal: " << isOrientationOfRetesselatedEqualToTriangle << endl;
  cerr << polygonsToProcess.front().second << endl;

  //label as already found..
  alreadyFound[polygonContainingSeedEdge] = true;


  int numEdges = 0;
  for(int id=0;id<numPolygons;id++) { //count number of edges
    BoundaryPolygon &polygon = *(itBeginPolygonsToOrient+id);
    vector<VertexId> &verticesPolygon = polygon.vertexSequence;
    numEdges += verticesPolygon.size()-1;
    cerr << "Size: "<< verticesPolygon.size() << endl;
    cerr << id << " id num vertices -1 ; " << verticesPolygon.size()-1 << endl;
  }

  array<int,3> edges[numEdges]; //for each edge, in what polygon is it?  (edges[i][0]=first vertex,edges[i][1]=second vertex, edges[i][2] = polygon id (starting in 0))

  cerr << "Number of edges: " << numEdges << endl;
  //copy the edges (and polygon id) to the edges array...
  numEdges = 0;
  for(int id=0;id<numPolygons;id++) {
    BoundaryPolygon &polygon = *(itBeginPolygonsToOrient+id);
    vector<VertexId> &verticesPolygon = polygon.vertexSequence;
    int numVerticesPolygon = verticesPolygon.size()-1; //the last vertex is equal to the first...
    //cerr << id << " id num vertices -1 ; " << numVerticesPolygon << endl;
    //cerr << "Curr num edges: " << numEdges << endl;
    for(int i=0;i<numVerticesPolygon;i++) {
      edges[numEdges][0] = verticesPolygon[i];
      edges[numEdges][1] = verticesPolygon[i+1];
      edges[numEdges][2] = id;
      //cerr << edges[numEdges][0] << " " << edges[numEdges][1] << " polygon: " << id << endl;
      numEdges++;
    }
  }
 // cerr << "Number of edges: " << numEdges << endl;

  //sort the edges for binary search...
  sort(edges,edges+numEdges,[](const array<int,3> &a,const array<int,3> &b) { //sort by first element, then by second element..
                                  if(a[0]!=b[0]) return a[0]<b[0]; 
                                  else return a[1]<b[1]; 
                                } ); //for binary search...


  cerr << "Orientation of polygons in set: " << endl;
  for(int id=0;id<numPolygons;id++) {
    BoundaryPolygon &polygon = *(itBeginPolygonsToOrient+id);
    polygon.isInClockwiseDirection(vertices,meshIdToProcess);
  }
  cerr << endl;

  //while queue is not empty
  //remove a polygon, set as processed, add neighbors, re-orient
  while(!polygonsToProcess.empty()) {
    const pair<int,bool> &polygonToProcess = polygonsToProcess.front();
    //cerr << "Polygon to process: " << polygonToProcess.first << endl;
    BoundaryPolygon &polygon = *(itBeginPolygonsToOrient+polygonToProcess.first);
    //cerr << "Get polygon" << endl;

    vector<VertexId> &verticesPolygon = polygon.vertexSequence;
    int numVerticesPolygon = verticesPolygon.size()-1; //the last vertex is equal to the first...
    for(int i=0;i<numVerticesPolygon;i++) {
      int u = verticesPolygon[i];
      int v = verticesPolygon[i+1];


      int vuPolygon = edges[binarySearch(edges,numEdges,v,u)][2] ;//what polygon contains v,u?
      //cerr << "Vu polygon: " << vuPolygon << " " << numPolygons << endl;
      if(vuPolygon>=0 && !alreadyFound[vuPolygon]) { //if vuPolygon <0 --> vuPolygon=-1 --> no polygon contain vu --> this happens for the boundary! (the neighbor of a polygon close to the boundary is the exterior! --> we, obviously, do not re-orient the exterior polygon)
        alreadyFound[vuPolygon] = true;
        polygonsToProcess.push(make_pair(vuPolygon,polygonToProcess.second)); //if this polygon is oriented clockwiselly, the neighbor will be in counterclockwise and vice versa...
      }
    }

    bool clockWise = false;
    //cerr << "Sum before: " << endl;
    clockWise = polygon.isInClockwiseDirection(vertices,meshIdToProcess);
  
    if(!polygonToProcess.second) {
      cerr << "Reverting... to ...\n";
      polygon.reverseVerticesOrder();
      clockWise =polygon.isInClockwiseDirection(vertices,meshIdToProcess);
    }
    cerr << endl;
    assert(clockWise);
    
    polygonsToProcess.pop();    
  }
}





int signCrossProduct2D(const Point &v11,const Point &v12,const Point &v21,const Point &v22,int whatPlaneProjectTo, VertCoord tempVars[]) {
  if(whatPlaneProjectTo == PLANE_X0) {
    VertCoord component1 = (v12[1]-v11[1])*(v22[2]-v21[2]); //v1.y * v2.z 
    VertCoord component2 = (v22[1]-v21[1])*(v12[2]-v11[2]); //v2.y*v1.z
    if(component1==component2) return 0;
    if(component1>component2) return 1;
    return -1;
  }
  else if(whatPlaneProjectTo == PLANE_Y0) { //TODO: xz or zx ?
    VertCoord component1 = (v12[2]-v11[2])*(v22[0]-v21[0]); //v1.z * v2.x 
    VertCoord component2 = (v22[2]-v21[2])*(v12[0]-v11[0]); //v2.z*v1.x
    if(component1==component2) return 0;
    if(component1>component2) return 1;
    return -1;
  }
  VertCoord component1 = (v12[0]-v11[0])*(v22[1]-v21[1]); //v1.x * v2.y 
  VertCoord component2 = (v22[0]-v21[0])*(v12[1]-v11[1]); //v2.x*v1.y
  if(component1==component2) return 0;
  if(component1>component2) return 1;
  return -1;  
}

//we consider a polygon with a straight edge is not convex.. (this is the same definition used in math..)
//example of polygon: 2,-3,5,9,2
//firstElement is the position of the first element in the vector polygons
//lastElement is the position of the last element
bool isConvex(const vector<int> &polygons,int firstElement,int lastElement, int whatPlaneProjectTo, int meshWhereTriangleIs, VertCoord tempVars[]) {
  //what is the sign of the cross product formed by the last vertex (notice that the last element is the first vertex, not the last one since polygons are in the format 2,3,5,7,8,2)
  //the first vertex and the second vertex?
  int signAllVerticesShouldHave = signCrossProduct2D(*getPointFromVertexId(polygons[lastElement-1], meshWhereTriangleIs)
                                    ,*getPointFromVertexId(polygons[firstElement], meshWhereTriangleIs)
                                    ,*getPointFromVertexId(polygons[firstElement], meshWhereTriangleIs)
                                    ,*getPointFromVertexId(polygons[firstElement+1], meshWhereTriangleIs)
                                    ,whatPlaneProjectTo, tempVars);
  if (signAllVerticesShouldHave==0) return false; //we shouldn't have straight edges..
  
  //if all cross products are the same, the polygon is convex..
  //example of polygon: 2,-3,5,9,2
  // now we will process vertices: -3,5 and 9 (in each step we consider the vertex is the center of the wedge)
  for(int i=firstElement+1;i<lastElement;i++) {
    int sign = signCrossProduct2D(*getPointFromVertexId(polygons[i-1], meshWhereTriangleIs)
                                    ,*getPointFromVertexId(polygons[i], meshWhereTriangleIs)
                                    ,*getPointFromVertexId(polygons[i], meshWhereTriangleIs)
                                    ,*getPointFromVertexId(polygons[i+1], meshWhereTriangleIs)
                                    ,whatPlaneProjectTo, tempVars);
    if (sign != signAllVerticesShouldHave) return false;
  }
  return true;
}

int doesPolygonContainAllEdges(const vector<int> &polygons,const int firstElementPolygon,const int lastElementPolygon,
                                const set<pair<int,int> > & edgesFromTriangleBoundary) {

  if(edgesFromTriangleBoundary.count(pair<int,int>(polygons[firstElementPolygon],polygons[firstElementPolygon+1]))) {
    //check if all edges in same order are there...
    for(int i=firstElementPolygon+1;i<lastElementPolygon;i++) {
      if(edgesFromTriangleBoundary.count(pair<int,int>(polygons[i],polygons[i+1]))==0) return 0;
    } 
    return 1;
  }
  else 
    if(edgesFromTriangleBoundary.count(pair<int,int>(polygons[firstElementPolygon+1],polygons[firstElementPolygon]))) {
      //check if all edges in same order are there...
      for(int i=firstElementPolygon+1;i<lastElementPolygon;i++) {
        if(edgesFromTriangleBoundary.count(pair<int,int>(polygons[i+1],polygons[i]))==0) return 0;
      } 
      return -1;
    }
    else return 0;
}




//tempVars should have at least 8 slots
//the new polygons generated from the retesselation will be added to the back of "newPolygonsGeneratedFromRetesselation"
void retesselateTriangleUsingWedgeSorting(const vector<pair<int,int> > &edgesUsingVertexId, const vector<int> &edgesFromIntersection,const Triangle &t, const int meshWhereTriangleIs, vector<BoundaryPolygon> &newPolygonsGeneratedFromRetesselation,TempVarsRetesselateTriangleFunction &tempVars, StatisticsAboutRetesseation &statistics) {
  //A triangle (except if it is degenerate to a single point) cannot be perpendicular to all 3 planes (x=0,y=0,z=0)
  //Thus, we chose one of these planes to project the triangle during computations...
  const int whatPlaneProjectTriangleTo = getPlaneTriangleIsNotPerpendicular(vertices[meshWhereTriangleIs][t.p[0]], vertices[meshWhereTriangleIs][t.p[1]],vertices[meshWhereTriangleIs][t.p[2]], tempVars.tempVertCoords);

  /*cerr << "Number of edges from intersection: " << edgesFromIntersection.size() << endl;
  for(int e: edgesFromIntersection) {
    Point &a = *getPointFromVertexId(edgesUsingVertexId[e].first,meshWhereTriangleIs);
    Point &b = *getPointFromVertexId(edgesUsingVertexId[e].second,meshWhereTriangleIs);
    cerr << edgesUsingVertexId[e].first << " ( " << a[0].get_d() << " , " << a[1].get_d() << " ) ";
    cerr << edgesUsingVertexId[e].second << " ( " << b[0].get_d() << " , " << b[1].get_d() << " ) " << endl;
  }*/

  /*
  //First we need to determine, for each vertex v, the edges incident in v
  //each vertex will be represented twice (as two directed edges): (u,v) and (v,u)
  //The edges incident in each vertex will be sorted based on their slope (supposing they are projected to the plane "whatPlaneProjectTriangleTo")
  map<int, vector<int>> edgesIncidentToEachVertex; // if v is in edgesIncidentToEachVertex[u] --> there is an edge (u,v)



  //extractPolygonsFromEdgeList(const vector<pair<int,int> > &edgesUsingVertexId,const int meshIdToProcess, const int whatPlaneProjectTo, vector<int> &polygons,vector<array<int,3> > &wedgesTemp, vector<pair<int,int> > &edgesTemp, vector<bool> &usedWedgesTemp, VertCoord tempCoords[2] )

  */

  //we need to "retesselate" t and orient all the new triangles properly
  //this set will store the edges we used in the retesselated triangle
  //we need to choose what edges to create and, then, use these new edges to reconstruct the triangulation..
  //This set will have the edges that will form the retriangulation of ts..
  vector<pair<int,int> > edgesUsedInThisTriangle; 

  int ct =0;
  //The edges from intersection will, necessarelly, be in the triangulation...
  for(int edgeId:edgesFromIntersection) {
    //cerr << edgesUsingVertexId[edgeId].first << " " << edgesUsingVertexId[edgeId].second << endl;
    assert(edgesUsingVertexId[edgeId].first != edgesUsingVertexId[edgeId].second);
    edgesUsedInThisTriangle.push_back(edgesUsingVertexId[edgeId]);   

    #ifdef COLLECT_STATISTICS
      #pragma omp atomic
      statistics.ctEdgesInsertedBecauseOfConstraint++;   
    #endif
  }


  //Now, we need to add more edges to fill the parts of ts that still do not form triangle
  vector<int> verticesToTesselate;
  //what vertices will we have in the new triangulation of ts? (the original vertices + the vertices of the edges from intersections)
  for(const auto &p:edgesUsedInThisTriangle) {
    verticesToTesselate.push_back(p.first);
    verticesToTesselate.push_back(p.second);
  }
  verticesToTesselate.push_back(t.p[0]);
  verticesToTesselate.push_back(t.p[1]);
  verticesToTesselate.push_back(t.p[2]);
  sort(verticesToTesselate.begin(),verticesToTesselate.end());
  auto newEnd = unique(verticesToTesselate.begin(),verticesToTesselate.end());
  verticesToTesselate.resize(newEnd-verticesToTesselate.begin());
  int numVTriangle = verticesToTesselate.size();

  #ifdef COLLECT_STATISTICS
    #pragma omp atomic
    statistics.ctVerticesInsertedTrianglesToRetesselate += numVTriangle;
  #endif

  

  const Point &triangleVertex0 = vertices[meshWhereTriangleIs][t.p[0]];
  const Point &triangleVertex1 = vertices[meshWhereTriangleIs][t.p[1]];
  const Point &triangleVertex2 = vertices[meshWhereTriangleIs][t.p[2]];
  //cerr << triangleVertex0[0].get_d( ) << " " << triangleVertex0[1].get_d() << " " << triangleVertex0[2].get_d() << endl;;
  //cerr << triangleVertex1[0].get_d( ) << " " << triangleVertex1[1].get_d() << " " << triangleVertex1[2].get_d() << endl;;
  //cerr << triangleVertex2[0].get_d( ) << " " << triangleVertex2[1].get_d() << " " << triangleVertex2[2].get_d() << endl;;
  //cerr << "Project to: " << whatPlaneProjectTriangleTo << endl;


  //Computes the vertices incidnt to each edge of the original triangle...
  //verticesIncidentEachEdgeOriginalTriangle[0] is for vertex t.p[0]-t.p[1]
  //verticesIncidentEachEdgeOriginalTriangle[1] is for vertex t.p[1]-t.p[2]
  //verticesIncidentEachEdgeOriginalTriangle[2] is for vertex t.p[2]-t.p[0]
  vector<int> verticesIncidentEachEdgeOriginalTriangle[3];
  for(int i=0;i<numVTriangle;i++) {
    int v = verticesToTesselate[i];
    if(v==t.p[0] || v==t.p[1] || v==t.p[2]) continue;
    //for each vertex v that is not one of the original vertices of the triangle, let's see if v is in one of the edges of t...

    const Point &vertex = *getPointFromVertexId(v, meshWhereTriangleIs);
    //cerr << "V: " << vertex[0].get_d( ) << " " << vertex[1].get_d() << " " << vertex[2].get_d() << endl;;
    int o1 = orientation(triangleVertex0, triangleVertex1,  vertex,whatPlaneProjectTriangleTo,tempVars.tempVertCoords);
    
    
    //cerr << o1 << o2 << o3 << endl;
    if(o1==0) {
      verticesIncidentEachEdgeOriginalTriangle[0].push_back(v);
      continue;
    }
    int o2 = orientation(triangleVertex1, triangleVertex2,  vertex,whatPlaneProjectTriangleTo,tempVars.tempVertCoords);
    if(o2==0) {
      //cerr << "12" << endl;
      verticesIncidentEachEdgeOriginalTriangle[1].push_back(v);
      continue;
    }
    int o3 = orientation(triangleVertex2, triangleVertex0,  vertex,whatPlaneProjectTriangleTo,tempVars.tempVertCoords);
    if(o3==0) {
     // cerr << 20 << endl;
      verticesIncidentEachEdgeOriginalTriangle[2].push_back(v);
      continue;
    }
  }
  //cerr << "In edges: " << numVerticesInEdges << "\n";
 // cerr << "Vertices: " << numVTriangle << "\n";
  //cerr << "Constraint edges: " << edgesFromIntersection.size() << "\n\n";

  set<pair<int,int> > edgesFromTriangleBoundary;

  bool allSimple = true;
  //no vertex intersect the edge t.p[0]-t.p[1] --> it will be in the triangulation!!!
  int verticesIncidentEdge;
  verticesIncidentEdge = verticesIncidentEachEdgeOriginalTriangle[0].size();
  if(verticesIncidentEdge==0) {
    int v1 = t.p[0];
    int v2 = t.p[1];
            
    pair<int,int> candidateEdge(v1,v2);  

    #ifdef COLLECT_STATISTICS
      #pragma omp atomic
      statistics.ctEdgesActuallyInsertedInRetesselation++;
    #endif

    assert(candidateEdge.first != candidateEdge.second);
    edgesUsedInThisTriangle.push_back(candidateEdge);
    edgesFromTriangleBoundary.insert(candidateEdge);
  } else {
    if(verticesIncidentEdge==1) {
      int v1 = t.p[0];
      int v2 = t.p[1]; 
      int v =  *(verticesIncidentEachEdgeOriginalTriangle[0].begin());          
      
      pair<int,int> candidateEdge;    

      candidateEdge.first = v1; 
      candidateEdge.second = v;


      #ifdef COLLECT_STATISTICS
        #pragma omp atomic
        statistics.ctEdgesActuallyInsertedInRetesselation++;
      #endif

      assert(candidateEdge.first != candidateEdge.second);
      edgesUsedInThisTriangle.push_back(candidateEdge);
      edgesFromTriangleBoundary.insert(candidateEdge);

      candidateEdge.first = v; 
      candidateEdge.second = v2;
      

      #ifdef COLLECT_STATISTICS
        #pragma omp atomic
        statistics.ctEdgesActuallyInsertedInRetesselation++;
      #endif
      assert(candidateEdge.first != candidateEdge.second);
      edgesUsedInThisTriangle.push_back(candidateEdge);
      edgesFromTriangleBoundary.insert(candidateEdge);
    } else {
      allSimple= false;
    }
  }

  verticesIncidentEdge = verticesIncidentEachEdgeOriginalTriangle[1].size();
  if(verticesIncidentEdge==0) {
    int v1 = t.p[1];
    int v2 = t.p[2];
            
    pair<int,int> candidateEdge(v1,v2);    



    #ifdef COLLECT_STATISTICS
        #pragma omp atomic
        statistics.ctEdgesActuallyInsertedInRetesselation++;
    #endif

    assert(candidateEdge.first != candidateEdge.second);
    edgesUsedInThisTriangle.push_back(candidateEdge);
    edgesFromTriangleBoundary.insert(candidateEdge);
  } else {
    if(verticesIncidentEdge==1) {
      int v1 = t.p[1];
      int v2 = t.p[2]; 
      int v =  *(verticesIncidentEachEdgeOriginalTriangle[1].begin());          
      
      pair<int,int> candidateEdge;    

      candidateEdge.first = v1; 
      candidateEdge.second = v;

      #ifdef COLLECT_STATISTICS
        #pragma omp atomic
        statistics.ctEdgesActuallyInsertedInRetesselation++;
      #endif
      assert(candidateEdge.first != candidateEdge.second);
      edgesUsedInThisTriangle.push_back(candidateEdge);
      edgesFromTriangleBoundary.insert(candidateEdge);

      candidateEdge.first = v; 
      candidateEdge.second = v2;

      #ifdef COLLECT_STATISTICS
        #pragma omp atomic
        statistics.ctEdgesActuallyInsertedInRetesselation++;
      #endif
      assert(candidateEdge.first != candidateEdge.second);
      edgesUsedInThisTriangle.push_back(candidateEdge);
      edgesFromTriangleBoundary.insert(candidateEdge);
    } else {
      allSimple= false;
    }
  }

  verticesIncidentEdge = verticesIncidentEachEdgeOriginalTriangle[2].size();
  if(verticesIncidentEdge==0) {
    int v1 = t.p[2];
    int v2 = t.p[0];
            
    pair<int,int> candidateEdge(v1,v2);

    #ifdef COLLECT_STATISTICS
      #pragma omp atomic
      statistics.ctEdgesActuallyInsertedInRetesselation++;
    #endif
    assert(candidateEdge.first != candidateEdge.second);
    edgesUsedInThisTriangle.push_back(candidateEdge);
    edgesFromTriangleBoundary.insert(candidateEdge);
  } else {
    if(verticesIncidentEdge==1) {
      int v1 = t.p[2];
      int v2 = t.p[0]; 
      int v =  *(verticesIncidentEachEdgeOriginalTriangle[2].begin());          
      
      pair<int,int> candidateEdge;    

      candidateEdge.first = v1; 
      candidateEdge.second = v;
      #ifdef COLLECT_STATISTICS
        #pragma omp atomic
        statistics.ctEdgesActuallyInsertedInRetesselation++;
      #endif
      assert(candidateEdge.first != candidateEdge.second);
      edgesUsedInThisTriangle.push_back(candidateEdge);
      edgesFromTriangleBoundary.insert(candidateEdge);

      candidateEdge.first = v; 
      candidateEdge.second = v2;

      #ifdef COLLECT_STATISTICS
        #pragma omp atomic
        statistics.ctEdgesActuallyInsertedInRetesselation++;
      #endif
      assert(candidateEdge.first != candidateEdge.second);
      edgesUsedInThisTriangle.push_back(candidateEdge);
      edgesFromTriangleBoundary.insert(candidateEdge);
    } else {
      allSimple= false;
    }
  }
  
  #ifdef COLLECT_STATISTICS        
      if(allSimple) {
        #pragma omp atomic
        statistics.simple++;
      }
      else {
        #pragma omp atomic
        statistics.difficult++;
      }
  #endif
  //Now we have a set of internal edges and possibly some boundary edges
  //let's split the boundary edges (with more than 2 vertices) at the intersection points
  for(int i=0;i<3;i++) {
    int v1 = t.p[(i)]; //the two vertices of the edge i...
    int v2 = t.p[(i+1)%3]; 

    vector<int> &verticesToCreateEdges = verticesIncidentEachEdgeOriginalTriangle[i];
    //cerr << "Num of vertices in edge: " << i << " : " << verticesToCreateEdges.size() << endl;
    if(verticesToCreateEdges.size()<=1) continue; //edges were already inserted in the previous step...
    //here we need to deal only with the situation where 2 or more vertices are incident to an edge.

    //cerr << "Ids of candidate vertices: " << endl;
    //for(int v: verticesToCreateEdges) cerr << v << endl;

    //sort the vertices basing on the distance from v1
    sortVerticesBasingOnDistance(verticesToCreateEdges,v1,meshWhereTriangleIs,tempVars);
    
    verticesToCreateEdges.push_back(v2);
    //we will add edges:
    //v1 to verticesToCreateEdges[0]
    //verticesToCreateEdges[0] to verticesToCreateEdges[1]
    //....... verticesToCreateEdges[n-2] to verticesToCreateEdges[n-1] (=v2, that was pushed_back above..)

    pair<int,int> candidateEdge;   
    candidateEdge.first = v1; 
    candidateEdge.second = verticesToCreateEdges[0];
    #ifdef COLLECT_STATISTICS
      #pragma omp atomic
      statistics.ctEdgesActuallyInsertedInRetesselation++;
    #endif
    //cerr << "Candidate edge: " << candidateEdge.first << " " << candidateEdge.second << endl;
    //cerr <<"V1 v2 : " << v1 << " " << v2 << endl;
    assert(candidateEdge.first != candidateEdge.second);
    //assert(edgesUsedInThisTriangle.count(candidateEdge)==0);
    //cerr << "Inserting edge: " << candidateEdge.first << " " << candidateEdge.second << endl;
    //printVertexForDebugging(*getPointFromVertexId(candidateEdge.first,meshWhereTriangleIs));
    //printVertexForDebugging(*getPointFromVertexId(candidateEdge.second,meshWhereTriangleIs));

    edgesUsedInThisTriangle.push_back(candidateEdge);
    edgesFromTriangleBoundary.insert(candidateEdge);

    const int numVerticesToCreateEdges = verticesToCreateEdges.size();
    for(int i=0;i<numVerticesToCreateEdges-1;i++) {
      candidateEdge.first = verticesToCreateEdges[i]; 
      candidateEdge.second = verticesToCreateEdges[i+1];      
      #ifdef COLLECT_STATISTICS
        #pragma omp atomic
        statistics.ctEdgesActuallyInsertedInRetesselation++;
      #endif
      assert(candidateEdge.first != candidateEdge.second);
      //assert(edgesUsedInThisTriangle.count(candidateEdge)==0);

      //cerr << "Inserting edge: " << candidateEdge.first << " " << candidateEdge.second << endl;
      //printVertexForDebugging(*getPointFromVertexId(candidateEdge.first,meshWhereTriangleIs));
      //printVertexForDebugging(*getPointFromVertexId(candidateEdge.second,meshWhereTriangleIs));
      edgesUsedInThisTriangle.push_back(candidateEdge);
      edgesFromTriangleBoundary.insert(candidateEdge);
    }
  }

  // (at least) three of the last edges added to the triangle are, necessarely, in the boundary of the triangle
  //they also have the same orientation of the triangle.
  //so, we can use any of them as the "seed edge"
  //we need a seed output edge that will be oriented in the same way t is oriented..
  //we will use this seed to reorient all the triangles resulting from the retesselation of t
  pair<int,int> seedEdge = edgesUsedInThisTriangle.back(); 

  //Now, edgesUsedInThisTriangle will contain the original edges of the triangle (split because of intersections) and new edges
  //created because of intersections with other triangles. Thus, edgesUsedInThisTriangle will define a planar partitioning..

 // cerr << "End..." << allSimple << endl;
  
  vector<pair<int,int> > edgesUsedInThisTriangleV(edgesUsedInThisTriangle.begin(),edgesUsedInThisTriangle.end());
  
  sortEdgesAndExtractPolygonsFromEdgeListUsingWedges(edgesUsedInThisTriangleV,edgesFromIntersection.size(),meshWhereTriangleIs,whatPlaneProjectTriangleTo,   tempVars );
  vector<int> &polygons = tempVars.polygons; 


  //at the end of the previous function raggedArrayEdges will contain the directed edges sorted by their first vertex...
  bool arePolygonsConnected = checkIfPolygonsAreConnected(tempVars.raggedArraySortedEdges);

  //at the end of the previous function, raggedArrayEdges will be filled with the edges sorted by the first vertex of each edge
  //(if there is a tie the edges are sorted by the slope w.r.t. the plane we projected the triangle)
  int numConvexPolygonsFound=0;
  int numConcavePolygonsFound = 0;

  #ifdef COLLECT_STATISTICS
  if(arePolygonsConnected) {
    #pragma omp atomic
    statistics.numConnectedPolygonSubdivisionOfTriangles++;
  } else {
    #pragma omp atomic
    statistics.numDisconnectedPolygonSubdivisionOfTriangles++;
  }

  #endif

  const int sizePolygonsVector = polygons.size();

  int numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes = newPolygonsGeneratedFromRetesselation.size();
  
  if(arePolygonsConnected) {
    //since the graph is connected, all polygons in "polygons" are oriented similarly
    //we need to check if this orientation is equal to the orientation of the original triangle
    //if it is not, all the polygons will have an orientation that is contrary to the triangle's orientation

    //Let's copy the polygons (except the exterior polygon) to the newPolygonsGeneratedFromRetesselation vector 
    //once we find the exterior polygon, we check if the orientation is correct

    /*cerr << "Number of vertices in polygons in triangle (including exterior): " << sizePolygonsVector << endl;
    cerr << "Triangle: " << t.p[0] << " " << t.p[1] << " " << t.p[2] << endl;
    for(int v:polygons) {
      Point &a = *getPointFromVertexId(v,meshWhereTriangleIs);
      cerr << v << " ( " << a[0].get_d() << " , " << a[1].get_d() << " ) \n";
    }
    cout << endl;
    */

    bool isOrientationOfRetesselatedEqualToTriangle = true;
    //newPolygonsGeneratedFromRetesselation.reserve(numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes+sizePolygonsVector-1);

    bool exteriorPolygonAlreadyFound = false; //have we already found the exterior polygon?

    //for each polygon p inside the triangle, we will triangulate p
    for(int firstElementPolygon=0;firstElementPolygon<sizePolygonsVector;) { 
      int lastElementPolygon = firstElementPolygon+1;
      while(polygons[lastElementPolygon]!= polygons[firstElementPolygon]) lastElementPolygon++;

      //let's process the polygon defined by the vertices in polygons[firstElementPolygon..lastElementPolygon]

      //First, let's check if this polygon is the exterior polygon (since all polygons are connected,
      //it will be iff the three vertices of the original triangle are in the polygon)
      if(!exteriorPolygonAlreadyFound) {
        int numVerticesOriginalTriangleInThisPolygon = 0;
        for(int i=firstElementPolygon;i<lastElementPolygon;i++) {
          if(polygons[i]==t.p[0] || polygons[i]==t.p[1] || polygons[i]==t.p[2]) numVerticesOriginalTriangleInThisPolygon++;
        }
        if(numVerticesOriginalTriangleInThisPolygon>=3) {          
          //this will be 0 if false, 1 if has all edges and in same orientation, -1 if has all edges but in reverse orientation..
          int polygonContainAllEdges = doesPolygonContainAllEdges(polygons,firstElementPolygon,lastElementPolygon,edgesFromTriangleBoundary);
          if(polygonContainAllEdges!=0) { //we found the boundary!
            exteriorPolygonAlreadyFound = true;
            isOrientationOfRetesselatedEqualToTriangle = polygonContainAllEdges==-1; //if the exterior triangle is oriented in the countrary way --> the polygons will be oriented similarly to the triangle
            firstElementPolygon = lastElementPolygon+1;
            continue; //we will not add the exterior polygon to the list of polygons in retesselation...
          }     
        }
      }

      //let's copy the polygon to the list of polygons from retesselation...
      newPolygonsGeneratedFromRetesselation.push_back(BoundaryPolygon(whatPlaneProjectTriangleTo));
      BoundaryPolygon &polygon = newPolygonsGeneratedFromRetesselation.back();
      polygon.vertexSequence.assign(polygons.begin()+firstElementPolygon,polygons.begin()+lastElementPolygon+1);
      polygon.whatPlaneProjectTriangleTo = whatPlaneProjectTriangleTo;

      firstElementPolygon = lastElementPolygon+1;
    }

    int numberPolygonsFromRetesselationNow = newPolygonsGeneratedFromRetesselation.size();
    //let's assign the polygon orientation...
    assert(exteriorPolygonAlreadyFound);
    if(isOrientationOfRetesselatedEqualToTriangle)
      for(int i=numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes;i<numberPolygonsFromRetesselationNow;i++) { //the polygons we are adding are only in this interval (the vector may have other polygons before we added new ones...)
        newPolygonsGeneratedFromRetesselation[i].above = t.above;
        newPolygonsGeneratedFromRetesselation[i].below = t.below;
      }
    else 
      for(int i=numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes;i<numberPolygonsFromRetesselationNow;i++) {
        newPolygonsGeneratedFromRetesselation[i].above = t.below;
        newPolygonsGeneratedFromRetesselation[i].below = t.above;
      }

    //if the cross product of the vectors: (vertex1-vertex0) and (vertex2-vertex0) is negative --> the triangle is oriented
    //in clockwise orientation --> the seed is also in clockwise orientation...
    bool isBoundaryOriginalTriangleClockwisedOriented = signCrossProduct2D(triangleVertex0,triangleVertex1,
                                                            triangleVertex0,triangleVertex2,whatPlaneProjectTriangleTo,tempVars.tempVertCoords)<0;

    bool areInteriorPolygonsClockwiseOriented = true;
    if(isBoundaryOriginalTriangleClockwisedOriented && !isOrientationOfRetesselatedEqualToTriangle) areInteriorPolygonsClockwiseOriented = false;
    if(!isBoundaryOriginalTriangleClockwisedOriented && isOrientationOfRetesselatedEqualToTriangle) areInteriorPolygonsClockwiseOriented = false;

   // cerr << "Polygons before and after: " << numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes << " " << newPolygonsGeneratedFromRetesselation.size() << endl;
    //cerr << "Orienting polygon.." << endl;
    //Now, we will reorient the polygons such that all polygons will be oriented clockwise
    //orientPolygonsClockwise(newPolygonsGeneratedFromRetesselation.begin()+numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes, 
   //                        newPolygonsGeneratedFromRetesselation.end(),seedEdge,isSeedEdgeClockwisedOriented,
    //                        isOrientationOfRetesselatedEqualToTriangle,meshWhereTriangleIs);

    if(!areInteriorPolygonsClockwiseOriented) {
      //we need to reverse the orientation of all polygons...
      for(int i=numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes;i<numberPolygonsFromRetesselationNow;i++) 
        newPolygonsGeneratedFromRetesselation[i].reverseVerticesOrder();
    }

    #ifdef SANITY_CHECKS

    
    for(int i=numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes;i<numberPolygonsFromRetesselationNow;i++) {
      bool isOrientedCorrectly = newPolygonsGeneratedFromRetesselation[i].isInClockwiseDirection(vertices, meshWhereTriangleIs);
      if(!isOrientedCorrectly) {
        #pragma omp critical
          {
          cerr << "Error... polygon " << i << " generated from triangle "  << " is not oriented correctly..\n";
          cerr << numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes << " " << numberPolygonsFromRetesselationNow << endl;
          cerr << "isBoundaryOrigintalTriOrientedClockwise " << isBoundaryOriginalTriangleClockwisedOriented << endl;
          cerr << "isOrientationOfRetesselatedEqualToTriangle " << isOrientationOfRetesselatedEqualToTriangle << endl;
          cerr << "Original triangle: " << endl;
          cerr << t.p[0] << " " << t.p[1] << " " << t.p[2] << endl;
          printVertexForDebugging(getPointFromVertexId(t.p[0],meshWhereTriangleIs)->data());
          printVertexForDebugging(getPointFromVertexId(t.p[1],meshWhereTriangleIs)->data());
          printVertexForDebugging(getPointFromVertexId(t.p[2],meshWhereTriangleIs)->data());
          cerr << "Polygon: " << endl;
          for(int v:newPolygonsGeneratedFromRetesselation[i].vertexSequence) cerr << v << " "; cerr << endl;
          for(int v:newPolygonsGeneratedFromRetesselation[i].vertexSequence) printVertexForDebugging(getPointFromVertexId(v,meshWhereTriangleIs)->data());
          cerr << "What plane: " << whatPlaneProjectTriangleTo << endl;
          cerr << "Edges from planar graph " << endl;
          for(auto a:edgesUsedInThisTriangleV) {
            cerr << a.first << " " << a.second << endl;
            printVertexForDebugging(getPointFromVertexId(a.first,meshWhereTriangleIs)->data());
            printVertexForDebugging(getPointFromVertexId(a.second,meshWhereTriangleIs)->data());
            cerr << endl;
          }
          cerr << "Polygons extracted " << endl;
          for(int j=numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes;j<numberPolygonsFromRetesselationNow;j++) {
            for(int v:newPolygonsGeneratedFromRetesselation[j].vertexSequence) cerr << v << " "; cerr << endl;
            for(int v:newPolygonsGeneratedFromRetesselation[j].vertexSequence) printVertexForDebugging(getPointFromVertexId(v,meshWhereTriangleIs)->data());
          }
        }
      }
      assert(isOrientedCorrectly);
    }
    #endif
   // cerr << "Oriented \n";
  } else {
    //TODO
    //TODO: remember to set whatPlaneProjectTriangleTo of the polygon...

    //let's print the triangles with disconnected polygons in its interior for debugging purposes...
    //edgesUsedInThisTriangle
    cerr << "Disconnected triangle:::: this situation is not treated yet..." << endl << endl << endl;
    int numDisconnectedTrianglesNow = 0;
    numDisconnectedTrianglesNow = statistics.numDisconnectedPolygonSubdivisionOfTriangles;
    stringstream trianglesDisconnectedPath;
    trianglesDisconnectedPath << "triangleDisconnectedPath_" << numDisconnectedTrianglesNow << ".gts";
    string path = trianglesDisconnectedPath.str();
    saveEdgesAsGTS(edgesUsedInThisTriangle, meshWhereTriangleIs,  path);


  }

  #ifdef COLLECT_STATISTICS
    #pragma omp atomic
    statistics.ctConvexPolygonsInTriangleRetesselations += numConvexPolygonsFound;
    #pragma omp atomic
    statistics.ctConcavePolygonsInTriangleRetesselations += numConcavePolygonsFound;
  #endif

  
}





//edges represent the edges generated from the intersection of the intersecting triangles.
//intersectingTrianglesThatGeneratedEdges[i] contains the pair of triangles whose intersection generated edges[i]

void retesselateIntersectingTriangles(const vector< pair< array<VertCoord,3>,array<VertCoord,3> > > &edges, const vector< pair<Triangle *,Triangle *> > &intersectingTrianglesThatGeneratedEdges) {
	assert(edges.size()==intersectingTrianglesThatGeneratedEdges.size());
	timespec t0,t1;
  cerr << "Retesselating triangles..." << "\n";
  clock_gettime(CLOCK_REALTIME, &t0);

  /*long long ctTrianglesRetesselate = 0;
  long long ctVerticesInsertedTrianglesToRetesselate = 0;
  long long ctEdgesTestedToInsertInRetesselation = 0;
  long long ctEdgesActuallyInsertedInRetesselation = 0; //besides the constraint edges (the ones from intersection)
  long long ctEdgesInsertedBecauseOfConstraint = 0; //constraint edges inserted in retesselation*/

  StatisticsAboutRetesseation statisticsAboutRetesseation;

	unordered_map<const Triangle *, vector<int> > intersectingEdgesInEachTriangle[2];
  //trianglesIntersectingEachTriangle[0] --> for each triangle t from map 0 (that intersects other triangles), 
  //trianglesIntersectingEachTriangle[0][t] is a vector of triangles from mesh 1 intersecting t...

  map< Point, int > vertexIdPlusOne; //lets add 1 to the ids (this will simplify the creation of new entries...)
  
  int numEdges = edges.size();

  clog << "Number of edges: " << numEdges << "\n";
  vector<pair<int,int> > edgesUsingVertexId(numEdges); //if we have an edge (a,b)  --> a<b
  for(int i=0;i<numEdges;i++) {
  	//cerr << "i " << i << " " << numEdges << endl;
  	const auto &e = edges[i];
  	int idA = vertexIdPlusOne[e.first];
  	if(idA==0) {
  		idA = vertexIdPlusOne[e.first] = vertices[2].size()+1; //when we create a new vertex we start counting from 1 (only here!!!)
  		vertices[2].push_back(e.first);
  	}
  	int idB = vertexIdPlusOne[e.second];
  	if(idB==0) {
  		idB = vertexIdPlusOne[e.second] = vertices[2].size()+1;  //when we create a new vertex we start counting from 1 (only here!!!)
  		vertices[2].push_back(e.second);
  	}
  	idA = -idA; //for simplicity, let's use a negative id to referrence a "shared" vertex (from layer 2)
  	idB = -idB;
  	if(idA > idB) swap(idA,idB); 

  	edgesUsingVertexId[i] = make_pair(idA,idB); 

  	//the intersection of tA (from mesh 0) with tB generated the i-th edge...
  	const Triangle *tA = intersectingTrianglesThatGeneratedEdges[i].first;
  	const Triangle *tB = intersectingTrianglesThatGeneratedEdges[i].second;

  	intersectingEdgesInEachTriangle[0][tA].push_back(i);
  	intersectingEdgesInEachTriangle[1][tB].push_back(i);
  }
  cerr << "Number of edges pushed back: " << vertices[2].size() << endl;
  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to extract the edges intersecting each triangle and create the new vertices: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 
 

  //edgesUsingVertexId : each edge contains the ids of its two vertices (idA,idB)
  //idA,idB start with 0 and they represent the position of the vertex in vertices[2].
  //i.e., if an edge is = (idA,idB) --> the two vertices that form this edge are vertices[2][idA],vertices[2][idB]

  //for each triangle t from mesh meshId, intersectingEdgesInEachTriangle[meshId][t] contains a vector of ids of edges formed by the intersection
  //of t with other triangles...


  //TODO: special cases, collinear, etc..

  vector<pair<int,int> > newTriEdgesFromEachMap[2];

  for(int meshIdToProcess=0;meshIdToProcess<2;meshIdToProcess++) {
  	vector< const pair<const Triangle * const, vector<int> >* > intersectingEdgesInEachTriangleToProcess;
  	for(const auto &elem:intersectingEdgesInEachTriangle[meshIdToProcess]) intersectingEdgesInEachTriangleToProcess.push_back(&elem);

  	int numTrianglesToProcess = intersectingEdgesInEachTriangleToProcess.size();

  	cerr << "Number of triangles to retesselate in this mesh: " << numTrianglesToProcess << "\n";

    statisticsAboutRetesseation.ctTrianglesRetesselate += numTrianglesToProcess;

  	#pragma omp parallel 
  	{
  		//vector<pair<int,int> > myNewTriEdgesFromEachMap;

  		vector<BoundaryPolygon> myNewPolygonsFromRetesselation;
      //myNewPolygonsFromRetesselation.reserve(numTrianglesToProcess);

  		TempVarsRetesselateTriangleFunction tempVars;

      const int percentShowLog = 1;
      int onePercentNumTriangles = (percentShowLog*numTrianglesToProcess)/100;
      if(onePercentNumTriangles==0) onePercentNumTriangles = 1;

  		#pragma omp for
		  for(int i=0;i<numTrianglesToProcess;i++) {  		
        if((i%onePercentNumTriangles)==0) clog << "Retesselating " << i << " of " << numTrianglesToProcess << " Percent= " << (i*100)/numTrianglesToProcess << "\n";  		
        //for each triangle ts, let's process the triangles intersecting t and the edges
        //formed by the intersection of ts and these triangles

		  	const auto &ts  = *intersectingEdgesInEachTriangleToProcess[i];
		  	//ts.first = a triangle from map 0
		  	//ts.second = list of edges formed by the intersection of ts.first with other triangles...
		  	const auto &t = *ts.first;
		  	const auto &edgesFromIntersection =  ts.second;

        
        #ifdef COLLECT_STATISTICS_PRINT_TRIANGLES_INTERSECTIONS
          #pragma omp critical
          //cout << "Triangle from mesh " << meshIdToProcess << " intersects: " << ts.second.size() << "\n";
          if(ts.second.size()==6) {
            #pragma omp critical
            printTriangleWithManyIntersection(t,meshIdToProcess,ts.second,intersectingTrianglesThatGeneratedEdges);
          }
          if(ts.second.size()==6) {
            #pragma omp critical
            saveEdgesAsGTS(edges,ts.second,"edgesFromTriangleIntersetingManyTriangles6.gts");
          }
         if(meshIdToProcess==0) {
            #pragma omp critical
            saveEdgesAsGTS(edges,ts.second,"edgesFromTriangleIntersetingMesh0.gts");
          }
        #endif

        //retesselateTriangle(edgesUsingVertexId, edgesFromIntersection,t, meshIdToProcess, myNewTrianglesFromRetesselation, tempVars, statisticsAboutRetesseation);
		    retesselateTriangleUsingWedgeSorting(edgesUsingVertexId, edgesFromIntersection,t, meshIdToProcess, myNewPolygonsFromRetesselation, tempVars, statisticsAboutRetesseation);

      }

		  #pragma omp critical 
		  {
        cerr << "Number of polygons in vector of retesselated: " << myNewPolygonsFromRetesselation.size() << endl;
        cerr << "Capacity: " << myNewPolygonsFromRetesselation.capacity() << endl;
		  	//newTriEdgesFromEachMap[meshIdToProcess].insert( newTriEdgesFromEachMap[meshIdToProcess].end(),myNewTriEdgesFromEachMap.begin(),myNewTriEdgesFromEachMap.end());
				//trianglesFromRetesselation[meshIdToProcess].insert(trianglesFromRetesselation[meshIdToProcess].end(),myNewTrianglesFromRetesselation.begin(),myNewTrianglesFromRetesselation.end());
			  polygonsFromRetesselation[meshIdToProcess].insert(polygonsFromRetesselation[meshIdToProcess].end(),myNewPolygonsFromRetesselation.begin(),myNewPolygonsFromRetesselation.end() );
      }

			//break;
		}

    
    for(BoundaryPolygon p:polygonsFromRetesselation[meshIdToProcess]) {
      //assert(p.isInClockwiseDirection(vertices,meshIdToProcess));
    }
	}


  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to retesselate creating new polygons: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n";  

  clock_gettime(CLOCK_REALTIME, &t0);
  cerr << "Triangulating polygons from retesselation...\n";
  for(int meshIdToProcess=0;meshIdToProcess<2;meshIdToProcess++) {
    int numPolygons = polygonsFromRetesselation[meshIdToProcess].size();

    const int percentShowLog = 10;
    int onePercentNumPolygons = (numPolygons*percentShowLog)/100;
    if(onePercentNumPolygons==0) onePercentNumPolygons = 1;

    #pragma omp parallel
    {
      VertCoord tempCoords[2];
      #pragma omp for
      for(int i=0;i<numPolygons;i++) {
        //cerr << "Triangulating " << i << " of " << numPolygons << " Percent= " << i*100/numPolygons << endl;
        if((i%onePercentNumPolygons)==0) clog << "Triangulating " << i << " of " << numPolygons << " Percent= " << (i*100)/numPolygons << "\n";
        polygonsFromRetesselation[meshIdToProcess][i].triangulatePolygon(vertices,meshIdToProcess,tempCoords);
      }
    }
  }

  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to triangulate polygons: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 

  
  cerr << "Counts computed during retesselation: " << "\n";
  statisticsAboutRetesseation.printStatsSoFar();
  cerr << "Number of inters. tests (for insertin tris.)that are true    : " << ctEdgeIntersect << "\n";
  cerr << "Number of inters. tests (for insertin tris.)that are false    : " << ctEdgeDoNotIntersect << "\n";
  
  cerr << "Number of triangles created from retesselation of tris from mesh 0: " << trianglesFromRetesselation[0].size() << "\n";
  cerr << "Number of triangles created from retesselation of tris from mesh 1: " << trianglesFromRetesselation[1].size() << "\n";




  /*
  for(int meshIdToProcess=0;meshIdToProcess<2;meshIdToProcess++) {
  	vector<pair<array<double,3>,array<double,3>> > edgesToStore;

  	for(const pair<int,int> &edge:newTriEdgesFromEachMap[meshIdToProcess]) {
  		const Point *v1 = (edge.first>=0)?&vertices[meshIdToProcess][edge.first]:&vertices[2][-edge.first -1];
  		const Point *v2 = (edge.second>=0)?&vertices[meshIdToProcess][edge.second]:&vertices[2][-edge.second -1];
  		edgesToStore.push_back(make_pair(toDoublePoint(*v1),toDoublePoint(*v2)));
  	}

  	
  	if(meshIdToProcess==0)
  		storeEdgesAsGts("out/retesselatedMesh0.gts",edgesToStore);
  	else 
  		storeEdgesAsGts("out/retesselatedMesh1.gts",edgesToStore);
  }
  */
}





//returns the number of intersections found
//the sets are filled with the triangles (from the corresponding mesh) that intersect
unsigned long long  computeIntersections(const Nested3DGridWrapper *uniformGrid, unordered_set<const Triangle *> trianglesThatIntersect[2]) {
  timespec t0,t1;
  clock_gettime(CLOCK_REALTIME, &t0);
   

  unsigned long long totalIntersections = 0;

  vector<pair<Triangle *,Triangle *> > vtPairsTrianglesToProcess;
  
  cerr << "Getting pairs of triangles from grid cells" << "\n";
  getPairsTrianglesInSameUnifGridCells(uniformGrid,vtPairsTrianglesToProcess);
  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time creating list of pairs of triangles to process (intersection): " << convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 
  //pairsTrianglesToProcess.reserve(1804900);
  

  //vector<pair<Triangle *,Triangle *> > vtPairsTrianglesToProcess(pairsTrianglesToProcess.begin(),pairsTrianglesToProcess.end());

  int numPairsToTest = vtPairsTrianglesToProcess.size();
  cerr << "Num pairs to test: " << numPairsToTest << endl;
  vector<bool> pairsIntersect(numPairsToTest,false); //true if the corresponding pair intersects..

  unsigned long long totalTests = 0;

  //TODO: try to reduce amount of computation here...

  vector< pair< array<VertCoord,3>,array<VertCoord,3> > > edges;
  vector< pair<Triangle *,Triangle *> >  intersectingTrianglesThatGeneratedEdges;
  #pragma omp parallel
  {
    VertCoord p0InterLine[3],p1InterLine[3];
    VertCoord tempRationals[100];

    vector< pair< array<VertCoord,3>,array<VertCoord,3> > > edgesTemp;
    vector< pair<Triangle *,Triangle *> >  intersectingTrianglesTemp;
    //edgesTemp.reserve(numPairsToTest/10);
    
    unsigned long long totalIntersectionsTemp = 0;
    unsigned long long totalTestsTemp = 0;

    #pragma omp for
    for(int i=0;i<numPairsToTest;i++) {   
      pair<Triangle *,Triangle *> &pairTriangles = vtPairsTrianglesToProcess[i];
      int coplanar;
      Triangle &a = *pairTriangles.first;
      Triangle &b = *pairTriangles.second;

      
      int ans = tri_tri_intersect_with_isectline(vertices[0][a.p[0]].data(),vertices[0][a.p[1]].data(),vertices[0][a.p[2]].data()     ,     
                                                  vertices[1][b.p[0]].data(),vertices[1][b.p[1]].data(),vertices[1][b.p[2]].data(),
                                                  &coplanar, p0InterLine ,p1InterLine,tempRationals);

      totalTestsTemp++;
      if (ans && !coplanar) {
        edgesTemp.push_back( pair<array<VertCoord,3>,array<VertCoord,3>>( {p0InterLine[0],p0InterLine[1],p0InterLine[2]} , {p1InterLine[0],p1InterLine[1],p1InterLine[2]}   ) );
      	
        #ifdef DEBUGGING_MODE
      	if(edgesTemp.back().first == edgesTemp.back().second) {
      		#pragma omp critical
      		{
      		cerr << "Triangle 1: " << endl;
      		printVertexForDebugging(vertices[0][a.p[0]].data());
      		printVertexForDebugging(vertices[0][a.p[1]].data());
      		printVertexForDebugging(vertices[0][a.p[2]].data());

      		cerr << "Triangle 2: " << endl;
      		printVertexForDebugging(vertices[1][a.p[0]].data());
      		printVertexForDebugging(vertices[1][a.p[1]].data());
      		printVertexForDebugging(vertices[1][a.p[2]].data());

      		cerr << "Inters:" << endl;
      		printVertexForDebugging(edgesTemp.back().first);
      		printVertexForDebugging(edgesTemp.back().second);
      		assert(false);
      		}
      	}
      	#endif
      	intersectingTrianglesTemp.push_back(pairTriangles);
      }

      if(ans) {
        totalIntersectionsTemp++;        
        pairsIntersect[i] = true;
      } 
    }

    #pragma omp critical
    {
      totalIntersections += totalIntersectionsTemp;
      totalTests += totalTestsTemp;
      edges.insert(edges.end(), edgesTemp.begin(), edgesTemp.end());
      intersectingTrianglesThatGeneratedEdges.insert(intersectingTrianglesThatGeneratedEdges.end(), intersectingTrianglesTemp.begin(), intersectingTrianglesTemp.end());
    }
  }
  
  clock_gettime(CLOCK_REALTIME, &t1);
  timeDetectIntersections = convertTimeMsecs(diff(t0,t1))/1000; 




  clock_gettime(CLOCK_REALTIME, &t0);
  retesselateIntersectingTriangles(edges,intersectingTrianglesThatGeneratedEdges);
  clock_gettime(CLOCK_REALTIME, &t1);
  timeRetesselate = convertTimeMsecs(diff(t0,t1))/1000; 



  //Some statistics...

  for(int i=0;i<numPairsToTest;i++)  
  	if(pairsIntersect[i]) {
  		trianglesThatIntersect[0].insert(vtPairsTrianglesToProcess[i].first);
  		trianglesThatIntersect[1].insert(vtPairsTrianglesToProcess[i].second);
  	}



  //TODO: remove this from timing data (this is just a statistic)
  map<const Triangle *,int> ctIntersectionsEachTriangleFromMap[2];
  for(int i=0;i<numPairsToTest;i++)  
  	if(pairsIntersect[i]) {
  		ctIntersectionsEachTriangleFromMap[0][vtPairsTrianglesToProcess[i].first]++;
  		ctIntersectionsEachTriangleFromMap[1][vtPairsTrianglesToProcess[i].second]++;
  	}  

  int ctTrianglesIntersect[2]  = {0,0};	
  int sumIntersections[2] = {0,0};	
  int maxIntersections[2] = {0,0};
  for(int meshId=0;meshId<2;meshId++) {
  	for(auto &a:ctIntersectionsEachTriangleFromMap[meshId]) {
  		if(a.second>maxIntersections[meshId])
  			maxIntersections[meshId] = a.second;
  		ctTrianglesIntersect[meshId]++;
  		sumIntersections[meshId]+= a.second;
  	}
  	cerr << "Mesh: " << meshId << endl;
  	cerr << "Max intersections                 : " << maxIntersections[meshId] << endl;
  	cerr << "Total intersections               : " << sumIntersections[meshId] << endl;
  	cerr << "Mum triangles intersecting        : " << ctTrianglesIntersect[meshId] << endl;
  	cerr << "Average intersections per triangle: " << sumIntersections[meshId]*1.0/ctTrianglesIntersect[meshId] << endl;
  }


  cerr << "Total tests: " << totalTests << endl;
  cerr << "Total intersections: " << totalIntersections << endl;


  //for debugging/visualization, let's store in a file for each triangle the original edges of that triangle and 
  //the edges formed by the intersection of other triangles...

  
  /*vector<pair<Triangle *,Triangle *> > pairsIntersectingTriangles;
  for(int i=0;i<numPairsToTest;i++)  
  	if(pairsIntersect[i]) {
  		pairsIntersectingTriangles.push_back(vtPairsTrianglesToProcess[i]);
  	}
  storeTriangleIntersections(pairsIntersectingTriangles); // TODO: use the edges vector (we already have the edges from intersections!!!)
  */

           
  return totalIntersections;
}





void readInputMesh(int meshId, string path) {
	assert(path.size()>=4);
  string inputFileExtension = path.substr(path.size()-3,3);
  assert(inputFileExtension=="gts" || inputFileExtension=="ium");
  bool isGtsFile = inputFileExtension=="gts";
  if(isGtsFile)
    TIME(readGTSFile(path,vertices[meshId],triangles[meshId],bBox[meshId]));
  else 
    TIME(readLiumFile(path,vertices[meshId],triangles[meshId],bBox[meshId]));
}




// TODO: if(a*b > 0 ) --> do not multiply!!! use test a>>0 && b>0 || a<0 && b<0

int main(int argc, char **argv) {
  if (argc!=7) {
      cerr << "Error... use ./3dIntersection inputMesh0 inputMesh1 gridSizeLevel1 gridSizeLevel2 triggerSecondLevel outputFile.gts" << endl;
      cerr << "The mesh file may be in the gts format or in the lium format (multimaterial)" << endl;
      exit(1);
  }
  clock_gettime(CLOCK_REALTIME, &t0BeginProgram);

  string mesh0Path = argv[1];
  string mesh1Path = argv[2];
  int gridSizeLevel1 = atoi(argv[3]);
  int gridSizeLevel2 = atoi(argv[4]);
  int triggerSecondLevel = atoi(argv[5]);
  int maxTreeDepth = 2;

  Print_Current_Process_Memory_Used();
  cerr << "Reading meshes..." << endl;

  timespec t0,t1;
  clock_gettime(CLOCK_REALTIME, &t0); 

  readInputMesh(0,mesh0Path);
  readInputMesh(1,mesh1Path);

  clock_gettime(CLOCK_REALTIME, &t1);
  t0AfterDatasetRead = t1;
  cerr << "Time to read: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  timeReadData = convertTimeMsecs(diff(t0,t1))/1000; 
  Print_Current_Process_Memory_Used();


  boundingBoxTwoMeshesTogetter[0] = bBox[0][0];
  boundingBoxTwoMeshesTogetter[1] = bBox[0][1];

  for(int i=0;i<3;i++) {
  	if(bBox[1][0][i] < boundingBoxTwoMeshesTogetter[0][i])
  		boundingBoxTwoMeshesTogetter[0][i] = bBox[1][0][i];
  	if(bBox[1][1][i] > boundingBoxTwoMeshesTogetter[1][i])
  		boundingBoxTwoMeshesTogetter[1][i] = bBox[1][1][i];
  }
  

  

  

  

  for(int meshId=0;meshId<2;meshId++)
  	cerr << "Bounding box mesh " << meshId << ": " << setprecision(18) << std::fixed << bBox[meshId][0][0].get_d() << " " << bBox[meshId][0][1].get_d() <<  " " << bBox[meshId][0][2].get_d() << " -- " << bBox[meshId][1][0].get_d() << " " << bBox[meshId][1][1].get_d() << " " << bBox[meshId][1][2].get_d() <<endl;
  cerr << "Bounding box two meshes togetter " << ": " << setprecision(18) << std::fixed << boundingBoxTwoMeshesTogetter[0][0].get_d() << " " << boundingBoxTwoMeshesTogetter[0][1].get_d() <<  " " << boundingBoxTwoMeshesTogetter[0][2].get_d() << " -- " << boundingBoxTwoMeshesTogetter[1][0].get_d() << " " << boundingBoxTwoMeshesTogetter[1][1].get_d() << " " << boundingBoxTwoMeshesTogetter[1][2].get_d() <<endl;


  
  cerr <<"Creating nested grid..." << endl;

  clock_gettime(CLOCK_REALTIME, &t0); 
  vector<Triangle *> trianglesPointers[2];

  int sz = triangles[0].size();
  trianglesPointers[0].resize(sz);
  for(int i=0;i<sz;i++) trianglesPointers[0][i]  = & (triangles[0][i]);

  sz = triangles[1].size();
  trianglesPointers[1].resize(sz);
  for(int i=0;i<sz;i++) trianglesPointers[1][i]  = & (triangles[1][i]);	


//  void init(const vector<Triangle *> trianglesInsert[2], const vector<Point> vertices[2], const int gridSizeLevel1, const int gridSizeLevel2, const Point &p0, const Point &p1,const long long prodThreshold);

  TIME(uniformGrid.init(trianglesPointers, vertices, gridSizeLevel1,gridSizeLevel2,boundingBoxTwoMeshesTogetter[0],boundingBoxTwoMeshesTogetter[1],triggerSecondLevel));


  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to create and refine grid: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  timeCreateGrid = convertTimeMsecs(diff(t0,t1))/1000; 

  //After the uniform grid is initialized, let's compute the intersection between the triangles...
  cerr << "Detecting intersections..." << endl;
  clock_gettime(CLOCK_REALTIME, &t0); 
  unordered_set<const Triangle *> trianglesThatIntersect[2];
  unsigned long long numIntersectionsDetected = computeIntersections(&uniformGrid,trianglesThatIntersect);

  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to detect intersections (includes time for computing statistics and for saving intersections for debugging): " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  Print_Current_Process_Memory_Used();



  clock_gettime(CLOCK_REALTIME, &t0); 

  ofstream outputStream(argv[6]);
  assert(outputStream);

  classifyTrianglesAndGenerateOutput(&uniformGrid, 
                                        trianglesThatIntersect,
                                        polygonsFromRetesselation,
                                        vertices,
                                        triangles,
                                        outputStream);

  clock_gettime(CLOCK_REALTIME, &t1);
  timeClassifyTriangles = convertTimeMsecs(diff(t0,t1))/1000;
  cerr << "Total time to classify triangles and generate output: " << timeClassifyTriangles << endl;
  Print_Current_Process_Memory_Used();


  cerr << "----------------------------------------------------" << endl;
  cerr << "Summary of ACTUAL times (excluding the times to compute statistics, to write edges for debugging, etc): " << endl;
  cerr << "Time to read the data         : " << timeReadData << endl;
  cerr << "Time to create and refine grid: " << timeCreateGrid << endl;
  cerr << "Time to detect intersections  : " << timeDetectIntersections << endl;  
  cerr << "Time to retesselate trinagles : " << timeRetesselate << endl;
  cerr << "Time to classify the triangles: " << timeClassifyTriangles << endl;
  cerr << "----------------------------------------------------" << endl;
/*
  if(isGtsFile)
    for(const ObjectId& id:pointIds) {
      cout << id << "\n";
    }
  else
    for(const ObjectId& id:pointIds) {  //when we deal with LIUM files the -1 represents the outside object and objects start with id 0.... (we convert from -1 to 0 and from 0 to 1 in the input... so, we need to revert this change...s)
        cout << id-1 << "\n";
    } 
*/
}