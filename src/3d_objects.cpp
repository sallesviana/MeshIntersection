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

#include <string>
#include <fstream>
#include <cstdio>
#include <cassert>
#include "3d_objects.h"







const Point * getPointFromVertexId(int vertexId, int meshId,const vector<Point> vertices[3]) {
 if(vertexId>=0) assert(vertexId < vertices[meshId].size()); //TODO: remover assertions...
 else assert(-vertexId -1 < vertices[2].size());

 return (vertexId>=0)?&vertices[meshId][vertexId]:&vertices[2][-vertexId -1];
}

void BoundaryPolygon::reverseVerticesOrder() {
	std::reverse(vertexSequence.begin(),vertexSequence.end());
	swap(above,below);
}


bool BoundaryPolygon::isInClockwiseDirection(const vector<Point> vertices[3],const int meshIdWherePolygonIs) const { //for debugging purposes
	int numV = vertexSequence.size();
	//supposing the triangle will be projected to z=0...
	int coordX = 0;
  int coordY = 1;  
  if(whatPlaneProjectTrianglesTo==PLANE_X0) { //if the triangle is projected to X=0 --> we need to use coordinates y,z (instead of x,y)
  	coordX = 1;
  	coordY = 2;
  } else if(whatPlaneProjectTrianglesTo ==PLANE_Y0) { //if the triangle is projected to Y=0 --> we need to use coordinates z,x (instead of x,y)
  	coordX = 2;
  	coordY = 0;
  }


	double sum =0;
	for(int i=0;i<numV-1;i++) {
		const Point &u = *getPointFromVertexId(vertexSequence[i],meshIdWherePolygonIs,vertices);
		const Point &v = *getPointFromVertexId(vertexSequence[i+1],meshIdWherePolygonIs,vertices);
		double x2 = v[coordX].get_d();
		double x1 = u[coordX].get_d();
		double y2 = v[coordY].get_d();
		double y1 = u[coordY].get_d();
		//cerr << x1 << " " << y1 << endl;
		//cerr << x2 << " " << y2 << endl;
		sum += (x2-x1)*(y2+y1);
	}

	//cerr << "What plane: "  << whatPlaneProjectTrianglesTo << " "<< sum << endl;

	return sum>=-1e-10;
}

/*
Implementation based on: TODO: license...
// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2016
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
*/

//this function will triangulate the polygon, storing the triangles in the vector triangulatedPolygon
/*
Algorithm: https://www.geometrictools.com/Documentation/TriangulationByEarClipping.pdf


Store polygon as a double connected linked list
Store convex vertices in a list
Store reflex vertices in a list
Find initial ears and store their tips in a list
For each convex vertex, check if a reflex vertex is inside it (if none --> ear)
Repeat:
	Remove an ear creating a triangle
	Update the status of neighbor vertices
*/


/*
struct TriangulationVertex {
	int index; //index of the vertex in the "vertexSequence" vector
	int pNext,pPrev; //linked list defining the polygon to triangulate, circular LL
	int eNext,ePrev; //linked list defining the ear tips, circular LL
	int cNext,cPrev; //linked list defining  the convex vertices
	int rNext,rPrev; //linked list defining the reflex vertices	
	bool convex; //convex or reflex?
	bool ear; //is ear?
};
*/
//does not initialize the ears..

bool BoundaryPolygon::pointInTriangleProj(const Point &p0, const Point &p1, const Point &p2, const Point &p)  {
  if ( p==p0 || p==p1 || p==p2) return false; //is the point directly above a vertex of the triangle?

  //supposing the triangle will be projected to z=0...
  int coordY = 1;
  int coordX = 0;
  if(whatPlaneProjectTrianglesTo==PLANE_X0) { //if the triangle is projected to X=0 --> we need to use coordinates y,z (instead of x,y)
  	coordX = 1;
  	coordY = 2;
  } else if(whatPlaneProjectTrianglesTo ==PLANE_Y0) { //if the triangle is projected to Y=0 --> we need to use coordinates z,x (instead of x,y)
  	coordX = 2;
  	coordY = 0;
  }
  
  VertCoord denominator = ((p1[coordY] - p2[coordY])*(p0[coordX] - p2[coordX]) + (p2[coordX] - p1[coordX])*(p0[coordY] - p2[coordY]));
  if (denominator==0) { //TODO: check this.... degenerate triangles or vertical triangles (a segment never intersects a vertical triangle...)
    return false;
  }
  VertCoord a = ((p1[coordY] - p2[coordY])*(p[coordX] - p2[coordX]) + (p2[coordX] - p1[coordX])*(p[coordY] - p2[coordY])) / denominator;
  if ( a<0 || a >1) return false;
  
  VertCoord b = ((p2[coordY] - p0[coordY])*(p[coordX] - p2[coordX]) + (p0[coordX] - p2[coordX])*(p[coordY] - p2[coordY])) / denominator;

  if (b<0 || b>1) return false;
  VertCoord c = 1 - a - b;
 
  //if ( (fabs(a) <= EPS && (fabs(b) <= EPS) || fabs(c) <= EPS)  || (fabs(b) <= EPS && fabs(c) <= EPS) ) return false; // the point is one of the 3 triangle vertices...
  //if ( (fabs(a) <= EPS && fabs(b) <= EPS) || (fabs(a) <= EPS && fabs(c) <= EPS)  || (fabs(b) <= EPS && fabs(c) <= EPS) ) return false; // the point is one of the 3 triangle vertices...
  //return 0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1; //maybe we should perform some tests using EPSILON to avoid floating point errors...
  
  
  return 0 <= c && c <= 1; //maybe we should perform some tests using EPSILON to avoid floating point errors...
}



bool BoundaryPolygon::isEar(int vertexId,TriangulationVertex listVerticesToProcess[], const vector<Point> vertices[3],const int meshIdWherePolygonIs,const int rBegin) {
	//check if vertex is convex and no reflex vertex is inside triangle formed with neighbors...	
	if(rBegin!=-1) {
		//if there is at least a reflex vertex:	
		listVerticesToProcess[vertexId].ear = true;
		//check if there is a reflex vertex inside the triangle defined by the neighbor vertices...
		const int v1 = listVerticesToProcess[listVerticesToProcess[vertexId].pPrev].index;
		const int v2 = listVerticesToProcess[vertexId].index;
		const int v3 = listVerticesToProcess[listVerticesToProcess[vertexId].pNext].index;

		const Point &vertex1 = *getPointFromVertexId(v1,meshIdWherePolygonIs,vertices);//vertices[meshIdWherePolygonIs][v1]; //TODO: shared vertices, etc...
		const Point &vertex2 = *getPointFromVertexId(v2,meshIdWherePolygonIs,vertices);//vertices[meshIdWherePolygonIs][v2];
		const Point &vertex3 = *getPointFromVertexId(v3,meshIdWherePolygonIs,vertices);//vertices[meshIdWherePolygonIs][v3];

		int reflexId = rBegin;
		while(reflexId!=-1) {
			const Point &reflex = *getPointFromVertexId(reflexId,meshIdWherePolygonIs,vertices);//vertices[meshIdWherePolygonIs][reflexId]; //TODO: consider shared vertices, etc..
			//if reflex is inside the triangle --> answer is false!
			if(reflexId == listVerticesToProcess[vertexId].pPrev || reflexId==vertexId || reflexId==listVerticesToProcess[vertexId].pNext) {
				reflexId = listVerticesToProcess[reflexId].crNext; //next reflex..
				continue;
			}

			if(pointInTriangleProj(vertex1,vertex2,vertex3,reflex)) {
				listVerticesToProcess[vertexId].ear = false; 
				break;
			}
			reflexId = listVerticesToProcess[reflexId].crNext; //next reflex..
		}

	} else {
		listVerticesToProcess[vertexId].ear = true; //no reflex vertex --> can only be ear (we obviously do not need to test if it is convex since there is no reflex vertex! --> all are convex)
	}
	return listVerticesToProcess[vertexId].ear;
}

//TempCoords should have at least 2 coordinates
bool BoundaryPolygon::isConvex(int vertexId,TriangulationVertex listVerticesToProcess[], const vector<Point> vertices[3],const int meshIdWherePolygonIs, VertCoord tempCoords[]) {
	const int v1 = listVerticesToProcess[listVerticesToProcess[vertexId].pPrev].index;
	const int v2 = listVerticesToProcess[vertexId].index;
	const int v3 = listVerticesToProcess[listVerticesToProcess[vertexId].pNext].index;

	const Point &vertex1 = *getPointFromVertexId(v1,meshIdWherePolygonIs,vertices);//vertices[meshIdWherePolygonIs][v1]; //TODO: shared vertices, etc...
	const Point &vertex2 = *getPointFromVertexId(v2,meshIdWherePolygonIs,vertices);//vertices[meshIdWherePolygonIs][v2];
	const Point &vertex3 = *getPointFromVertexId(v3,meshIdWherePolygonIs,vertices);//vertices[meshIdWherePolygonIs][v3];

	//supposing the triangle will be projected to z=0...
  int coordY = 1;
  int coordX = 0;
  if(whatPlaneProjectTrianglesTo==PLANE_X0) { //if the triangle is projected to X=0 --> we need to use coordinates y,z (instead of x,y)
  	coordX = 1;
  	coordY = 2;
  } else if(whatPlaneProjectTrianglesTo ==PLANE_Y0) { //if the triangle is projected to Y=0 --> we need to use coordinates z,x (instead of x,y)
  	coordX = 2;
  	coordY = 0;
  }

  /*
 	VertCoord twoArea = vertex1[coordX]*(vertex2[coordY]-vertex3[coordY]) +
 									    vertex2[coordX]*(vertex3[coordY]-vertex1[coordY]) +
 									    vertex3[coordX]*(vertex1[coordY]-vertex2[coordY]) ;

 	listVerticesToProcess[vertexId].convex = twoArea<0;
	return listVerticesToProcess[vertexId].convex;
	*/
	tempCoords[0] = vertex2[coordY]; //VertCoord twoArea = vertex1[coordX]*(vertex2[coordY]-vertex3[coordY]);
	tempCoords[0] -= vertex3[coordY];
	tempCoords[0] *= vertex1[coordX];

	tempCoords[1] = vertex3[coordY];
	tempCoords[1] -= vertex1[coordY];
	tempCoords[1] *= vertex2[coordX]; 
	tempCoords[0] += tempCoords[1]; //twoArea += vertex2[coordX]*(vertex3[coordY]-vertex1[coordY]);
 	
 	tempCoords[1] = vertex1[coordY];
	tempCoords[1] -= vertex2[coordY];
	tempCoords[1] *= vertex3[coordX]; 
	tempCoords[0] += tempCoords[1]; //twoArea += vertex3[coordX]*(vertex1[coordY]-vertex2[coordY]);
 	

 	listVerticesToProcess[vertexId].convex = sgn(tempCoords[0])<0;
	return listVerticesToProcess[vertexId].convex;
}

void BoundaryPolygon::initializeLinkedList(TriangulationVertex listVerticesToProcess[],int numVerticesPolygon,const vector<Point> vertices[3], 
																						const int meshIdWherePolygonIs, int &eBegin,int &eEnd,
																						int &cBegin,int &cEnd,
																						int &rBegin,int &rEnd, VertCoord tempCoords[]) {
	//if all vertices are convex --> no need to look for ears (all vertices are ear...)
	eBegin = eEnd = cBegin = cEnd = rBegin = rEnd = -1;
	for(int i = 0;i<numVerticesPolygon;i++) {
		TriangulationVertex &v = listVerticesToProcess[i];
		v.index = vertexSequence[i]; //maybe (probably...) unnecessary...
		v.pNext = (i+1);
		if(v.pNext>=numVerticesPolygon) v.pNext = 0;
		v.pPrev = (i-1);
		if(v.pPrev<0) v.pPrev = numVerticesPolygon-1;
	}

	eBegin = eEnd = -1;

	//now,let's update the list fo convex and reflex edges...
	cBegin = cEnd = rBegin = rEnd = -1;
	for(int i=0;i<numVerticesPolygon;i++) {
		if(isConvex(i,listVerticesToProcess, vertices, meshIdWherePolygonIs,tempCoords)) {
			if(cBegin==-1) {
				cBegin = cEnd = i;
				listVerticesToProcess[i].crNext = -1;
				listVerticesToProcess[i].crPrev = -1;
			} else {
				listVerticesToProcess[cEnd].crNext = i;
				listVerticesToProcess[i].crPrev = cEnd;
				cEnd = i;
			}
		} else {
			if(rBegin==-1) {
				rBegin = rEnd = i;
				listVerticesToProcess[i].crNext = -1;
				listVerticesToProcess[i].crPrev = -1;
			} else {
				listVerticesToProcess[rEnd].crNext = i;
				listVerticesToProcess[i].crPrev = rEnd;
				rEnd = i;
			}
		}
	}
}


//updates the status of vertices after one ear is removed
//the vertexId must be a neighbor of the removed ear...
void BoundaryPolygon::updateStatusVertex(int vertexId,TriangulationVertex listVerticesToProcess[], const vector<Point> vertices[3],
																					const int meshIdWherePolygonIs,int &rBegin,int &rEnd, int &eBegin, VertCoord tempCoords[]) {
	int wasEar = listVerticesToProcess[vertexId].ear;
	if(wasEar) { //if the vertex was an ear... it may not be an ear now...
		if(!isEar(vertexId,listVerticesToProcess, vertices, meshIdWherePolygonIs,rBegin)) {
			//remove the neighbor from the list of ears...
			int eNext = listVerticesToProcess[vertexId].eNext;
			int ePrev = listVerticesToProcess[vertexId].ePrev;
			listVerticesToProcess[eNext].ePrev = ePrev;
			listVerticesToProcess[ePrev].eNext = eNext; 
		}
	} else { //maybe it is now an ear...
		//if it was a reflex --> it may still be reflex (don't change), converx or an ear...
		//if it was convex --> still convex...
		//So, what matters now is: is it convex now? if it is, it may be a new ear 
		int wasConvex = listVerticesToProcess[vertexId].convex;
		if(isConvex(vertexId,listVerticesToProcess, vertices, meshIdWherePolygonIs,tempCoords)) {
			if(!wasConvex) { //the reflex vertex became convex --> remove from reflex list (this accelerates the ear tests...)
				//we need to keep track of reflex vertices list --> we need to update the rBegin, rEnd if necessary...
				if(vertexId ==rBegin) {
					rBegin = listVerticesToProcess[vertexId].crNext;
					if(rBegin!=-1) { 
						listVerticesToProcess[rBegin].crPrev = -1;
					}
					listVerticesToProcess[vertexId].crNext = -1;
				}
				else if(vertexId == rEnd) {
					rEnd = listVerticesToProcess[vertexId].crPrev;
					if(rEnd!=-1) { 
						listVerticesToProcess[rBegin].crNext = -1;
					}
					listVerticesToProcess[vertexId].crPrev = -1;

				} else { //in the middle of the list --> easy
					int rNext = listVerticesToProcess[vertexId].crNext;
					int rPrev = listVerticesToProcess[vertexId].crPrev;
					listVerticesToProcess[rNext].crPrev = rPrev;
					listVerticesToProcess[rPrev].crNext = rNext; 

					listVerticesToProcess[vertexId].crPrev = -1;
					listVerticesToProcess[vertexId].crNext = -1;
				}
			} 
			if(isEar(vertexId,listVerticesToProcess, vertices, meshIdWherePolygonIs,rBegin)) { //a reflex became an ear...
				//insert it before the beginning of the circular list of ears (i.e., at the end)...
				int firstEarPrev = listVerticesToProcess[eBegin].ePrev;
				listVerticesToProcess[vertexId].ePrev = firstEarPrev;
				listVerticesToProcess[vertexId].eNext = eBegin;
				listVerticesToProcess[eBegin].ePrev = vertexId;
				listVerticesToProcess[firstEarPrev].eNext = vertexId;
			}

		} //if is reflex --> it was reflex before (a convex cannot become reflex --> the else does not matter...)

	}
}

void BoundaryPolygon::triangulatePolygon(const vector<Point> vertices[3],const int meshIdWherePolygonIs, VertCoord tempCoords[]) {
	int numVerticesPolygon = vertexSequence.size()-1; //the vertices always end with the first vertex...
	TriangulationVertex listVerticesToProcess[numVerticesPolygon]; //we will store here the linked lists

	int eBegin,eEnd;
	int cBegin,cEnd;
	int rBegin,rEnd;
	//initialize the linked list
	initializeLinkedList(listVerticesToProcess,numVerticesPolygon,vertices,meshIdWherePolygonIs,eBegin,eEnd,cBegin,cEnd,rBegin,rEnd,tempCoords);

	/*cerr << "Triangulating: " << endl;
	for(int i=0;i<vertexSequence.size();i++) {
		const Point &vertex1 = *getPointFromVertexId(vertexSequence[i],meshIdWherePolygonIs,vertices);
		cerr << i << " " << vertexSequence[i] << " " << vertex1[0].get_d() << " " << vertex1[1].get_d() << " " << vertex1[2].get_d() << endl;
		cerr << i << " " << vertexSequence[i] << " " << vertex1[0] << " " << vertex1[1] << " " << vertex1[2] << endl;
	}*/


	//first, check if the polygon is convex... TODO
	/*if (rBegin==-1) {
		cerr << "TODO: faster processing convex polygons..." << endl;
	}*/
	/*if(rBegin==-1) {
		//the polygon is convex...

		return;
	}*/
	/*cerr <<"What plane project to: " << whatPlaneProjectTriangleTo << endl;
	cerr << "cBegin, cEnd : " << cBegin << " " << cEnd << endl;
	cerr << "rBegin, rEnd : " << rBegin << " " << rEnd << endl;
	for(int i=0;i<numVerticesPolygon;i++) {
		listVerticesToProcess[i].print();
	}*/

	//if the polygon is not convex, let's initialize the list of ears
	int vertexPtr = cBegin;
	assert(vertexPtr!=-1); //we should have at least one convex vertex..
	while(vertexPtr!=-1) {
		//cerr << "Processing: " << vertexPtr << endl;
		if(isEar(vertexPtr,listVerticesToProcess, vertices, meshIdWherePolygonIs,rBegin)) {
			//cerr << "Is ear!" << endl;
			//let's add it to the list of ears...
			if(eBegin==-1) {
				eBegin = eEnd = vertexPtr;
			} else {
				listVerticesToProcess[eEnd].eNext = vertexPtr;
				listVerticesToProcess[vertexPtr].ePrev = eEnd;

				eEnd = vertexPtr;
			}
		}
		vertexPtr = listVerticesToProcess[vertexPtr].crNext; //get the next convex vertex...
	}
	assert(eBegin !=-1); //we should have at least one ear...

	listVerticesToProcess[eEnd].eNext = eBegin; //let's make this list circular...
	listVerticesToProcess[eBegin].ePrev = eEnd;

	/*
	for(int i=0;i<numVerticesPolygon;i++) {
		listVerticesToProcess[i].print();
	}*/

	//triangulatedPolygon.reserve... TODO

	//them, we remove each year at a time (always updating the neighbor vertices...)
	while(true) {
		//create a triangle with the first ear in the list of ears..
		int prevVertexEar = listVerticesToProcess[eBegin].pPrev;
		int nextVertexEar = listVerticesToProcess[eBegin].pNext;

		//cerr << "Removed an ear: " << prevVertexEar << " " << eBegin << " " << nextVertexEar << endl;
		triangulatedPolygon.push_back({listVerticesToProcess[prevVertexEar].index,listVerticesToProcess[eBegin].index,listVerticesToProcess[nextVertexEar].index});


		//when we remove an ear (a necessarely convex vertex) from the list of ears we:
		//- do not need to remove it from the list of convex (this does not matter b/c we do not use this list directly.. we only need list of reflex vertices to test for earness)
		//- need to remove it from the polygon.... (this matters because when an ear is removed we check if status of adjacent vertices (in polygon) changed...
		//we also do not need to update pBegin, pEnd (this will not be needed anymore...)
		{
			int pNext = listVerticesToProcess[eBegin].pNext;
			int pPrev = listVerticesToProcess[eBegin].pPrev;
			listVerticesToProcess[pNext].pPrev = pPrev;
			listVerticesToProcess[pPrev].pNext = pNext; 
		}

		numVerticesPolygon--; //we removed one of the vertices...
		if(numVerticesPolygon<=3) {
			if(numVerticesPolygon==3) { //if we have only 3 vertices --> extract this last triangle (this saves computation...)
				//remove the ear, and set the first ear as the next one...
				{
					int eNext = listVerticesToProcess[eBegin].eNext;
					int ePrev = listVerticesToProcess[eBegin].ePrev;
					listVerticesToProcess[eNext].ePrev = ePrev;
					listVerticesToProcess[ePrev].eNext = eNext; 
					eBegin = eNext; //makes the first new ear the next one...
				}
				prevVertexEar = listVerticesToProcess[eBegin].pPrev;
				nextVertexEar = listVerticesToProcess[eBegin].pNext;


				triangulatedPolygon.push_back({listVerticesToProcess[prevVertexEar].index,listVerticesToProcess[eBegin].index,listVerticesToProcess[nextVertexEar].index});
				break;
			} 
			else break; //we started with 3 vertices (a triangle) and removed the only ear
		}

		//update status of prevVertexEar and nextVertexEar
		updateStatusVertex(prevVertexEar,listVerticesToProcess,vertices,meshIdWherePolygonIs,rBegin,rEnd, eBegin,tempCoords);
		updateStatusVertex(nextVertexEar,listVerticesToProcess,vertices,meshIdWherePolygonIs,rBegin,rEnd, eBegin,tempCoords);

		//remove the ear, and set the first ear as the next one...
		{
			int eNext = listVerticesToProcess[eBegin].eNext;
			int ePrev = listVerticesToProcess[eBegin].ePrev;
			listVerticesToProcess[eNext].ePrev = ePrev;
			listVerticesToProcess[ePrev].eNext = eNext; 
			eBegin = eNext; //makes the first new ear the next one...
		}
	}

}
