#include "3dGeometry.h"




void BoundaryPolygon::reverseVerticesOrder() {
	std::reverse(vertexSequence.begin(),vertexSequence.end());
	swap(above,below);
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

bool BoundaryPolygon::pointInTriangleProj(MeshIntersectionGeometry &geometry,int p0,int p1, int p2, int queryPoint,
																						TempVarsTriangulatePolygon &tempVars)  {

  return geometry.isVertexInTriangleProjection(*vertexSequence[p0],*vertexSequence[p1], *vertexSequence[p2], *vertexSequence[queryPoint],whatPlaneProjectTrianglesTo,tempVars.tempVarsIsVertexTriangleProjection);
}


//TempCoords should have at least 2 coordinates
bool BoundaryPolygon::isConvex(MeshIntersectionGeometry &geometry,int vertexId,	TriangulationVertex listVerticesToProcess[], 
																TempVarsTriangulatePolygon &tempVars) {
	const Vertex* v1 = listVerticesToProcess[listVerticesToProcess[vertexId].pPrev].index;
	const Vertex* v2 = listVerticesToProcess[vertexId].index;
	const Vertex* v3 = listVerticesToProcess[listVerticesToProcess[vertexId].pNext].index;

 	listVerticesToProcess[vertexId].convex = geometry.isVertexConvex(*v1,*v2, *v3, whatPlaneProjectTrianglesTo,tempVars.tempVarsIsVertexConvex);
	return listVerticesToProcess[vertexId].convex;
}




bool BoundaryPolygon::isEar(MeshIntersectionGeometry &geometry,int vertexId,TriangulationVertex listVerticesToProcess[], 
															const int rBegin,TempVarsTriangulatePolygon &tempVars) {
	//check if vertex is convex and no reflex vertex is inside triangle formed with neighbors...	
	if(rBegin!=-1) {
		//if there is at least a reflex vertex:	
		listVerticesToProcess[vertexId].ear = true;
		//check if there is a reflex vertex inside the triangle defined by the neighbor vertices...
		/*const int v1 = listVerticesToProcess[listVerticesToProcess[vertexId].pPrev].index;
		const int v2 = listVerticesToProcess[vertexId].index;
		const int v3 = listVerticesToProcess[listVerticesToProcess[vertexId].pNext].index;*/
	

		int reflexId = rBegin;
		while(reflexId!=-1) {
			//if reflex is inside the triangle --> answer is false!
			if(reflexId == listVerticesToProcess[vertexId].pPrev || reflexId==vertexId || reflexId==listVerticesToProcess[vertexId].pNext) {
				reflexId = listVerticesToProcess[reflexId].crNext; //next reflex..
				continue;
			}

			if(pointInTriangleProj(geometry,listVerticesToProcess[vertexId].pPrev,vertexId,listVerticesToProcess[vertexId].pNext,reflexId,tempVars)) {
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


void BoundaryPolygon::initializeLinkedList(MeshIntersectionGeometry &geometry,TriangulationVertex listVerticesToProcess[],int numVerticesPolygon, 
																						int &eBegin,int &eEnd,
																						int &cBegin,int &cEnd,
																						int &rBegin,int &rEnd,TempVarsTriangulatePolygon &tempVars) {
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
		if(isConvex(geometry,i,listVerticesToProcess, tempVars)) {
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
void BoundaryPolygon::updateStatusVertex(MeshIntersectionGeometry &geometry,int vertexId,TriangulationVertex listVerticesToProcess[], 
																					int &rBegin,int &rEnd, int &eBegin,TempVarsTriangulatePolygon &tempVars) {
	int wasEar = listVerticesToProcess[vertexId].ear;
	if(wasEar) { //if the vertex was an ear... it may not be an ear now...
		if(!isEar(geometry,vertexId,listVerticesToProcess,rBegin,tempVars)) {
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
		if(isConvex(geometry,vertexId,listVerticesToProcess, tempVars)) {
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
						listVerticesToProcess[rEnd].crNext = -1;
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
			if(isEar(geometry,vertexId,listVerticesToProcess,  rBegin,tempVars)) { //a reflex became an ear...
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

void BoundaryPolygon::triangulatePolygon(MeshIntersectionGeometry &geometry,TempVarsTriangulatePolygon &tempVars) {
	int numVerticesPolygon = vertexSequence.size()-1; //the vertices always end with the first vertex...
	TriangulationVertex listVerticesToProcess[numVerticesPolygon]; //we will store here the linked lists

	int eBegin,eEnd;
	int cBegin,cEnd;
	int rBegin,rEnd;
	//initialize the linked list
	initializeLinkedList(geometry, listVerticesToProcess,numVerticesPolygon,eBegin,eEnd,cBegin,cEnd,rBegin,rEnd,tempVars);


	//if the polygon is not convex, let's initialize the list of ears
	int vertexPtr = cBegin;
	assert(vertexPtr!=-1); //we should have at least one convex vertex..
	while(vertexPtr!=-1) {
		//cerr << "Processing: " << vertexPtr << endl;
		if(isEar(geometry, vertexPtr,listVerticesToProcess,  rBegin,tempVars)) {
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
		updateStatusVertex(geometry,prevVertexEar,listVerticesToProcess,rBegin,rEnd, eBegin,tempVars);
		updateStatusVertex(geometry,nextVertexEar,listVerticesToProcess,rBegin,rEnd, eBegin,tempVars);

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
