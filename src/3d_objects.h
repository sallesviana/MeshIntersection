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
#ifndef D3_OBJECTS_H
#define D3_OBJECTS_H

#include <string>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <array>
#include <vector>
#include <iostream>
#include "rationals.h"
#include "utils.h"


using namespace std;

//===============================================================
// Definitions...
//===============================================================

//typedef rational VertCoord;
//#define slightlyMoreThanOne VertCoord("10001/10000")



typedef rational VertCoord;
#define slightlyMoreThanOne 1.000001

typedef int ObjectId;
typedef int VertexId;
typedef double Area;
typedef array<int,3> CellNo;//Identifier for a cell

typedef array<VertCoord,3> Point;

const ObjectId OUTSIDE_OBJECT = 0;
const ObjectId DONT_KNOW_ID = -1;

//Ids of the planes..
const int PLANE_X0 =0;
const int PLANE_Y0 =1;
const int PLANE_Z0 =2;


//A triangle is oriented following the right hand rule
class Triangle {
public:
	VertexId p[3];
	ObjectId above, below; //ids of the objects above and below the triangle (considering the right hand rule)

	Triangle() {}

	Triangle(VertexId p0, VertexId p1, VertexId p2,ObjectId above, ObjectId below, vector<Point> &vertices) {
		//cerr << p0 << " " << p1 << " " << p2 << " " << above << " " << below << endl;
		p[0] = p0;
		p[1] = p1;
		p[2] = p2;
		this->above = above;
		this->below = below;

		for(int i=0;i<3;i++) {
			boundingBox[0][i] = p0;
			if (vertices[p1][i] < vertices[p0][i]) boundingBox[0][i] = 	p1;
			if (vertices[p2][i] < vertices[ boundingBox[0][i] ][i]) boundingBox[0][i] = 	p2;
		}
		for(int i=0;i<3;i++) {
			boundingBox[1][i] = p0;
			if (vertices[p1][i] > vertices[p0][i]) boundingBox[1][i] = 	p1;
			if (vertices[p2][i] > vertices[ boundingBox[1][i] ][i]) boundingBox[1][i] = 	p2;
		}		
	}


	//for efficiency purposes, we will store the bounding-box of the Triangles
	//In order to avoid storing the coordinates, each of the two points representing the bouding box (boundingBox[0] and [1]) will be represented by
	//the id of the vertex that has the corresponding coordinate as extreme.
	//For example, if boundingBox[0][2] = 6 --> the "smaller" vertex of the bouding box will have z coordinate similar to the z coordinate of the
	//vertex with id 6
	VertexId boundingBox[2][3];


	VertexId& operator[](const int i) { return this->p[i]; }
  const VertexId& operator[](const int i) const { return this->p[i]; }
};


//Triangle without bounding box...
class TriangleNoBB {
public:
	VertexId p[3];
	ObjectId above, below; //ids of the objects above and below the triangle (considering the right hand rule)

	TriangleNoBB() {}

	TriangleNoBB(VertexId p0, VertexId p1, VertexId p2,ObjectId above, ObjectId below) {
		p[0] = p0;
		p[1] = p1;
		p[2] = p2;
		this->above = above;
		this->below = below;	
	}

	VertexId& operator[](const int i) { return this->p[i]; }
  const VertexId& operator[](const int i) const { return this->p[i]; }
};


//used by the triangulation algorithm....
struct TriangulationVertex {
	int index; //index of the vertex in the "vertexSequence" vector
	int pNext,pPrev;
	int crNext,crPrev; //as implemented by [https://www.geometrictools.com/Documentation/TriangulationByEarClipping.pdf], the two lists are disjoint (no need to have separate pointers to both...)
	int eNext,ePrev;
	bool convex; //convex or reflex?
	bool ear; //is ear?

	TriangulationVertex(): pNext(-1),pPrev(-1),crNext(-1),crPrev(-1),eNext(-1),ePrev(-1),convex(false),ear(false) {}

	void print() { //for debugging purposes..
		cerr << index << " ; " << pNext << " " << pPrev << " ; " << crNext << " " << crPrev << " ; " << eNext << " " << ePrev << " ; " << convex << " " << ear << endl;
	}
};

class BoundaryPolygon {
public:
	vector<VertexId> vertexSequence; //the vertices in the vertexSequence always end with the first vertex... (e.g.: 1,5,9,4,1)
	vector<array<VertexId,3> > triangulatedPolygon; //stores a triangulated version of this polygon...
	ObjectId above, below;
	int whatPlaneProjectTriangleTo;

	BoundaryPolygon(const int whatPlaneProjectTrianglesTo_): whatPlaneProjectTrianglesTo(whatPlaneProjectTrianglesTo_) {}

	void triangulatePolygon(const vector<Point> vertices[3],const int meshIdWherePolygonIs, VertCoord tempCoords[]);

	void reverseVerticesOrder(); //reverse the order of the vertices in the vertex sequence (as a consequence, above and below are swapped)
	
	void printTriangles(const vector<Point> vertices) { //for debugging purposes..
		for(const array<VertexId,3>&t:triangulatedPolygon) {
			cout<< vertices[t[0]][0] << " " << vertices[t[0]][1] << "\n";
			cout<< vertices[t[1]][0] << " " << vertices[t[1]][1] << "\n";
			cout<< vertices[t[2]][0] << " " << vertices[t[2]][1] << "\n";
			cout<< vertices[t[0]][0] << " " << vertices[t[0]][1] << "\n}\n";			
		}
	}

	void printPolygon(const vector<Point> vertices) {
		for(VertexId v:vertexSequence) {
			cout << v << " | " << vertices[v][0] << " " << vertices[v][1] << " " << vertices[v][2] << "\n";
		}
	}

	bool isInClockwiseDirection(const vector<Point> vertices[3],const int meshIdWherePolygonIs) const; //for debugging purposes
private:
	int whatPlaneProjectTrianglesTo; //to what plane can we project this polygon without creating a degenerate polygon?

	bool pointInTriangleProj(const Point &p0, const Point &p1, const Point &p2, const Point &p);
	void initializeLinkedList(TriangulationVertex listVerticesToProcess[],int numVerticesPolygon,const vector<Point> vertices[3], 
																						const int meshIdWherePolygonIs, int &eBegin,int &eEnd,
																						int &cBegin,int &cEnd,
																						int &rBegin,int &rEnd, VertCoord tempCoords[]);
	void updateStatusVertex(int vertexId,TriangulationVertex listVerticesToProcess[], const vector<Point> vertices[3],const int meshIdWherePolygonIs,
													int &rBegin,int &rEnd, int &eBegin, VertCoord tempCoords[]);

	//vertexId is the position of the vertex in the "listVerticesToProcess" arrat...
	bool isEar(int vertexId,TriangulationVertex listVerticesToProcess[], const vector<Point> vertices[3],const int meshIdWherePolygonIs, const int rBegin);
	
	//TempCoords should have at least 2 coordinates
	bool isConvex(int vertexId,TriangulationVertex listVerticesToProcess[], const vector<Point> vertices[3],const int meshIdWherePolygonIs, VertCoord tempCoords[]);
};



//Reads a GTS file, fills the boundingBox with the boundingBox of the triangles read
void readGTSFile(string fileName, vector<Point> &vertices,vector<Triangle> &triangles, Point boundingBox[2]);


//Reads a Lium file, fills the boundingBox with the boundingBox of the triangles read
void readLiumFile(string fileName, vector<Point> &vertices,vector<Triangle> &triangles, Point boundingBox[2]);

#endif