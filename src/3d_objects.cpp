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

		/*cerr << "Triangle: " << endl;
		cerr << p0 << "  " << p1 << "  " << p2 << endl;
		for(int i=0;i<2;i++) {
			for(int j=0;j<3;j++)
				 cerr << "BBox: " << i << " " << j << "  " << boundingBox[i][j] << endl;
			cerr << endl;
		}*/
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







//Reads a GTS file, fills the boundingBox with the boundingBox of the triangles read
void readGTSFile(string fileName, vector<Point> &vertices,vector<Triangle> &triangles, Point boundingBox[2]) {
	FILE *inputFile = fopen(fileName.c_str(),"r");
	if (inputFile==NULL ) {
    cerr << "ERROR: failed to open file " << fileName << endl;
    exit(1);
  }
  cerr << "Reading file " << fileName << endl;
  int numVertices, numEdges, numTriangles;
	assert(fscanf(inputFile,"%d %d %d",&numVertices,&numEdges,&numTriangles)==3);
	vertices.resize(numVertices);
	triangles.resize(numTriangles);
	vector< array<VertexId,2> > edges(numEdges);

	//rational xr,yr,zr;
	for(int i=0;i<numVertices;i++) {
		double x,y,z;
		assert(fscanf(inputFile,"%lf %lf %lf",&x,&y,&z)==3);

		vertices[i][0] = x;
		vertices[i][1] = y;
		vertices[i][2] = z;

		if (i==0) { //initialize the bounding box in the first iteration...
			boundingBox[0] = vertices[0];
			boundingBox[1] = vertices[0];
		}

		accum_min(boundingBox[0][0],vertices[i][0]);
		accum_min(boundingBox[0][1],vertices[i][1]);
		accum_min(boundingBox[0][2],vertices[i][2]);

		accum_max(boundingBox[1][0],vertices[i][0]);
		accum_max(boundingBox[1][1],vertices[i][1]);
		accum_max(boundingBox[1][2],vertices[i][2]);
	}
	for(int i=0;i<numEdges;i++) {
		VertexId a,b;
		assert(fscanf(inputFile,"%d %d",&a,&b)==2);
		edges[i][0] = a-1; //GTS starts counting from 1... we will start the ids from 0
		edges[i][1] = b-1;
	}

	for(int i=0;i<numTriangles;i++) {
		int a,b,c;
		assert(fscanf(inputFile,"%d %d %d",&a,&b,&c)==3);
		a--;b--;c--; //we count from 0, not from 1...

		if (edges[a][0] == edges[b][0] || edges[a][0] == edges[b][1]) swap(edges[a][0], edges[a][1]); // we will ensure that the edges are in the format (a,b)-(b,c)-(d,e) or (a,b)-(c,b)-(d,e)
		if (edges[a][1] == edges[b][1]) swap(edges[b][0],edges[b][1]); //ensure that edges are in the format (a,b)-(b,c)-(d,e)
		if (edges[b][1] == edges[c][1]) swap(edges[c][0],edges[c][1]); //ensure that edges are in the format (a,b)-(b,c)-(c,a)

		assert(edges[a][1]==edges[b][0]);
		assert(edges[b][1]==edges[c][0]);
		assert(edges[c][1]==edges[a][0]);

		triangles[i] = Triangle(edges[a][0], edges[a][1],edges[b][1],OUTSIDE_OBJECT,1,vertices);
	}
}



//Reads a Lium file, fills the boundingBox with the boundingBox of the triangles read
void readLiumFile(string fileName, vector<Point> &vertices,vector<Triangle> &triangles, Point boundingBox[2]) {
	FILE *inputFile = fopen(fileName.c_str(),"r");
	if (inputFile==NULL ) {
    cerr << "ERROR: failed to open file " << fileName << endl;
    exit(1);
  }
  cerr << "Reading file " << fileName << endl;
  int numVertices, numTriangles,numObjects;
	assert(fscanf(inputFile,"%d %d %d",&numVertices,&numTriangles,&numObjects)==3);
	cerr << "Vertices, triangles, objects: " << numVertices << " " << numTriangles << " " << numObjects << endl;
 	vertices.resize(numVertices);
	triangles.resize(numTriangles);

	//rational xr,yr,zr;
	for(int i=0;i<numVertices;i++) {
		double x,y,z;
		assert(fscanf(inputFile,"%lf %lf %lf",&x,&y,&z)==3);

		vertices[i][0] = x;
		vertices[i][1] = y;
		vertices[i][2] = z;

		if (i==0) { //initialize the bounding box in the first iteration...
			boundingBox[0] = vertices[0];
			boundingBox[1] = vertices[0];
		}

		accum_min(boundingBox[0][0],vertices[i][0]);
		accum_min(boundingBox[0][1],vertices[i][1]);
		accum_min(boundingBox[0][2],vertices[i][2]);

		accum_max(boundingBox[1][0],vertices[i][0]);
		accum_max(boundingBox[1][1],vertices[i][1]);
		accum_max(boundingBox[1][2],vertices[i][2]);
	}
	
	for(int i=0;i<numTriangles;i++) {
		int a,b,c;
		int objContrNormal,objNormal;
		assert(fscanf(inputFile,"%d %d %d %d %d",&a,&b,&c,&objNormal,&objContrNormal)==5);
		//a--;b--;c--; //we count from 0, not from 1...
		
		objNormal++;
		objContrNormal++;
		//cerr << objNormal << " " << objContrNormal << endl;
		triangles[i] = Triangle(a,b,c,objNormal,objContrNormal,vertices);
	}
}
