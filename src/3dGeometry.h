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

//
class MeshIntersectionGeometry;

class Vertex {
	public:
		const int getId() const { //different vertices from the intersection (mesh 3) may have the same id
			return id;
		}
		const int getMeshId() const {
			return meshId;
		}
		Vertex(int meshId_,int id_): meshId(meshId_),id(id_){}
		Vertex() : meshId(-1), id(-1) {}

		friend class MeshIntersectionGeometry;
	private: 
		int id;
		int meshId;	
};

class InputVertex: public Vertex {
	public:		
		InputVertex(int meshId_,int id_): Vertex(meshId_,id_) {}
		InputVertex(): Vertex(-1,-1){}
	private:
			
};


//A triangle is oriented following the right hand rule
class Triangle {
public:	
	ObjectId above, below; //ids of the objects above and below the triangle (considering the right hand rule)

	Triangle() {}

	Triangle(ObjectId above, ObjectId below) {		
		this->above = above;
		this->below = below;				
	}

	const Vertex* getVertex(const int i) const;
};

class InputTriangle: public Triangle {
public:
	InputTriangle(const InputVertex &p0, const InputVertex &p1, const InputVertex &p2,ObjectId above, ObjectId below)
			:Triangle(above,below) {
		//cerr << p0 << " " << p1 << " " << p2 << " " << above << " " << below << endl;
		p[0] = p0;
		p[1] = p1;
		p[2] = p2;		
	}
	InputTriangle() {}

  const InputVertex* getInputVertex(const int i) const { 
  	return &p[i]; 
  }
  const Vertex* getVertex(const int i) const { return &p[i]; }
private:
	InputVertex p[3];
};


//Two vertices from the intersection may have the same ids (the ids are used only to help finding the positions of the
//vertex coordinates in arrays..)
class VertexFromIntersection: public Vertex {
	public:
		VertexFromIntersection(const InputVertex &edgeV1, const InputVertex &edgeV2, const InputTriangle triangleGeneratedVertex) {
			edge[0] = edgeV1;
			edge[1] = edgeV2;
			if(edge[0].getId()>edge[1].getId()) swap(edge[0],edge[1]);
			triangle = triangleGeneratedVertex;
		}
		VertexFromIntersection() {};
	//private:

		//this vertex is the intersection between the edge and the input triangle triangle 
		InputVertex edge[2];
		InputTriangle triangle;
};



//The same pool (even in different functions) cannot be used simultaneously by different threads...
class MeshIntersectionGeometry {
	public:
		MeshIntersectionGeometry(const string &pathMesh0, const string &pathMesh1);

		struct TempVarsGetGridCellContainingVertex { VertCoord tempVertCoords; big_int tempVarsInt[3];};
		array<int,3> getGridCellContainingVertex(const int meshId, const int vertexId, const VertCoord &cellScale, TempVarsGetGridCellContainingVertex &tempVars ) const;
		int getGridCellXContainingVertex(int meshId,const VertCoord &x, const VertCoord &cellScale, TempVarsGetGridCellContainingVertex &tempVars  ) const;
		int getGridCellYContainingVertex(int meshId,const VertCoord &y, const VertCoord &cellScale, TempVarsGetGridCellContainingVertex &tempVars  ) const;
		int getGridCellZContainingVertex(int meshId,const VertCoord &z, const VertCoord &cellScale, TempVarsGetGridCellContainingVertex &tempVars  ) const;


		vector<InputTriangle> inputTriangles[2];
		int getVertexIdInputTriangleWithMinCoord(int meshId, int triangleId, int coord) const { return inputTrianglesBoundingBox[meshId][triangleId][0][coord];}
		int getVertexIdInputTriangleWithMaxCoord(int meshId, int triangleId, int coord) const { return inputTrianglesBoundingBox[meshId][triangleId][1][coord];}
		int getNumVertices(int meshId) const {return verticesCoordinates[meshId].size();}
	

		//Given pairs of triangles to test for intersection, returns (fills) the edges from intersection and the pairs of triangles that intersect
		//verticesCoordinates[2] will be filled with the coordinates of vertices generated by intersections
		void computeIntersections(const vector<pair<InputTriangle *,InputTriangle *> > &inputTrianglesToConsider, vector< pair<InputTriangle *,InputTriangle *> >  &intersectingTrianglesThatGeneratedEdges, vector< pair<VertexFromIntersection, VertexFromIntersection> >  &edgesFromIntersection, unsigned long long &numIntersectionTests);	

		~MeshIntersectionGeometry() {cerr << "Delegint memory..." << endl;}


		void printBoundingBoxes();
		array<VertCoord,3> coordRangeMeshes() const; //"width" of the coordinates of the two meshes togetter (i.e., size of the common bounding-box in each coordinate)
	private:
		struct PlaneEquation {Point normal; VertCoord d;};
		vector<PlaneEquation> planeEquationsInputTriangles[2];
		vector<bool> isPlaneEquationInputTrianglesInitialized[2];

		void storeIntersectionVerticesCoordinatesAndUpdateVerticesIds(vector< pair<VertexFromIntersection, VertexFromIntersection> >  &edgesFromIntersection,const vector< pair<Point, Point> > &coordsVerticesOfEdges);
		
		struct TempVarsComputePlaneEquation {Point E1,E2; VertCoord temp;};
		void computePlaneEquation(PlaneEquation &equation, const Point &V0, const Point &V1, const Point &V2, TempVarsComputePlaneEquation &tempVars);
		const PlaneEquation &getPlaneEquationInputTriangle(int meshId, int triId,TempVarsComputePlaneEquation &tempVars);
		void initPlaneEquationInputTriangle(int meshId, int triId,TempVarsComputePlaneEquation &tempVars);

		vector<Point> verticesCoordinates[3]; //verticesCoordinates[0] are from mesh 0, verticesCoordinates[1] are from mesh 1, verticesCoordinates[2] are from intersections


		//for efficiency purposes, we will store the bounding-box of the Triangles
		//the bounding-boxes are defined as the extreme vertices of the triangles
		//For example, if inputTrianglesBoundingBox[0][0][0][2] = vertex with id 6 --> the "smaller" vertex of the bouding box will have z coordinate similar to the z coordinate of the
		//vertex with id 6
		vector<array<array<int,3>,2>> inputTrianglesBoundingBox[2]; 

		Point &getCoordinates(const Vertex &v) {
			return verticesCoordinates[v.getMeshId()][v.getId()];
		}

		void initTriangleBoundingBox(int meshId, int triangleId);

		Point meshBoundingBoxes[2][2]; //boundingbox of the two input meshes
		Point boundingBoxTwoMeshesTogetter[2];
		void loadInputMesh(int meshId, const string &path);

		struct TempVarsComputeIntersections {
			Point D,isectpointA1,isectpointA2,isectpointB1,isectpointB2; 
			TempVarsComputePlaneEquation tempVarsComputeEquation;
			array<VertCoord,2> isect1,isect2;
			VertCoord du0,du1,du2, dv0,dv1,dv2, vp0,vp1,vp2, up0,up1,up2,b,c,max,tmp,diff;
			VertCoord tempRationals[6];	
		};
		int intersectTwoTriangles(const InputTriangle &triMesh0,const InputTriangle &triMesh1,
				     Point &coordsPt1,VertexFromIntersection &vertexThatCreatedPt1, Point &coordsPt2,
             VertexFromIntersection &vertexThatCreatedPt2, TempVarsComputeIntersections &tempVars);



		//debugging purposes...
		Point computePointFromIntersectionVertex(VertexFromIntersection &v);



		



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
void readGTSFile(string fileName, vector<Point> &vertices,vector<InputTriangle> &triangles, Point boundingBox[2],const int meshId);


//Reads a Lium file, fills the boundingBox with the boundingBox of the triangles read
void readLiumFile(string fileName, vector<Point> &vertices,vector<InputTriangle> &triangles, Point boundingBox[2],const int meshId);

#endif
