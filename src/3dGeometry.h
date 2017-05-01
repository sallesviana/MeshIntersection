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
#include <parallel/algorithm>
#include "rationals.h"
#include "utils.h"


using namespace std;


#define COLLECT_GEOMETRY_STATISTICS


//if defined, we will double check the results with the one obtained by the original functions 
//created using mathematica (we tried to optimize these functions... thus, the double check is good
//to make sure we optimized correctly...)
#define DOUBLE_CHECK_SOS_PREDICATES_WITH_MATHEMATICA


//===============================================================
// Definitions...
//===============================================================

//typedef rational VertCoord;
//#define slightlyMoreThanOne VertCoord("10001/10000")

struct Nested3DGridWrapper ;

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
class SosPredicatesImpl;
class MeshIntersectionGeometry;

class VertexFromIntersection;
class InputVertex;

class Vertex {
	public:
		virtual bool isInputVertex() const = 0;

		const int getId() const { //different vertices from the intersection (mesh 3) may have the same id
			return id;
		}
		const int getMeshId() const {
			return meshId;
		}
		bool isFromIntersection() const {
			return getMeshId()==2; //mesh id 2 represents vertices generated from intersection.
		}


		Vertex(int meshId_,int id_): meshId(meshId_),id(id_){}
		Vertex() : meshId(-1), id(-1) {}

		friend class MeshIntersectionGeometry;

		virtual int compare(const Vertex &v) const = 0;
		virtual int compare(const InputVertex &v) const = 0;
		virtual int compare(const VertexFromIntersection &v) const = 0;

		virtual void print() const {
			cerr << meshId << " " << id ;
		}

	protected:
		int id;
		int meshId;	
};

class InputVertex: public Vertex {
	public:		
		bool isInputVertex() const { return true; }

		InputVertex(int meshId_,int id_): Vertex(meshId_,id_) {}
		InputVertex(): Vertex(-1,-1){}

		//a.compare(b) returns 0 if equal, <0 if a<b, >0 if a>b
		int compare(const InputVertex &v) const {
			//cerr << "Comparing InputVertex with InputVertex vertex" << endl;

			if(meshId!=v.meshId) return meshId-v.meshId;
			return id-v.id;
		}
		//let's suppose that all input vertices are > vertices from intersection...
		//vertices from intersection are smaller than input vertices...
		int compare(const VertexFromIntersection &v) const {
			//cerr << "Comparing InputVertex with VertexFromIntersection vertex" << endl;

			return 1;
		}
		//let's use polymorphism to determine what type of vertex v is... (double dispatch)
		int compare(const Vertex &v) const {
			//cerr << "Comparing InputVertex with Vertex vertex" << endl;

			return -v.compare(*this);
		}
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

	virtual const Vertex* getVertex(const int i) const = 0;
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


	/*
	//we assume the input triangles are unique
	bool operator==(const InputTriangle &t) const {
		if(above!=t.above || below!=t.below) return false;
		return p[0]==t.p[0] && p[1]==t.p[1] && p[2] == t.p[2];
	}

	bool operator!=(const InputTriangle &t) const {
		return !((*this)==t);
	}

	//for sorting purposes...
	bool operator<(const InputTriangle &t) const {
		if(above!=t.above) return above<t.above;
		if(below!=t.below) return below<t.below;
		if(p[0]!=t.p[0]) return p[0]<t.p[0];
		if(p[1]!=t.p[1]) return p[1]<t.p[1];
		return p[2]<t.p[2];
	}

	bool operator>=(const InputTriangle &t) const {
		return !(*this < t);
	}*/

	int compare(const InputTriangle &t) const {
		if(above!=t.above) return above-t.above;
		if(below!=t.below) return below-t.below;

		int compare0 = p[0].compare(t.p[0]);
		if(compare0!=0) return compare0;

		int compare1 = p[1].compare(t.p[1]);
		if(compare1!=0) return compare1;

		return p[2].compare(t.p[2]);
	}



	int getMeshId() const {
		return p[0].getMeshId();
	}

  const InputVertex* getInputVertex(const int i) const { 
  	return &p[i]; 
  }
  const Vertex* getVertex(const int i) const { return &p[i]; }

  virtual void print() const {
  	cerr << "Tri from: " << p[0].getMeshId() << " " << p[0].getId() << " "<< p[1].getId() << " "<< p[2].getId() << "";
  }
private:
	InputVertex p[3];
};


class RetesselationTriangle: public Triangle {
public:
	RetesselationTriangle(const Vertex &p0, const Vertex &p1, const Vertex &p2,ObjectId above, ObjectId below)
			:Triangle(above,below) {
		//cerr << p0 << " " << p1 << " " << p2 << " " << above << " " << below << endl;
		p[0] = &p0;
		p[1] = &p1;
		p[2] = &p2;		
	}

  const Vertex* getVertex(const int i) const { return p[i]; }
private:
	const Vertex *p[3];
};


//Two vertices from the intersection may have the same ids (the ids are used only to help finding the positions of the
//vertex coordinates in arrays..)
class VertexFromIntersection: public Vertex {
	public:
		bool isInputVertex() const { return false; }

		VertexFromIntersection(const InputVertex &edgeV1, const InputVertex &edgeV2, const InputTriangle triangleGeneratedVertex) {
			setEdges(edgeV1,edgeV2);
			triangle = triangleGeneratedVertex;
		}
		VertexFromIntersection() { };

		/*
		bool operator==(const VertexFromIntersection &v) const {
			cerr << "Using operator ==..." << endl;
			bool sameGeometry = Vertex::operator==(v);
			if(!sameGeometry) return false;
			return (edge[0]==v.edge[0] && edge[1]==v.edge[1]) && (triangle==v.triangle);
		}

		bool operator<(const VertexFromIntersection &v) const {
			if(getMeshId()!=v.getMeshId()) return getMeshId()<v.getMeshId();
      if(getId()!=v.getId()) return getId()<v.getId();
      bool diffEdges = !(edge[0]==v.edge[0]);
      if(diffEdges) return edge[0]<v.edge[0];
      diffEdges = !(edge[1]==v.edge[1]);
      if(diffEdges) return edge[1]<v.edge[1];
      return triangle<v.triangle;
		}

		bool operator>=(const VertexFromIntersection &v) const {
			return !(*this < v);
		}*/
		

		//let's suppose that all input vertices are > vertices from intersection...
		int compare(const InputVertex &v) const {
			//cerr << "Comparing vertex inter with input vertex" << endl;
			return -1;
		}		
		//a.compare(b) returns 0 if equal, <0 if a<b, >0 if a>b
		int compare(const VertexFromIntersection &v) const {
			//cerr << "Comparing vertex inter with VertexFromIntersection vertex" << endl;

			if(getMeshOfTriangleDefiningVertex()!=v.getMeshOfTriangleDefiningVertex()) 
				return getMeshOfTriangleDefiningVertex()-v.getMeshOfTriangleDefiningVertex();
			
      if(getId()!=v.getId()) return getId()-v.getId();


      int edgeCompare = (edge[0].compare(v.edge[0]));
      if(edgeCompare!=0) return edgeCompare;

      edgeCompare = (edge[1].compare(v.edge[1]));
      if(edgeCompare!=0) return edgeCompare;

      return triangle.compare(v.triangle);
		}


		//let's use polymorphism to determine what type of vertex v is...
		int compare(const Vertex &v) const {
			//cerr << "Comparing vertex inter with Vertex vertex, meshId: " << v.getMeshId() << endl;
			const VertexFromIntersection &thisObj = *this;
			return -v.compare(thisObj);
		}

		void print() const {
			assert(edge[0].getId()<=edge[1].getId());
			cerr << "Intersection: { (";
			Vertex::print();
			cerr << ") ";
			triangle.print();
			cerr << " -- ";
			edge[0].print();
			cerr << ",";
			edge[1].print();
			cerr << "} ";
		}


		int getMeshOfEdgeDefiningVertex() const {
			return edge[0].getMeshId();
		}

		int getMeshOfTriangleDefiningVertex() const {
			return 1-getMeshOfEdgeDefiningVertex();
		}

		void setEdges(const InputVertex &edgeV1, const InputVertex &edgeV2) {
			edge[0] = edgeV1;
			edge[1] = edgeV2;
			if(edge[0].getId()>edge[1].getId()) swap(edge[0],edge[1]);
			assert(edge[0].getId()<=edge[1].getId());
		}


	//private:

		//this vertex is the intersection between the edge and the input triangle triangle 
		InputVertex edge[2];
		InputTriangle triangle;
};



//The same pool (even in different functions) cannot be used simultaneously by different threads...
class MeshIntersectionGeometry {
	public:
		MeshIntersectionGeometry(const string &pathMesh0, const string &pathMesh1);

		//the meshId can be either 0 or 1
		struct TempVarsGetGridCellContainingVertex { VertCoord tempVertCoords; big_int tempVarsInt[3];};
		array<int,3> getGridCellContainingVertex(const int meshId, const int vertexId, const VertCoord &cellScale, TempVarsGetGridCellContainingVertex &tempVars ) const;
		

		vector<InputTriangle> inputTriangles[2];
		int getVertexIdInputTriangleWithMinCoord(int meshId, int triangleId, int coord) const { return inputTrianglesBoundingBox[meshId][triangleId][0][coord];}
		int getVertexIdInputTriangleWithMaxCoord(int meshId, int triangleId, int coord) const { return inputTrianglesBoundingBox[meshId][triangleId][1][coord];}
		int getNumVertices(int meshId) const {return verticesCoordinates[meshId].size();}
	

		~MeshIntersectionGeometry();



		struct TempVarsGetPlaneTriangleIsNotPerpendicular {VertCoord tempRationals[10];};
		int getPlaneTriangleIsNotPerpendicular(const InputTriangle &t, TempVarsGetPlaneTriangleIsNotPerpendicular &tempVars);

		void printBoundingBoxes();
		array<VertCoord,3> coordRangeMeshes() const; //"width" of the coordinates of the two meshes togetter (i.e., size of the common bounding-box in each coordinate)


		//Given pairs of triangles to test for intersection, returns (fills) the edges from intersection and the pairs of triangles that intersect
		//verticesCoordinates[2] will be filled with the coordinates of vertices generated by intersections
		void computeIntersections(const vector<pair<InputTriangle *,InputTriangle *> > &inputTrianglesToConsider, vector< pair<InputTriangle *,InputTriangle *> >  &intersectingTrianglesThatGeneratedEdges, vector< pair<VertexFromIntersection, VertexFromIntersection> >  &edgesFromIntersection, unsigned long long &numIntersectionTests);	

		
		struct TempVarsIsCloser {VertCoord tempVertCoords[4];};
		bool isCloser(const InputVertex &origV, const VertexFromIntersection &v1V, const VertexFromIntersection &v2V, TempVarsIsCloser &tempVars) const;




		//Supposing all vertices are projected to a plane
		//we will have two vectors (v1V-origV) and (v2V-origV).
		//Is the angle the second vector larger than the angle of the first one? (supposing the positive part of y=0 is the angle 0 (when the plane to project is the z=0) )
		struct TempVarsIsAngleWith0Greater {VertCoord v1x,v2x,v1y,v2y;};

		struct TempVarsIsTriangleClockwisedOriented { VertCoord tempCoords[0]; };
		bool isTriangleClockwisedOriented(const InputTriangle &t,const int whatPlaneProjectTo, TempVarsIsTriangleClockwisedOriented &tempVarsIsTriangleClockwisedOriented) const;
	
		struct TempVarsIsVertexTriangleProjection {};
		bool isVertexInTriangleProjection(const Vertex &v1,const Vertex &v2, const Vertex &v3, const Vertex &queryPoint,int whatPlaneProjectTo,TempVarsIsVertexTriangleProjection &tempVars);

		struct TempVarsIsVertexConvex { VertCoord tempCoords[2];};
		bool isVertexConvex(const Vertex &v1,const Vertex &queryVertex, const Vertex &v3,int whatPlaneProjectTo,TempVarsIsVertexConvex &tempVars);
	
		//given two vertices, do they intersect (except at endpoints) ?
		struct TempVarsDoIntersect { };
		bool doIntersect(const pair<const Vertex *,const Vertex *> &e1, const pair<const Vertex *,const Vertex *> &e2, int whatPlaneProjectTriangleTo, TempVarsDoIntersect &tempVars) const;
		bool onSegment(const Vertex & p, const Vertex & q, const Vertex & r, int whatPlaneProjectTo) const;
		

		//**************************************************************************************//
		//------------------ For PinMesh...
		
		//Checks if a vertex if projected to z=0 is in the projection of the triangle to z = 0
		struct TempVarsIsVertexTriangleProjectionZ0 { VertCoord tempVertCoords[6]; };
		bool isVertexInTriangleProjection(const InputTriangle &t, const InputVertex &queryPoint,TempVarsIsVertexTriangleProjectionZ0 &tempVars) const;

		// checks if a triangle is above a point
		//the projection of the point to z=0 is on the projection of the triangle to z=0
		//thus, we need to check if the projection of the point onto the triangle is above the point
		struct TempVarIsTriangleAbovePointSoS {};
		bool isTriangleAbovePointSoS(const InputTriangle &t, const InputVertex &p,TempVarIsTriangleAbovePointSoS &tempVars) const;

		struct TempVarIsTriangleNormalPointingPositiveZ {};
		bool isTriangleNormalPointingPositiveZ(const InputTriangle &t, TempVarIsTriangleNormalPointingPositiveZ &tempVars) const;

		//Given two triangles above a point, where the height above point is equal for both triangles, decide which one is lower according after SoS
		//the point may be even on the triangles (below, after SoS)
		struct TempVarGetBestTrianglePointInObjectSoS {};
		const InputTriangle * getBestTrianglePointInObjectSoS(const InputTriangle *candidateTriangle,const InputTriangle *bestTriangle, const InputVertex &p,TempVarGetBestTrianglePointInObjectSoS &tempVars) const;
		
		
		struct HeightPointInTriangleProjection {  private: VertCoord height; friend class MeshIntersectionGeometry;};
		


		// I think we do not need SoS here! if the point is exactly on the boundary we already consider it is in the cell above (never in the cell below)
		// If after SoS it should be in the cell below --> no problem! the result will be still correct! (will only take slightly more time to be computed)
		// Notice that the result of these functions may be wrong! the point may actually be in the cell below after SoS
		// However, the wrong result never make the algorithm wrong (only may make it slower since these two functions are only employed to determine when
		// PinMesh should stop when the cells are processed to find the lowest triangle above a point)
		struct TempVarZCellGlobalFromProjectionOfPoint {VertCoord tempVertCoord; big_int tempVarsInt[2];};
		int zCellGlobalFromProjectionOfPoint(const HeightPointInTriangleProjection &heightAbovePoint, const InputTriangle &triangle, const InputVertex &p, const Nested3DGridWrapper &uniformGrid, TempVarZCellGlobalFromProjectionOfPoint &tempVars) const;
		
		struct TempVarZCellFromProjectionOfPoint{VertCoord tempVertCoord; big_int tempVarsInt[2];};
		//In what uniform grid cell
		int zCellLevel1FromProjectionOfPoint(const HeightPointInTriangleProjection &heightAbovePoint, const InputTriangle &triangle, const InputVertex &p, const Nested3DGridWrapper &uniformGrid, TempVarZCellFromProjectionOfPoint &tempVars) const;

		
		//No need for SoS here (performance workaround...):
		//For consistency, we try to "hide" coordinates inside the HeightPointInTriangleProjection struct (so that PinMesh, like other classes, does not have access to coordinates)
		
		struct TempVarComputeHeightAbovePointNoSoS { VertCoord tempVertCoords[5]; Point vec[2]; };
		void computeHeightAbovePointNoSoS(HeightPointInTriangleProjection &height,const InputTriangle &triangle, const InputVertex &p, TempVarComputeHeightAbovePointNoSoS &tempVars) const;
		//1 if p is lower than the point with heightOtherPoint (considering the z-coordinate)
    //0 if equal
		int compareHeightWithPointHeightNoSoS(const InputVertex &queryPoint,const HeightPointInTriangleProjection &heightOtherPoint) const;
		
		//1 if height is smaller than the best, 0 if equal, -1 otherwise
		int compareHeightWithBestHeightNoSoS(const HeightPointInTriangleProjection &heightAbovePointObj,const HeightPointInTriangleProjection &bestHeightAbovePoint) const;


		//-------------------------------- End PinMesh

		struct TempVarsIsOnZeroPlusAxisNoSoS { VertCoord vecLen[2]; };
		int isOnZeroPlusAxisNoSoS(const Vertex &v1,const Vertex &v2,const int whatPlaneProjectTo, TempVarsIsOnZeroPlusAxisNoSoS &tempVars) const;


		bool isAngleWith0GreaterNoSoSNonZeroAngle(const Vertex &origV, const Vertex &v1V, const Vertex &v2V, const int planeToProject, TempVarsIsAngleWith0Greater &tempVars) const;


		//Sorting edges by angle with 0...
		struct TempVarsSortEdgesByAngle { TempVarsIsAngleWith0Greater tempVarsIsAngleWith0Greater; TempVarsIsOnZeroPlusAxisNoSoS tempVarsIsOnZeroPlusAxisNoSoS; };
		void sortEdgesSharingStartingVertexByAngle(vector<pair<const Vertex *,const Vertex *> >::iterator begin,
                                            vector<pair<const Vertex *,const Vertex *> >::iterator end,
                                            const int planeProjectTriangleTo,  TempVarsSortEdgesByAngle &tempVars) const;

		//

	
		void storeAllVertices(ostream &out); //the coordinates of vertices from mesh 0, mesh 1 and from intersection are stored in a stream
	//TODO: for debugging purposes...

		

		//for debugging purposes...
		void storeEdgesAsGts(const string &path,const vector<pair<array<double,3>,array<double,3>> > &edgesToStore) const;
		void saveEdgesAsGTS(const vector<pair<const Vertex *,const Vertex *>>  &edges,const string &path) const;

		void printPointForDebugging(const Vertex &v) const {
			array<double,3> coords = getCoordinatesForDebugging(v);
			cerr << "[" << coords[0] << "," << coords[1] << "," << coords[2] << "]";
		}
	private:
		friend class SosPredicatesImpl;
		friend class OriginalAlgFromMathematicaSosPredicatesImpl;

		struct PlaneEquation {Point normal; VertCoord d;};
		vector<PlaneEquation> planeEquationsInputTriangles[2];
		vector<int> isPlaneEquationInputTrianglesInitialized[2];

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

		const Point &getCoordinates(const Vertex &v) const {
			return verticesCoordinates[v.getMeshId()][v.getId()];
		}

		const array<double,3> getCoordinatesForDebugging(const Vertex &v) const {
			const Point &p = getCoordinates(v);
			array<double,3> pDouble;
			for(int i=0;i<3;i++) pDouble[i] = p[i].get_d();
		  return pDouble;
		}

		void initTriangleBoundingBox(int meshId, int triangleId);

		Point meshBoundingBoxes[2][2]; //boundingbox of the two input meshes
		Point boundingBoxTwoMeshesTogetter[2];
		void loadInputMesh(int meshId, const string &path);


		Point computePointFromIntersectionVertex(VertexFromIntersection &v);


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


		

		//We are actually not using these functions...
		/*
		int getGridCellXContainingVertex(int meshId,const VertCoord &x, const VertCoord &cellScale, TempVarsGetGridCellContainingVertex &tempVars  ) const;
		int getGridCellYContainingVertex(int meshId,const VertCoord &y, const VertCoord &cellScale, TempVarsGetGridCellContainingVertex &tempVars  ) const;
		int getGridCellZContainingVertex(int meshId,const VertCoord &z, const VertCoord &cellScale, TempVarsGetGridCellContainingVertex &tempVars  ) const;
		*/
		

		//We will have two versions of each predicate:
		//bool function() //the public function, no coincidency is allowed to happen
		//int functionMainImpl() //the (private) actual implementation , coincidencies (0) can happen, 1 is true, -1 is false
		//bool functionSoSImpl() //the (private) implementation using SoS...


		//this works only when the angle is >0 and there is no degenerate edge..
		int isAngleWith0GreaterNonZeroAngleMainImpl(const Vertex &origV, const Vertex &v1V, const Vertex &v2V, const int planeToProject, TempVarsIsAngleWith0Greater &tempVars) const;


		//Implementation of geometrical predicates (these implementations do not handle degeneracies --> they return 0 for degenerate cases)
		int isCloserMainImpl(const InputVertex &origV, const VertexFromIntersection &v1V, const VertexFromIntersection &v2V, TempVarsIsCloser &tempVars) const;
		int isAngleWith0GreaterMainImpl(const Vertex &origV, const Vertex &v1V, const Vertex &v2V, const int planeToProject, TempVarsIsAngleWith0Greater &tempVars) const;
		int isVertexInTriangleProjectionMainImpl(const Vertex &v1,const Vertex &v2, const Vertex &v3, const Vertex &queryPoint,int whatPlaneProjectTrianglesTo,TempVarsIsVertexTriangleProjection &tempVars) const;
		int isVertexInTriangleProjectionMainImpl(const InputTriangle &t, const InputVertex &queryPoint,TempVarsIsVertexTriangleProjectionZ0 &tempVars) const;
		int isTriangleNormalPointingPositiveZMainImpl(const InputTriangle &t, TempVarIsTriangleNormalPointingPositiveZ &tempVars) const;
		int isVertexConvexMainImpl(const Vertex &v1,const Vertex &queryVertex, const Vertex &v3,int whatPlaneProjectTrianglesTo,TempVarsIsVertexConvex &tempVars) const;
		int intersectTwoTrianglesMainImpl(const InputTriangle &triMesh0,const InputTriangle &triMesh1,
				     Point &coordsPt1,VertexFromIntersection &vertexThatCreatedPt1, Point &coordsPt2,
             VertexFromIntersection &vertexThatCreatedPt2, TempVarsComputeIntersections &tempVars);

		array<int,3> getGridCellContainingVertexOrig(const int meshId, const int vertexId, const VertCoord &cellScale, TempVarsGetGridCellContainingVertex &tempVars ) const ;
		bool isCloserOrig(const InputVertex &origV, const Vertex &v1V, const Vertex &v2V, TempVarsIsCloser &tempVars) const ;
		bool isAngleWith0GreaterOrig(const Vertex &origV, const Vertex &v1V, const Vertex &v2V, const int planeToProject, TempVarsIsAngleWith0Greater &tempVars) const ;
		bool isVertexInTriangleProjectionOrig(const Vertex &v1,const Vertex &v2, const Vertex &v3, const Vertex &queryPoint,int whatPlaneProjectTrianglesTo,TempVarsIsVertexTriangleProjection &tempVars) const;
		bool isVertexConvexOrig(const Vertex &v1,const Vertex &queryVertex, const Vertex &v3,int whatPlaneProjectTrianglesTo,TempVarsIsVertexConvex &tempVars) const;
		bool isVertexInTriangleProjectionOrig(const InputTriangle &t, const InputVertex &queryPoint,TempVarsIsVertexTriangleProjectionZ0 &tempVars) const ;
		//bool isTriangleAbovePointSoSOrig(const InputTriangle &t, const InputVertex &p,TempVarIsTriangleAbovePointSoS &tempVars) const ;
		bool isTriangleNormalPointingPositiveZOrig(const InputTriangle &t, TempVarIsTriangleNormalPointingPositiveZ &tempVars) const ;
		int zCellGlobalFromProjectionOfPointOrig(const HeightPointInTriangleProjection &heightAbovePoint, const InputTriangle &triangle, const InputVertex &p, const Nested3DGridWrapper &uniformGrid, TempVarZCellGlobalFromProjectionOfPoint &tempVars) const ;
		int zCellLevel1FromProjectionOfPointOrig(const HeightPointInTriangleProjection &heightAbovePoint, const InputTriangle &triangle, const InputVertex &p, const Nested3DGridWrapper &uniformGrid, TempVarZCellFromProjectionOfPoint &tempVars) const ;
		const InputTriangle * getBestTrianglePointInObjectSoSOrig(const InputTriangle *candidateTriangle,const InputTriangle *bestTriangle, const InputVertex &p,TempVarGetBestTrianglePointInObjectSoS &tempVars) const ;

		//SoS functions...
		int orientation(const Vertex &v1, const Vertex &v2, const Vertex &p, int whatPlaneProjectTrianglesTo) const ; 
		int orientation(const InputVertex &v1, const InputVertex &v2, const InputVertex &p, int whatPlaneProjectTrianglesTo) const ; 
		int orientation(const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &p, int whatPlaneProjectTrianglesTo) const ;
		int orientation(const InputVertex &v1, const VertexFromIntersection &v2, const VertexFromIntersection &p, int whatPlaneProjectTrianglesTo) const ;
		int orientation(const VertexFromIntersection &v1, const VertexFromIntersection &v2, const VertexFromIntersection &p, int whatPlaneProjectTrianglesTo) const ;	
		int orientation(const InputTriangle&t, const InputVertex &v) const;
		int orientation(const InputTriangle&t, const VertexFromIntersection &v) const;
		int orientation(const InputVertex&p1, const InputVertex&p2,const InputVertex&p3, const InputVertex &v) const;
  	int orientation(const InputVertex&p1, const InputVertex&p2,const InputVertex&p3, const VertexFromIntersection &v) const;
  

		//what is the signal of each coordinate the vector from orig to dest
		//cannot be 0 (SoS)
		int signalVectorCoord(const Vertex &orig, const Vertex &dest, int coord) const;
		int signalVectorCoordOnlyCallWhenCoincident(const InputVertex &orig, const InputVertex &dest, int coord) const;
		int signalVectorCoordOnlyCallWhenCoincident(const InputVertex &orig, const VertexFromIntersection &dest, int coord) const;
		int signalVectorCoordOnlyCallWhenCoincident(const VertexFromIntersection &orig, const VertexFromIntersection &dest, int coord) const;
		//int signalVectorCoordCanBe0(const Vertex &orig, const Vertex &dest, int coord) const;


		bool isCloserSoSImpl(const InputVertex &origV, const VertexFromIntersection &v1V, const VertexFromIntersection &v2V, TempVarsIsCloser &tempVars) const;
		
		bool isOrientationPositiveSoSImpl(const Vertex &origV, const Vertex &v1V, const Vertex &v2V, const int planeToProject, TempVarsIsAngleWith0Greater &tempVars) const;
		bool isAngleWith0GreaterSoSImpl(const Vertex &origV, const Vertex &v1V, const Vertex &v2V, const int planeToProject, TempVarsIsAngleWith0Greater &tempVars) const;
		bool isVertexInTriangleProjectionSoSImpl(const Vertex &v1,const Vertex &v2, const Vertex &v3, const Vertex &queryPoint,int whatPlaneProjectTrianglesTo,TempVarsIsVertexTriangleProjection &tempVars) const;
		bool isVertexConvexSoSImpl(const Vertex &v1,const Vertex &queryVertex, const Vertex &v3,int whatPlaneProjectTrianglesTo,TempVarsIsVertexConvex &tempVars) const;
		bool isVertexInTriangleProjectionSoSImpl(const InputTriangle &t, const InputVertex &queryPoint,TempVarsIsVertexTriangleProjectionZ0 &tempVars) const;
		bool isTriangleNormalPointingPositiveZSoSImpl(const InputTriangle &t, TempVarIsTriangleNormalPointingPositiveZ &tempVars) const;
		bool isTriangleAbovePointSoSImpl(const InputTriangle &t, const InputVertex &p,TempVarIsTriangleAbovePointSoS &tempVars) const;
		
		bool intersectTwoTrianglesSoSImpl(const InputTriangle &triMesh0,const InputTriangle &triMesh1,
				     VertexFromIntersection &vertexThatCreatedPt1, VertexFromIntersection &vertexThatCreatedPt2, TempVarsComputeIntersections &tempVars) const;

		//does edge (p1,p2) intersect the triangle?
		bool intersectEdgeWithTriangleSoSImpl(const InputTriangle &triangle, const InputVertex &p1, const InputVertex &p2, TempVarsComputeIntersections &tempVars) const;


		const InputTriangle * getBestTrianglePointInObjectSoSImpl(const InputTriangle *candidateTriangle,const InputTriangle *bestTriangle, const InputVertex &p,TempVarGetBestTrianglePointInObjectSoS &tempVars) const;

};











//used by the triangulation algorithm....
struct TriangulationVertex {
	const Vertex *  index; //index of the vertex in the "vertexSequence" vector
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


class BoundaryPolygon;

class BoundaryPolygon {
public:
	vector<const Vertex * > vertexSequence; //the vertices in the vertexSequence always end with the first vertex... (e.g.: 1,5,9,4,1)
	vector<array<const Vertex * ,3> > triangulatedPolygon; //stores a triangulated version of this polygon...
	
	vector<BoundaryPolygon *> boundaryPolygonOtherSideEdge;
	vector<pair<ObjectId,ObjectId> > objectsOtherMeshBoundedByThisEdge;

	ObjectId above, below;

	BoundaryPolygon(const int whatPlaneProjectTrianglesTo_): whatPlaneProjectTrianglesTo(whatPlaneProjectTrianglesTo_), polyhedronOfOtherMeshWherePolygonIs(DONT_KNOW_ID) {}

	struct TempVarsTriangulatePolygon {
		MeshIntersectionGeometry::TempVarsIsVertexTriangleProjection tempVarsIsVertexTriangleProjection;
		MeshIntersectionGeometry::TempVarsIsVertexConvex tempVarsIsVertexConvex;
	};
	void triangulatePolygon(MeshIntersectionGeometry &geometry,TempVarsTriangulatePolygon &tempVars);

	void reverseVerticesOrder(); //reverse the order of the vertices in the vertex sequence (as a consequence, above and below are swapped)
	
	/*void printTriangles(const vector<Point> &vertices) { //for debugging purposes..
		for(const array<const Vertex *,3>&t:triangulatedPolygon) {
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
	}*/


	void setPolyhedronWherePolygonIs(ObjectId polyId) {
		polyhedronOfOtherMeshWherePolygonIs = polyId;
	}
	ObjectId getPolyhedronWherePolygonIs() const {
		return polyhedronOfOtherMeshWherePolygonIs;
	}
private:
	ObjectId polyhedronOfOtherMeshWherePolygonIs;

	int whatPlaneProjectTrianglesTo; //to what plane can we project this polygon without creating a degenerate polygon?

	bool pointInTriangleProj(MeshIntersectionGeometry &geometry,int p0,int p1, int p2, int queryPoint,TempVarsTriangulatePolygon &tempVars);
	void initializeLinkedList(MeshIntersectionGeometry &geometry,TriangulationVertex listVerticesToProcess[],int numVerticesPolygon,
																						 int &eBegin,int &eEnd,
																						int &cBegin,int &cEnd,
																						int &rBegin,int &rEnd,TempVarsTriangulatePolygon &tempVars);
	void updateStatusVertex(MeshIntersectionGeometry &geometry,int vertexId,TriangulationVertex listVerticesToProcess[], 
													int &rBegin,int &rEnd, int &eBegin,TempVarsTriangulatePolygon &tempVars);

	//vertexId is the position of the vertex in the "listVerticesToProcess" arrat...
	bool isEar(MeshIntersectionGeometry &geometry,int vertexId,TriangulationVertex listVerticesToProcess[], const int rBegin,TempVarsTriangulatePolygon &tempVars);
	
	//TempCoords should have at least 2 coordinates
	bool isConvex(MeshIntersectionGeometry &geometry,int vertexId,TriangulationVertex listVerticesToProcess[], TempVarsTriangulatePolygon &tempVars);


};



//Reads a GTS file, fills the boundingBox with the boundingBox of the triangles read
void readGTSFile(string fileName, vector<Point> &vertices,vector<InputTriangle> &triangles, Point boundingBox[2],const int meshId);


//Reads a Lium file, fills the boundingBox with the boundingBox of the triangles read
void readLiumFile(string fileName, vector<Point> &vertices,vector<InputTriangle> &triangles, Point boundingBox[2],const int meshId);

#endif
