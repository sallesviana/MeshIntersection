#ifndef PINMESH_H
#define PINMESH_H

#include "utils.h"
#include "3dGeometry.h"
#include "nested3DGrid.h"



#include "floodFillScanline.h"

class PinMesh {
  //Returns true iff the projection of the point p into a horizontal plane is inside of the projection of triangle
  // (p0,p1,p2) into the horizontal plane...

  //returns true if the point p is vertically directly above/below the point p...
  bool pointInTriangleProj(const Point &p0, const Point &p1, const Point &p2, const Point &p)  ;

  //Returns true iff the projection of the point p into a horizontal plane is inside of the projection of triangle
  // (p0,p1,p2) into the horizontal plane...

  //returns true if the point p is vertically directly above/below the point p...

  // Uses SoS to treat the degenerate cases
  // We pretend the point (p) is sligtly translated translated: translation (epsilon, epsilon^2, epsilon^3) 

  
  // a_epsilon = a + (epsilon(p1y-p2y) + epsilon^2 (p2x-p1x) )/den
  // b_epsilon = b + (epsilon(p2y-p0y) + epsilon^2 (p0x-p2x) )/den
  // c_epsilon = 1 - a_epsilon - b_epsilon = 1 - a - (epsilon(p1y-p2y) + epsilon^2 (p2x-p1x) )/den - b - (epsilon(p2y-p0y) + epsilon^2 (p0x-p2x) )/den
  // c_epsilon = 1 - a - b - ( epsilon(p1y-p0y) + epsilon^2(p0x-p1x)  )/den
  bool pointInTriangleProjSoS(const Point &p0, const Point &p1, const Point &p2, const Point &p) ;

   

  //p0,p1,p2 are the vertices of the triangle...
  void getHeigthAbovePoint(VertCoord &heightAbovePoint,const Point &p0,const Point &p1,const Point &p2,const Point &p) ;

  //--------------------------------------------
  //--------------------------------------------
  //--------------------------------------------



  //Returns true iff the projection of the point p into a horizontal plane is inside of the projection of triangle
  // (p0,p1,p2) into the horizontal plane...
  //Uses SoS to treat the special cases (degeneracies...)
  //tempVertCoords shoudl have at least 5 elements (this variable is used to avoid allocating temporary memory).
  bool pointInTriangleProj(const Point &p0, const Point &p1, const Point &p2, const Point &p, VertCoord *tempVertCoords); 


  //returns true if the point p is vertically directly above/below the point p...
  

  //p0,p1,p2 are the vertices of the triangle...
  // vec is a temporary matrix of coordinates
  // tempVertCoords is a temporary array with size at least 4
  void getHeigthAbovePoint(VertCoord &heightAbovePoint,const Point &p0,const Point &p1,const Point &p2,const Point &p, VertCoord vec[2][3],VertCoord *tempVertCoords) ;

  void computeBarycentricCoordinates(const Point &p0,const Point &p1,const Point &p2, VertCoord &A, VertCoord &B, VertCoord &C) ;



  const Triangle * getBestTrianglePointInObjectSoS(int meshId,const Triangle *newTriangle,const Triangle *bestTriangleSoFar,const Point &p) ;

  bool triangleAbovePointSoS(int meshId,const Triangle &triangle,const Point &p) ;

  //lets suppose the uniform grid works with only one level...
  //given a point, the coordinates of the grid cell where p is and a 1-level uniform grid, computes the object where p is
  //if p is outside the objects, returns OUTSIDE_OBJECT
  //meshId (0 or 1) represents the map that we want to verify 


  //temp_big_ints should have size at least 2
  //tempVertCoords should have size at least 8
  struct TempVarsComputeObjectWherePointIs {
    MeshIntersectionGeometry::TempVarsIsVertexTriangleProjectionZ0 tempVarsIsVertexTriangleProjectionZ0;
    MeshIntersectionGeometry::TempVarComputeHeightAbovePointNoSoS tempVarComputeHeightAbovePointNoSoS;
    MeshIntersectionGeometry::TempVarZCellFromProjectionOfPoint tempVarZCellFromProjectionOfPoint;
    MeshIntersectionGeometry::TempVarZCellGlobalFromProjectionOfPoint tempVarZCellGlobalFromProjectionOfPoint;
    MeshIntersectionGeometry::TempVarGetBestTrianglePointInObjectSoS tempVarGetBestTrianglePointInObjectSoS;
    MeshIntersectionGeometry::TempVarIsTriangleAbovePointSoS tempVarIsTriangleAbovePointSoS;
    MeshIntersectionGeometry::TempVarIsTriangleNormalPointingPositiveZ tempVarIsTriangleNormalPointingPositiveZ;
    MeshIntersectionGeometry::HeightPointInTriangleProjection heightAbovePoint,heightAbovePointBestTriangle;
    TempVarsSoSPredicatesImpl tempVarsSoSPredicatesImpl;
  };
  ObjectId computeObjectWherePointIsTwoLevel(const InputVertex &p,int globalGridCoordX,int globalGridCoordY, int globalGridCoordZ, int meshId, TempVarsComputeObjectWherePointIs &tempVars, bool &foundUsingGrid);






  //#define PINMESH_VERBOSE




  const Nested3DGridWrapper *uniformGrid;
  MeshIntersectionGeometry *geometry;  //array with size at least 2 (one for each mesh)


public:
  PinMesh(const Nested3DGridWrapper *uniformGrid_, MeshIntersectionGeometry *geometry_): 
              uniformGrid(uniformGrid_),geometry(geometry_) {  
  }  

  //Locate a set of vertices in the objects in map 0...
  //if meshIdToLocate=0, this means that vertices will be located in mesh 0
  //#define PINMESH_VERBOSE
  void  locateVerticesInObject(const vector<InputVertex> &verticesToLocate,std::vector<ObjectId> &verticesIds,int meshIdToLocate);

};

#endif