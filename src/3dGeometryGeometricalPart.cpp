
array<VertCoord,3> MeshIntersectionGeometry::coordRangeMeshes() const {
  return {boundingBoxTwoMeshesTogetter[1][0]-boundingBoxTwoMeshesTogetter[0][0],boundingBoxTwoMeshesTogetter[1][1]-boundingBoxTwoMeshesTogetter[0][1],boundingBoxTwoMeshesTogetter[1][2]-boundingBoxTwoMeshesTogetter[0][2]};
} 

array<int,3> MeshIntersectionGeometry::getGridCellContainingVertex(const int meshId, const int vertexId, const VertCoord &cellScale, TempVarsGetGridCellContainingVertex &tempVars ) const {
	VertCoord &tempVar = tempVars.tempVertCoords;//tempVertCoords[poolToUse];
	const vector<Point> &points = verticesCoordinates[meshId];

  const Point &boundingBoxMin = boundingBoxTwoMeshesTogetter[0];

	tempVar = points[vertexId][0];
  tempVar -= boundingBoxMin[0];
  tempVar *= cellScale;
  const int x = convertToInt(tempVar,tempVars.tempVarsInt);

  tempVar = points[vertexId][1];
  tempVar -= boundingBoxMin[1];
  tempVar *= cellScale;
  const int y = convertToInt(tempVar,tempVars.tempVarsInt);

  tempVar = points[vertexId][2];
  tempVar -= boundingBoxMin[2];
  tempVar *= cellScale;
  const int z  = convertToInt(tempVar,tempVars.tempVarsInt);

  return {x,y,z};
}

int MeshIntersectionGeometry::getGridCellXContainingVertex(int meshId,const VertCoord &xCoord, const VertCoord &cellScale, TempVarsGetGridCellContainingVertex &tempVars ) const {
  VertCoord &tempVar = tempVars.tempVertCoords;
  const Point &boundingBoxMin = boundingBoxTwoMeshesTogetter[0];

  tempVar = xCoord;
  tempVar -= boundingBoxMin[0];
  tempVar *= cellScale;
  const int x = convertToInt(tempVar,tempVars.tempVarsInt);

  return x;
}
int MeshIntersectionGeometry::getGridCellYContainingVertex(int meshId,const VertCoord &yCoord, const VertCoord &cellScale, TempVarsGetGridCellContainingVertex &tempVars ) const {
  VertCoord &tempVar = tempVars.tempVertCoords;
  const Point &boundingBoxMin = boundingBoxTwoMeshesTogetter[0];

  tempVar = yCoord;
  tempVar -= boundingBoxMin[1];
  tempVar *= cellScale;
  const int y = convertToInt(tempVar,tempVars.tempVarsInt);

  return y;
}
int MeshIntersectionGeometry::getGridCellZContainingVertex(int meshId,const VertCoord &zCoord, const VertCoord &cellScale, TempVarsGetGridCellContainingVertex &tempVars ) const {
  VertCoord &tempVar = tempVars.tempVertCoords;
  const Point &boundingBoxMin = boundingBoxTwoMeshesTogetter[0];

  tempVar = zCoord;
  tempVar -= boundingBoxMin[2];
  tempVar *= cellScale;
  const int z  = convertToInt(tempVar,tempVars.tempVarsInt);

  return z;
}



//Given the vertices of a triangle, compute the plane equation ( N1.x + d = 0)
void MeshIntersectionGeometry::computePlaneEquation(PlaneEquation &equation, const Point &V0, const Point &V1, const Point &V2, TempVarsComputePlaneEquation &tempVars) {
  /* compute plane equation of triangle(V0,V1,V2) */
  SUB(tempVars.E1,V1,V0);
  SUB(tempVars.E2,V2,V0);
  CROSS(equation.normal,tempVars.E1,tempVars.E2,tempVars.temp);
  //d1=-DOT(N1,V0);
  MinusDOT(equation.d,equation.normal,V0,tempVars.temp);  

  /* plane equation 1: N1.X+d1=0 */
}


int MeshIntersectionGeometry::getPlaneTriangleIsNotPerpendicular(const InputTriangle &t, TempVarsGetPlaneTriangleIsNotPerpendicular &tempVariables) {

  VertCoord*tempVars = tempVariables.tempRationals;

  const Point &p1 = getCoordinates(*t.getInputVertex(0));
  const Point &p2 = getCoordinates(*t.getInputVertex(1));
  const Point &p3 = getCoordinates(*t.getInputVertex(2));

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

//Is v1 closer to origV than v2 is?
bool MeshIntersectionGeometry::isCloser(const InputVertex &origV, const Vertex &v1V, const Vertex &v2V, TempVarsIsCloser &tempVars) const {
  const Point &orig = getCoordinates(origV);
  const Point &v1 = getCoordinates(v1V);
  const Point &v2 = getCoordinates(v2V);
  VertCoord *tempVertCoords = tempVars.tempVertCoords;

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
}


bool MeshIntersectionGeometry::isAngleWith0Greater(const Vertex &origV, const Vertex &v1V, const Vertex &v2V, const int planeToProject, TempVarsIsAngleWith0Greater &tempVars) const {
  int xCoord = 0;
  int yCoord = 1;
  if(planeToProject==PLANE_Y0) {
    xCoord= 2;
    yCoord= 0;
  } else if(planeToProject==PLANE_X0) {
    xCoord= 1;
    yCoord= 2;    
  } 

  const Point &e10 = getCoordinates(origV);
  const Point &e11 = getCoordinates(v1V);

  //e20 = e10, since the two edges share the origin...
  const Point &e21 = getCoordinates(v2V);

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
  VertCoord &v1x = tempVars.v1x;
  VertCoord &v1y = tempVars.v1y;
  VertCoord &v2x = tempVars.v2x;
  VertCoord &v2y = tempVars.v2y;
  v1x = e11[xCoord];
  v1x -= e10[xCoord];

  v1y = e11[yCoord];
  v1y -= e10[yCoord];

  v2x = e21[xCoord];
  v2x -= e10[xCoord];

  v2y = e21[yCoord];
  v2y -= e10[yCoord];

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

bool MeshIntersectionGeometry::isTriangleClockwisedOriented(const InputTriangle &t,const int whatPlaneProjectTo, TempVarsIsTriangleClockwisedOriented &tempVars) const {
  return signCrossProduct2D(getCoordinates(*t.getInputVertex(0)),getCoordinates(*t.getInputVertex(1)),getCoordinates(*t.getInputVertex(0)),getCoordinates(*t.getInputVertex(2)),whatPlaneProjectTo,tempVars.tempCoords)<0;
}

bool MeshIntersectionGeometry::isVertexInTriangleProjection(const Vertex &v1,const Vertex &v2, const Vertex &v3, const Vertex &queryPoint,int whatPlaneProjectTrianglesTo,TempVarsIsVertexTriangleProjection &tempVars) {
  const Point &p = getCoordinates(queryPoint);
  const Point &p0 =  getCoordinates(v1);
  const Point &p1 =  getCoordinates(v2);
  const Point &p2 =  getCoordinates(v3);


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
  
  
  return 0 <= c && c <= 1; 
}

bool MeshIntersectionGeometry::isVertexConvex(const Vertex &v1,const Vertex &queryVertex, const Vertex &v3,int whatPlaneProjectTrianglesTo,TempVarsIsVertexConvex &tempVars) {
  const Point &vertex1 =  getCoordinates(v1);
  const Point &vertex2 =  getCoordinates(queryVertex);
  const Point &vertex3 =  getCoordinates(v3);

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

  /* //REMOVE....
  VertCoord twoArea = vertex1[coordX]*(vertex2[coordY]-vertex3[coordY]) +
                      vertex2[coordX]*(vertex3[coordY]-vertex1[coordY]) +
                      vertex3[coordX]*(vertex1[coordY]-vertex2[coordY]) ;

  listVerticesToProcess[vertexId].convex = twoArea<0;
  return listVerticesToProcess[vertexId].convex;
  */
  VertCoord *tempCoords  = tempVars.tempCoords;
  
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
  

  return sgn(tempCoords[0])<0;
}





/*************** PinMesh Part... ****************/

bool MeshIntersectionGeometry::isVertexInTriangleProjection(const InputTriangle &t, const InputVertex &queryPoint,TempVarsIsVertexTriangleProjectionZ0 &tempVars) const {
  const Point &p0 = getCoordinates(*t.getInputVertex(0));
  const Point &p1 = getCoordinates(*t.getInputVertex(1));
  const Point &p2 = getCoordinates(*t.getInputVertex(2));

  const Point &p = getCoordinates(queryPoint);

  if ( p==p0 || p==p1 || p==p2)  {
    cerr << "TODO: SoS vertex in triangle projection 1" << endl;
    return false;
    //return pointInTriangleProjSoS(p0, p1, p2, p) ;//return 0; // degenerate case --> SoS
  }
  
  VertCoord * tempVertCoords = tempVars.tempVertCoords;
  tempVertCoords[0] = p1[1];
  tempVertCoords[0] -= p2[1];
  tempVertCoords[1] = p0[0];
  tempVertCoords[1] -= p2[0];
  tempVertCoords[0] *= tempVertCoords[1];

  tempVertCoords[1]  = p2[0];
  tempVertCoords[1]  -= p1[0];
  tempVertCoords[2]  = p0[1];
  tempVertCoords[2]  -= p2[1];
  tempVertCoords[1] *= tempVertCoords[2];

  tempVertCoords[0] += tempVertCoords[1]; //denominator

  //VertCoord denominator = ((p1[1] - p2[1])*(p0[0] - p2[0]) + (p2[0] - p1[0])*(p0[1] - p2[1]));
  //assert(denominator == tempVertCoords[0]);

  if (sgn(tempVertCoords[0])==0) { //if (denominator==0) { //TODO: check this....
    //cerr << "TODO: SoS vertex in triangle projection 2" << endl;
    //Because of SoS no point should be below a vertical triangle....
    return false; //pointInTriangleProjSoS(p0, p1, p2, p) ; //SoS
    //return 0;
  }
  tempVertCoords[1] = p1[1];
  tempVertCoords[1] -= p2[1];

  tempVertCoords[2] = p[0];
  tempVertCoords[2] -= p2[0];

  tempVertCoords[1] *= tempVertCoords[2];

  tempVertCoords[3] = p2[0];
  tempVertCoords[3] -= p1[0];

  tempVertCoords[4] =  p[1];
  tempVertCoords[4] -= p2[1];

  tempVertCoords[3] *= tempVertCoords[4];

  tempVertCoords[1] += tempVertCoords[3];

  tempVertCoords[1] /= tempVertCoords[0]; // a

  //VertCoord a = ((p1[1] - p2[1])*(p[0] - p2[0]) + (p2[0] - p1[0])*(p[1] - p2[1])) / denominator;
  //assert(a == tempVertCoords[1]);

  if (sgn(tempVertCoords[1])<0 || tempVertCoords[1]>1) return false; //if ( a<0 || a >1) return false;
  
  tempVertCoords[2] = p2[1];
  tempVertCoords[2] -= p0[1];

  tempVertCoords[3] = p[0];
  tempVertCoords[3] -= p2[0];

  tempVertCoords[2] *= tempVertCoords[3];

  tempVertCoords[3] = p0[0];
  tempVertCoords[3] -= p2[0];

  tempVertCoords[4] = p[1];
  tempVertCoords[4] -= p2[1];

  tempVertCoords[3] *= tempVertCoords[4];
  tempVertCoords[2] += tempVertCoords[3];

  tempVertCoords[2] /= tempVertCoords[0]; //b
  //VertCoord b = ((p2[1] - p0[1])*(p[0] - p2[0]) + (p0[0] - p2[0])*(p[1] - p2[1])) / denominator;
  //assert(b == tempVertCoords[2]);

  if (sgn(tempVertCoords[2]) <0 || tempVertCoords[2]>1) return false;//if (b<0 || b>1) return false;
  //VertCoord c = 1 - a - b;
  tempVertCoords[3] = 1;
  tempVertCoords[3] -= tempVertCoords[1];
  tempVertCoords[3] -= tempVertCoords[2];


  if (sgn(tempVertCoords[3]) < 0 || tempVertCoords[3] >1) return false;

  //All coordinates are >=0 && <=1... so we have either a degenerate case or the point is in the triangle...
  if (sgn(tempVertCoords[3])==0 || sgn(tempVertCoords[2])==0 || sgn(tempVertCoords[1])==0) {
    cerr << "TODO: SoS vertex in triangle projection 3" << endl;
    return false;
    //bool ans =  pointInTriangleProjSoS(p0, p1, p2, p) ; //SoS
    //return ans;
  }
  return 1;
}



// given two triangles with same real height above a point, use SoS to break tie...
bool MeshIntersectionGeometry::isTriangleAbovePointSoS(const InputTriangle &t, const InputVertex &p,TempVarIsTriangleAbovePointSoS &tempVars) const {
  cerr << "TODO: SoS is triangle above point" << endl;
  return false;
}





struct TempVarIsTriangleNormalPointingPositiveZ {};
//TODO: use pre-computed normals...
bool MeshIntersectionGeometry::isTriangleNormalPointingPositiveZ(const InputTriangle &t, TempVarIsTriangleNormalPointingPositiveZ &tempVars) const {
    const Point &p0 = getCoordinates(*t.getInputVertex(0));
    const Point &p1 = getCoordinates(*t.getInputVertex(1));
    const Point &p2 = getCoordinates(*t.getInputVertex(2));

    VertCoord vec[2][3];  //First, lets compute the cross product between the vectors representing the triangle
    vec[0][0] = p1[0]-p0[0];
    vec[0][1] = p1[1]-p0[1];
    //vec[0][2] = triangle[1][2]-triangle[0][2];

    vec[1][0] = p2[0]-p0[0];
    vec[1][1] = p2[1]-p0[1];
    //vec[1][2] = triangle[2][2]-triangle[0][2];

    //lets compute the values A,B,C basing on the two 3D vectors vec[0] and vec[1]
    //  | i             j           k     |
    //  | vec[0][0]  vec[0][1]  vec[0][2] |
    //  | vec[1][0]  vec[1][1]  vec[1][2] |

    //VertCoord A = vec[0][1]*vec[1][2] - vec[0][2]*vec[1][1];
    //VertCoord B = vec[0][0]*vec[1][2] - vec[0][2]*vec[1][0];
    VertCoord C = vec[0][0]*vec[1][1] - vec[0][1]*vec[1][0];

    //Now we have the triangle normal (A,B,C) ... it points to the "positive" direction of the triangle (that is, it has the same orientation used to define
    //what is above and below the triangle according to the right hand rule)

    //If the dot product between the triangle's normal and the z+'s vector ( (0,0,1) ) is positive  --> the point is in the volume "below" the triangle
    // if it is 0 --> error (the triangle is vertical -- TODO treat this... maybe with SoS)
    //if is negative the point is in the volume above the triangle 
    
    //we actually only need to test if C is negative, 0 or positive...

    if (C==0) {
      cerr << "Error... vertical triangle (TODO: fix this with SoS)" << endl;
      exit(1);
    }
    return C>0; 
}



#include "nested3DGrid.h"
struct HeightPointInTriangleProjection {};
struct TempVarZCellGlobalFromProjectionOfPoint {};

//TODO: use input point/triangle to correctly classify according to SoS in case of coincidences...
int MeshIntersectionGeometry::zCellGlobalFromProjectionOfPoint(const HeightPointInTriangleProjection &heightAbovePoint, const InputTriangle &triangle, const InputVertex &p, const Nested3DGridWrapper &uniformGrid, TempVarZCellGlobalFromProjectionOfPoint &tempVars) const {
  VertCoord &tempVar = tempVars.tempVertCoord;
  const Point &boundingBoxMin = boundingBoxTwoMeshesTogetter[0];

  tempVar = heightAbovePoint.height;
  tempVar -= boundingBoxMin[2];
  tempVar *= uniformGrid.cellScale2Levels;
  const int z  = convertToInt(tempVar,tempVars.tempVarsInt);

  return z;
}




int MeshIntersectionGeometry::zCellLevel1FromProjectionOfPoint(const HeightPointInTriangleProjection &heightAbovePoint, const InputTriangle &triangle, const InputVertex &p, const Nested3DGridWrapper &uniformGrid, TempVarZCellFromProjectionOfPoint &tempVars) const {
  VertCoord &tempVar = tempVars.tempVertCoord;
  const Point &boundingBoxMin = boundingBoxTwoMeshesTogetter[0];

  tempVar = heightAbovePoint.height;
  tempVar -= boundingBoxMin[2];
  tempVar *= uniformGrid.cellScaleLevel1;
  const int z  = convertToInt(tempVar,tempVars.tempVarsInt);

  return z;
}


//Given two triangles above a point, where the height above point is equal for both triangles, decide which one is lower according after SoS
const InputTriangle * MeshIntersectionGeometry::getBestTrianglePointInObjectSoS(const InputTriangle *candidateTriangle,const InputTriangle *bestTriangle, const InputVertex &p,TempVarGetBestTrianglePointInObjectSoS &tempVars) const {
  cerr << "TODO: SoS get best triangle point in object SoS" << endl;
  return candidateTriangle;
}


//No need for SoS here (performance workaround...):
//For consistency, we try to "hide" coordinates inside the HeightPointInTriangleProjection struct (so that PinMesh, like other classes, does not have access to coordinates)

//No need for SoS here: the triangle will be non-degenerate (input triangle) and, thus, the height of the projection
//onto the plane is well defined...
void MeshIntersectionGeometry::computeHeightAbovePointNoSoS(HeightPointInTriangleProjection &heightAbovePoint, const InputTriangle &triangle, const InputVertex &queryPoint, TempVarComputeHeightAbovePointNoSoS &tempVars) const {
  //TODO: re-use plane equations, pre-computed normals, etc..
  const Point &p0 = getCoordinates(*triangle.getInputVertex(0));
  const Point &p1 = getCoordinates(*triangle.getInputVertex(1));
  const Point &p2 = getCoordinates(*triangle.getInputVertex(2));

  const Point &p = getCoordinates(queryPoint);

  Point *vec = tempVars.vec;

  // lets compute the cross product of two vectors of the triangle in order to get the plane equation A(x-x0)+ B(y-y0)+C(z-z0)= 0
  //more info: http://tutorial.math.lamar.edu/Classes/CalcIII/EqnsOfPlanes.aspx
  //VertCoord vec[2][3];
  vec[0][0] = p1[0]; vec[0][0] -= p0[0]; //vec[0][0] = p1[0]-p0[0];
  vec[0][1] = p1[1]; vec[0][1] -= p0[1]; //vec[0][1] = p1[1]-p0[1];
  vec[0][2] = p1[2]; vec[0][2] -= p0[2];  //vec[0][2] = p1[2]-p0[2];

  vec[1][0] = p2[0]; vec[1][0] -= p0[0]; //vec[1][0] = p2[0]-p0[0];
  vec[1][1] = p2[1]; vec[1][1] -= p0[1]; //vec[1][1] = p2[1]-p0[1];
  vec[1][2] = p2[2]; vec[1][2] -= p0[2]; //vec[1][2] = p2[2]-p0[2];

  VertCoord *tempVertCoords = tempVars.tempVertCoords;

  //lets compute the values A,B,C basing on the two 3D vectors vec[0] and vec[1]
  //  | i             j           k     |
  //  | vec[0][0]  vec[0][1]  vec[0][2] |
  //  | vec[1][0]  vec[1][1]  vec[1][2] |

  //VertCoord A = vec[0][1]*vec[1][2] - vec[0][2]*vec[1][1];
  tempVertCoords[0] = vec[0][1];
  tempVertCoords[0] *= vec[1][2];
  tempVertCoords[1] =  vec[0][2];
  tempVertCoords[1] *= vec[1][1];
  tempVertCoords[0] -= tempVertCoords[1]; //this will be A


  //wrong--->  VertCoord B = vec[0][0]*vec[1][2] - vec[0][2]*vec[1][0];
  //correct --> VertCoord B = vec[0][2]*vec[1][0] - vec[0][0]*vec[1][2];
  tempVertCoords[1] = vec[0][2];
  tempVertCoords[1] *= vec[1][0];
  tempVertCoords[2] =  vec[0][0];
  tempVertCoords[2] *= vec[1][2];
  tempVertCoords[1] -= tempVertCoords[2]; // this will be B

  //VertCoord C = vec[0][0]*vec[1][1] - vec[0][1]*vec[1][0];
  tempVertCoords[2] = vec[0][0];
  tempVertCoords[2] *= vec[1][1];
  tempVertCoords[3] = vec[0][1];
  tempVertCoords[3] *= vec[1][0];
  tempVertCoords[2] -= tempVertCoords[3]; //this will be C

  //assert(C!=0);
  assert(tempVertCoords[2]!=0);

  // now we have an equation A(x-triangle[0][0]) +  B(y-triangle[0][1]) + C(z-triangle[0][2]) = 0

  //heightAbovePoint = ( A*(p0[0]-p[0]) +  B*(p0[1]-p[1]) +  C*p0[2])/C; //TODO: check this...
  tempVertCoords[3] = p0[0];
  tempVertCoords[3] -= p[0];
  tempVertCoords[3] *= tempVertCoords[0];
  heightAbovePoint.height = tempVertCoords[3];

  tempVertCoords[3] = p0[1];
  tempVertCoords[3] -= p[1];
  tempVertCoords[3] *= tempVertCoords[1];
  heightAbovePoint.height += tempVertCoords[3];

  tempVertCoords[3] = p0[2];
  tempVertCoords[3] *= tempVertCoords[2];
  heightAbovePoint.height += tempVertCoords[3];

  heightAbovePoint.height /= tempVertCoords[2];
}

//1 if p < heightAbovePoint  (considering the z-coordinate)
//0 if equal
//-1 otherwise
int MeshIntersectionGeometry::compareHeightWithPointHeightNoSoS(const InputVertex &queryPoint,const HeightPointInTriangleProjection &heightOtherPoint) const {
  const Point &p = getCoordinates(queryPoint);
  if(heightOtherPoint.height > p[2]) return 1;
  if(heightOtherPoint.height == p[2]) return 0;
  return -1;
}

//1 if height is smaller than the best, 0 if equal, -1 otherwise
int MeshIntersectionGeometry::compareHeightWithBestHeightNoSoS(const HeightPointInTriangleProjection &heightAbovePointObj,const HeightPointInTriangleProjection &bestHeightAbovePoint) const {
  if(heightAbovePointObj.height < bestHeightAbovePoint.height) return 1;
  if(heightAbovePointObj.height == bestHeightAbovePoint.height) return 0;
  return -1;
}
