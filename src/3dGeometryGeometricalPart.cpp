
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
  cerr << "TODO: use pre-computed normals!\n" << endl;
  cerr << "TODO: pre-compute!\n" << endl;
  assert(false); //just to remember...

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