


/*
Predicates and functions needing SoS


*/


int MeshIntersectionGeometry::intersectTwoTriangles(const InputTriangle &triMesh0,const InputTriangle &triMesh1,
             Point &coordsPt1,VertexFromIntersection &vertexThatCreatedPt1, Point &coordsPt2,
             VertexFromIntersection &vertexThatCreatedPt2, TempVarsComputeIntersections &tempVars) {

  //TODO: detect coincidencies properly here...
  //co-planar is not necessarely a coincidence...
  //TODO
  //TODO: intersection at edge/vertex --> SoS...
  int ans = intersectTwoTrianglesMainImpl(triMesh0,triMesh1,coordsPt1,vertexThatCreatedPt1, coordsPt2, vertexThatCreatedPt2, tempVars);
  
  /*VertexFromIntersection p1SoS,p2SoS;
  bool ansSoS = intersectTwoTrianglesSoSImpl(triMesh0,triMesh1,coordsPt1,p1SoS, coordsPt2, p2SoS, tempVars);

  #pragma omp critical
  {
    //cerr << ans << " " << ansSoS << endl;
    if(ans!=ansSoS) {
      cerr << "@@@@@@@@@@@@@@@@@ DIFF ans in intersect two triangles" << endl;
    }
    if(ans && ansSoS) {
      Point pOrig1 = computePointFromIntersectionVertex(vertexThatCreatedPt1);
      Point pOrig2 = computePointFromIntersectionVertex(vertexThatCreatedPt2);
      if(pOrig1>pOrig2) swap(pOrig1,pOrig2);

      Point pSoS1 = computePointFromIntersectionVertex(p1SoS);
      Point pSoS2 = computePointFromIntersectionVertex(p2SoS);
      if(pSoS1>pSoS2) swap(pSoS1,pSoS2);

      assert(pOrig1==pSoS1);
      assert(pOrig2==pSoS2);
    }
  }*/

  #ifdef COLLECT_GEOMETRY_STATISTICS
    if(ans==0) {
      #pragma omp atomic
      geometryStatisticsDegenerateCases.ctDegeneraciesIntersectTwoTriangles++;
    } else {
      #pragma omp atomic
      geometryStatisticsNonDegenerateCases.ctDegeneraciesIntersectTwoTriangles++;
    }
  #endif

  return ans==1;
}


//Is v1 closer to origV than v2 is?
bool MeshIntersectionGeometry::isCloser(const InputVertex &origV, const VertexFromIntersection &v1V, const VertexFromIntersection &v2V, TempVarsIsCloser &tempVars) const {
  int ans = isCloserMainImpl(origV, v1V, v2V, tempVars);
  

  #ifdef COLLECT_GEOMETRY_STATISTICS
    if(ans==0) {
      #pragma omp atomic
      geometryStatisticsDegenerateCases.ctDegeneraciesIsCloser++;
    } else {
      #pragma omp atomic
      geometryStatisticsNonDegenerateCases.ctDegeneraciesIsCloser++;
    }
  #endif

  if(ans==0) {
    bool ansSoS = isCloserSoSImpl(origV, v1V, v2V, tempVars);
    bool ansOrig = isCloserOrig(origV, v1V, v2V, tempVars);

    //cerr << "Coincidency in isAngleGreater...: " << ans << " " << ansSoS << endl;
    if(ansOrig!=ansSoS) {
      #pragma omp critical
      {
        cerr << "@@@ Is closer..." << endl;
        cerr << ansOrig << " " << ansSoS  << endl;
        for(int i=0;i<3;i++) { VertCoord c = getCoordinates(v1V)[i]-getCoordinates(origV)[i]; cerr << (c).get_d() << " "; } cerr << endl;
        for(int i=0;i<3;i++) { VertCoord c = getCoordinates(v2V)[i]-getCoordinates(origV)[i]; cerr << (c).get_d() << " "; } cerr << endl;
      }
    }

    return ansSoS;
  } else {
    //assert( ((ans==1)&&(ansSoS)) || ((ans==-1)&&(!ansSoS)) ); //if SoS answer is true --> ans have to be 1
    return ans==1;
  }
}


bool MeshIntersectionGeometry::isAngleWith0Greater(const Vertex &origV, const Vertex &v1V, const Vertex &v2V, const int planeToProject, TempVarsIsAngleWith0Greater &tempVars) const {
  int ans = isAngleWith0GreaterMainImpl(origV, v1V, v2V, planeToProject, tempVars);
  

  #ifdef COLLECT_GEOMETRY_STATISTICS
    if(ans==0) {
      #pragma omp atomic
      geometryStatisticsDegenerateCases.ctDegeneraciesIsAngleWith0Greater++;
    } else {
      #pragma omp atomic
      geometryStatisticsNonDegenerateCases.ctDegeneraciesIsAngleWith0Greater++;
    }
  #endif

  if(ans==0) {
    bool ansSoS = isAngleWith0GreaterSoSImpl(origV, v1V, v2V, planeToProject, tempVars);
    bool ansOrig = isAngleWith0GreaterOrig(origV, v1V, v2V, planeToProject, tempVars);
    //cerr << "Coincidency in isAngleGreater...: " << ans << " " << ansSoS << endl;
    if(ansOrig!=ansSoS) {
      #pragma omp critical
      {
        cerr << "@@@ Is angle with 0 greater..." << endl;
        cerr << ansOrig << " " << ansSoS << " " << planeToProject << endl;
        for(int i=0;i<3;i++) { VertCoord c = getCoordinates(v1V)[i]-getCoordinates(origV)[i]; cerr << (c).get_d() << " "; } cerr << endl;
        for(int i=0;i<3;i++) { VertCoord c = getCoordinates(v2V)[i]-getCoordinates(origV)[i]; cerr << (c).get_d() << " "; } cerr << endl;
      }
    }

    return ansSoS;
  } else {
    //assert( ((ans==1)&&(ansSoS)) || ((ans==-1)&&(!ansSoS)) ); //if SoS answer is true --> ans have to be 1
    //assert( ((ans==1)&&(ansOrig)) || ((ans==-1)&&(!ansOrig)) ); //if SoS answer is true --> ans have to be 1
    return ans==1;
  }

  
  
}

bool MeshIntersectionGeometry::isVertexInTriangleProjection(const Vertex &v1,const Vertex &v2, const Vertex &v3, const Vertex &queryPoint,int whatPlaneProjectTrianglesTo,TempVarsIsVertexTriangleProjection &tempVars) {
  int ans = isVertexInTriangleProjectionMainImpl(v1,v2,v3, queryPoint, whatPlaneProjectTrianglesTo,tempVars);

  #ifdef COLLECT_GEOMETRY_STATISTICS
    if(ans==0) {
      #pragma omp atomic
      geometryStatisticsDegenerateCases.ctDegeneraciesIsVertexInTriangleProjection++;
    } else {
      #pragma omp atomic
      geometryStatisticsNonDegenerateCases.ctDegeneraciesIsVertexInTriangleProjection++;
    }
  #endif

  if(ans==0) {
    bool ansSoS = isVertexInTriangleProjectionSoSImpl(v1,v2,v3, queryPoint, whatPlaneProjectTrianglesTo,tempVars);
    bool ansOrig = isVertexInTriangleProjectionOrig(v1,v2,v3, queryPoint, whatPlaneProjectTrianglesTo,tempVars);

    if(ansOrig!=ansSoS) {
      #pragma omp critical
      {
        cerr << "@@@ Is vertex in projection v1 v2 v3..." << endl;
        cerr << ansOrig << " " << ansSoS << " " << whatPlaneProjectTrianglesTo << endl;
        //for(int i=0;i<3;i++) { VertCoord c = getCoordinates(v1V)[i]-getCoordinates(origV)[i]; cerr << (c).get_d() << " "; } cerr << endl;
        //for(int i=0;i<3;i++) { VertCoord c = getCoordinates(v2V)[i]-getCoordinates(origV)[i]; cerr << (c).get_d() << " "; } cerr << endl;
      }
    }

    return ansSoS;
  } else {
    //assert( ((ans==1)&&(ansSoS)) || ((ans==-1)&&(!ansSoS)) ); //if SoS answer is true --> ans have to be 1
    //assert( ((ans==1)&&(ansOrig)) || ((ans==-1)&&(!ansOrig)) );
    return ans==1;
  }
}

bool MeshIntersectionGeometry::isVertexConvex(const Vertex &v1,const Vertex &queryVertex, const Vertex &v3,int whatPlaneProjectTrianglesTo,TempVarsIsVertexConvex &tempVars) {
  int ans = isVertexConvexMainImpl(v1,queryVertex,v3,whatPlaneProjectTrianglesTo,tempVars);
  

  #ifdef COLLECT_GEOMETRY_STATISTICS
    if(ans==0) {
      #pragma omp atomic
      geometryStatisticsDegenerateCases.ctDegeneraciesIsVertexConvex++;
    } else {
      #pragma omp atomic
      geometryStatisticsNonDegenerateCases.ctDegeneraciesIsVertexConvex++;
    }
  #endif

  if(ans==0) {
    //cerr << "Coincidency in isAngleGreater...: " << ans << " " << ansSoS << endl;
    bool ansSoS = isVertexConvexSoSImpl(v1,queryVertex,v3,whatPlaneProjectTrianglesTo,tempVars);
    bool ansOrig  = isVertexConvexOrig(v1,queryVertex,v3,whatPlaneProjectTrianglesTo,tempVars);

    assert(ansSoS==ansOrig);
    return ansSoS;
  } else {
    //assert( ((ans==1)&&(ansSoS)) || ((ans==-1)&&(!ansSoS)) ); //if SoS answer is true --> ans have to be 1
    //assert( ((ans==1)&&(ansOrig)) || ((ans==-1)&&(!ansOrig)) );
    return ans==1;
  }
}


/*************** PinMesh Part... ****************/
bool MeshIntersectionGeometry::isVertexInTriangleProjection(const InputTriangle &t, const InputVertex &queryPoint,TempVarsIsVertexTriangleProjectionZ0 &tempVars) const {
  int ans = isVertexInTriangleProjectionMainImpl(t, queryPoint,tempVars);
  

  #ifdef COLLECT_GEOMETRY_STATISTICS    
    if(ans==0) {
      #pragma omp atomic
      geometryStatisticsDegenerateCases.ctDegeneraciesIsVertexInInputTriangleProjection++;
    } else {
      #pragma omp atomic
      geometryStatisticsNonDegenerateCases.ctDegeneraciesIsVertexInInputTriangleProjection++;
    }
  #endif

  if(ans==0) {
    bool ansSoS = isVertexInTriangleProjectionSoSImpl(t, queryPoint,tempVars);
    bool ansOrig = isVertexInTriangleProjectionOrig(t, queryPoint,tempVars);
    //cerr << "Coincidency in isAngleGreater...: " << ans << " " << ansSoS << endl;
    if(ansOrig!=ansSoS) {
      #pragma omp critical
      {
        cerr << "@@@ Vertex in triangle projection z 0..." << endl;
        cerr << ansOrig << " " << ansSoS << " " <<  endl;
        //for(int i=0;i<3;i++) { VertCoord c = getCoordinates(v1V)[i]-getCoordinates(origV)[i]; cerr << (c).get_d() << " "; } cerr << endl;
        //for(int i=0;i<3;i++) { VertCoord c = getCoordinates(v2V)[i]-getCoordinates(origV)[i]; cerr << (c).get_d() << " "; } cerr << endl;
      }
    }

    return ansSoS;
  } else {
    //assert( ((ans==1)&&(ansSoS)) || ((ans==-1)&&(!ansSoS)) ); //if SoS answer is true --> ans have to be 1
    //assert( ((ans==1)&&(ansOrig)) || ((ans==-1)&&(!ansOrig)) );
    return ans==1;
  }
}

//TODO: use pre-computed normals here...
//If all points from the same mesh are translated equally, we will not need SoS here
//because we will never chose an input triangle that is vertical 
bool MeshIntersectionGeometry::isTriangleNormalPointingPositiveZ(const InputTriangle &t, TempVarIsTriangleNormalPointingPositiveZ &tempVars) const {
  int ans = isTriangleNormalPointingPositiveZMainImpl(t, tempVars);
  
  //bool ansOrig = isTriangleNormalPointingPositiveZOrig(t, tempVars);

  #ifdef COLLECT_GEOMETRY_STATISTICS
    if(ans==0) {
      #pragma omp atomic
      geometryStatisticsDegenerateCases.ctDegeneraciesIsTriangleNormalPointingPositiveZ++;
    } else {
      #pragma omp atomic
      geometryStatisticsNonDegenerateCases.ctDegeneraciesIsTriangleNormalPointingPositiveZ++;
    }
  #endif

  if(ans==0) {
    bool ansSoS = isTriangleNormalPointingPositiveZSoSImpl(t, tempVars);
    cerr << "@@@@ Coincidency in isTriangleNormalPoitingPositiveZ...: " << ans << " " << ansSoS << endl;
    return ansSoS;
  } else {
    //test disabled because SoS not implemented yet...
    //assert( ((ans==1)&&(ansSoS)) || ((ans==-1)&&(!ansSoS)) ); //if SoS answer is true --> ans have to be 1
    return ans==1;
  }
}


//These functions are only called when we have a coincidence....
//The only option is to call the SoS implementation...

//Given a vertex p, is p below the triangle t ? (we know p projected to z=0 is on t projected to z=0...)
bool MeshIntersectionGeometry::isTriangleAbovePointSoS(const InputTriangle &t, const InputVertex &p,TempVarIsTriangleAbovePointSoS &tempVars) const {
  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.ctDegeneraciesIsTriangleAbovePointSoS++;
  #endif

  return isTriangleAbovePointSoSOrig(t, p,tempVars);
}

//Given two triangles above a point, where the height above point is equal for both triangles, decide which one is lower according after SoS
//Make sure the two triangles are not the same... (this could happen, but we avoid that by checkin in PinMesh... )
const InputTriangle * MeshIntersectionGeometry::getBestTrianglePointInObjectSoS(const InputTriangle *candidateTriangle,const InputTriangle *bestTriangle, const InputVertex &p,TempVarGetBestTrianglePointInObjectSoS &tempVars) const {
  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.ctDegeneraciesGetBestTrianglePointInObjectSoS++;
  #endif

  return getBestTrianglePointInObjectSoSOrig(candidateTriangle,bestTriangle, p,tempVars);
}




//These are not predicates....


//SoS is easy here: just round (perturbation is always positive...)
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






/**********************************************************************************************************************************/
/**********************************************************************************************************************************/
/**********************************************************************************************************************************/
/**********************************************************************************************************************************/
/**********************************************************************************************************************************/
/**********************************************************************************************************************************/
/**********************************************************************************************************************************/
/**********************************************************************************************************************************/
/**********************************************************************************************************************************/
/**********************************************************************************************************************************/


/*
Maybe SoS? --> TODO: yes, for vertical triangles!

*/

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

//Input triangles should not be degenerated
//The projection of the triangle should not be degenerated (we always choose a plane whose projection is non degenerated)
//Thus, the sign will be either negative or positive (with or without SoS)
//SoS (infinitesimals) should not change the orientation of the triangles
//TODO: reuse normals --> ans: no need for now... only 0.11% of time!
bool MeshIntersectionGeometry::isTriangleClockwisedOriented(const InputTriangle &t,const int whatPlaneProjectTo, TempVarsIsTriangleClockwisedOriented &tempVars) const {
  return signCrossProduct2D(getCoordinates(*t.getInputVertex(0)),getCoordinates(*t.getInputVertex(1)),getCoordinates(*t.getInputVertex(0)),getCoordinates(*t.getInputVertex(2)),whatPlaneProjectTo,tempVars.tempCoords)<0;
}


//No need for SoS here (performance workaround...):
//Maybe we need SoS... what if the point is below a vertical triangle?
//Maybe we should have two versions of PinMesh: if point below vertical triangle --> use second version
//For consistency, we try to "hide" coordinates inside the HeightPointInTriangleProjection struct (so that PinMesh, like other classes, does not have access to coordinates)
//No need for SoS here: the triangle will be non-degenerate (input triangle) and, thus, the height of the projection
//onto the plane is well defined...
//TODO: what if triangle is vertical???? we still need to compute...
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


/*
Possible no need for SoS...

*/


array<VertCoord,3> MeshIntersectionGeometry::coordRangeMeshes() const {
  return {boundingBoxTwoMeshesTogetter[1][0]-boundingBoxTwoMeshesTogetter[0][0],boundingBoxTwoMeshesTogetter[1][1]-boundingBoxTwoMeshesTogetter[0][1],boundingBoxTwoMeshesTogetter[1][2]-boundingBoxTwoMeshesTogetter[0][2]};
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


//We are actually not using functions below...
/*
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
*/