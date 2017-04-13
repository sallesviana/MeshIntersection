
int MeshIntersectionGeometry::orientation(const InputVertex &v1, const InputVertex &v2, const InputVertex &queryPoint,int whatPlaneProjectTrianglesTo) const { 
  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.orientation2DOOO++;
  #endif

  const Point &p0 =  getCoordinates(v1);
  const Point &p1 =  getCoordinates(v2);
  const Point &p = getCoordinates(queryPoint);

  int coordY = 1;
  int coordX = 0;
  if(whatPlaneProjectTrianglesTo==PLANE_X0) { //if the triangle is projected to X=0 --> we need to use coordinates y,z (instead of x,y)
    coordX = 1;
    coordY = 2;
  } else if(whatPlaneProjectTrianglesTo ==PLANE_Y0) { //if the triangle is projected to Y=0 --> we need to use coordinates z,x (instead of x,y)
    coordX = 2;
    coordY = 0;
  }

  //a.x * b.y - a.y * b.x;
  return sgn( (p1[coordX]-p0[coordX])*(p[coordY]-p0[coordY]) -  (p1[coordY]-p0[coordY])*(p[coordX]-p0[coordX]) );
}

int MeshIntersectionGeometry::orientation(const InputVertex &v1,const InputVertex &v2, const VertexFromIntersection &queryPoint, int whatPlaneProjectTrianglesTo) const { 
  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.orientation2DOOI++;
  #endif

  const Point &p0 =  getCoordinates(v1);
  const Point &p1 =  getCoordinates(v2);
  const Point &p = getCoordinates(queryPoint);

  int coordY = 1;
  int coordX = 0;
  if(whatPlaneProjectTrianglesTo==PLANE_X0) { //if the triangle is projected to X=0 --> we need to use coordinates y,z (instead of x,y)
    coordX = 1;
    coordY = 2;
  } else if(whatPlaneProjectTrianglesTo ==PLANE_Y0) { //if the triangle is projected to Y=0 --> we need to use coordinates z,x (instead of x,y)
    coordX = 2;
    coordY = 0;
  }

  //a.x * b.y - a.y * b.x;
  return sgn( (p1[coordX]-p0[coordX])*(p[coordY]-p0[coordY]) -  (p1[coordY]-p0[coordY])*(p[coordX]-p0[coordX]) );
}

int MeshIntersectionGeometry::orientation(const InputVertex &v1, const VertexFromIntersection &v2, const VertexFromIntersection &queryPoint, int whatPlaneProjectTrianglesTo) const { 
  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.orientation2DOII++;
  #endif

  const Point &p0 =  getCoordinates(v1);
  const Point &p1 =  getCoordinates(v2);
  const Point &p = getCoordinates(queryPoint);

  int coordY = 1;
  int coordX = 0;
  if(whatPlaneProjectTrianglesTo==PLANE_X0) { //if the triangle is projected to X=0 --> we need to use coordinates y,z (instead of x,y)
    coordX = 1;
    coordY = 2;
  } else if(whatPlaneProjectTrianglesTo ==PLANE_Y0) { //if the triangle is projected to Y=0 --> we need to use coordinates z,x (instead of x,y)
    coordX = 2;
    coordY = 0;
  }

  //a.x * b.y - a.y * b.x;
  return sgn( (p1[coordX]-p0[coordX])*(p[coordY]-p0[coordY]) -  (p1[coordY]-p0[coordY])*(p[coordX]-p0[coordX]) );
}

int MeshIntersectionGeometry::orientation(const VertexFromIntersection &v1, const VertexFromIntersection &v2, const VertexFromIntersection &queryPoint, int whatPlaneProjectTrianglesTo) const { 
  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.orientation2DIII++;
  #endif


  const Point &p0 =  getCoordinates(v1);
  const Point &p1 =  getCoordinates(v2);
  const Point &p = getCoordinates(queryPoint);

  int coordY = 1;
  int coordX = 0;
  if(whatPlaneProjectTrianglesTo==PLANE_X0) { //if the triangle is projected to X=0 --> we need to use coordinates y,z (instead of x,y)
    coordX = 1;
    coordY = 2;
  } else if(whatPlaneProjectTrianglesTo ==PLANE_Y0) { //if the triangle is projected to Y=0 --> we need to use coordinates z,x (instead of x,y)
    coordX = 2;
    coordY = 0;
  }

  //a.x * b.y - a.y * b.x;
  return sgn( (p1[coordX]-p0[coordX])*(p[coordY]-p0[coordY]) -  (p1[coordY]-p0[coordY])*(p[coordX]-p0[coordX]) );
}

int cts[20] = {0};
//TODO: remove cts++ for performance...

//1 if make a left turn, -1 if make a right turn, 0 --> degeneracy (shouldn't happen...)
int MeshIntersectionGeometry::orientation(const Vertex &v1, const Vertex &v2, const Vertex &p, int whatPlaneProjectTrianglesTo) const { 
  bool isV1InputVertex = (&v1)->isInputVertex();
  bool isV2InputVertex = (&v2)->isInputVertex();  
  bool isPInputVertex = (&p)->isInputVertex();

  int numInputVertices = isV1InputVertex+isV2InputVertex+isPInputVertex;
  if(numInputVertices==3) {
    cts[0]++;
    return orientation(*static_cast<const InputVertex*>(&v1),*static_cast<const InputVertex*>(&v2),*static_cast<const InputVertex*>(&p),whatPlaneProjectTrianglesTo);
  } else if(numInputVertices==0) {
    cts[1]++;
    return orientation(*static_cast<const VertexFromIntersection*>(&v1),*static_cast<const VertexFromIntersection*>(&v2),*static_cast<const VertexFromIntersection*>(&p),whatPlaneProjectTrianglesTo);
  } else if(numInputVertices==1) {
    if(isV1InputVertex) {
      cts[2]++;
      return orientation(*static_cast<const InputVertex*>(&v1),*static_cast<const VertexFromIntersection*>(&v2),*static_cast<const VertexFromIntersection*>(&p),whatPlaneProjectTrianglesTo);
    } else if(isV2InputVertex) {
      cts[3]++;
      return -orientation(*static_cast<const InputVertex*>(&v2),*static_cast<const VertexFromIntersection*>(&v1),*static_cast<const VertexFromIntersection*>(&p),whatPlaneProjectTrianglesTo);
    } else { //p is the input vertex...
      cts[4]++;
      return -orientation(*static_cast<const InputVertex*>(&p),*static_cast<const VertexFromIntersection*>(&v2),*static_cast<const VertexFromIntersection*>(&v1),whatPlaneProjectTrianglesTo);
    }
  } else { //numInputVertices is 2...
    if(!isV1InputVertex) { // v2 and p are input vertices
      cts[5]++;
      return orientation(*static_cast<const InputVertex*>(&v2),*static_cast<const InputVertex*>(&p),*static_cast<const VertexFromIntersection*>(&v1),whatPlaneProjectTrianglesTo);
    } else if(!isV2InputVertex) { // v1 and p are input vertices
      cts[6]++;
      return orientation(*static_cast<const InputVertex*>(&p),*static_cast<const InputVertex*>(&v1),*static_cast<const VertexFromIntersection*>(&v2),whatPlaneProjectTrianglesTo);
    } else { // v1 and v2 are input vertices...
      cts[7]++;
      return orientation(*static_cast<const InputVertex*>(&v1),*static_cast<const InputVertex*>(&v2),*static_cast<const VertexFromIntersection*>(&p),whatPlaneProjectTrianglesTo);
    }
  }
}

/*
bool MeshIntersectionGeometry::isVertexInTriangleProjectionSoS(const Vertex &v1,const Vertex &v2, const Vertex &v3, const Vertex &queryPoint,int whatPlaneProjectTrianglesTo) {
  #pragma omp critical
  for(int i=0;i<8;i++) cerr << i << " --> " << cts[i] << endl;

  int o1 = orientation(v1,v2,queryPoint,whatPlaneProjectTrianglesTo);
  int o2 = orientation(v2,v3,queryPoint,whatPlaneProjectTrianglesTo);
  if(o1!=o2) return false;
  int o3 = orientation(v3,v1,queryPoint,whatPlaneProjectTrianglesTo);
  return o2==o3;
}*/


//SOS ORIENTATION of vertex with respect to a triangle and other basic geometric operations...

int signDeterminant4(const Point &p1,const Point &p2,const Point &p3,const Point &p) {
  const VertCoord &p1x = p1[0];
  const VertCoord &p1y = p1[1];
  const VertCoord &p1z = p1[2];

  const VertCoord &p2x = p2[0];
  const VertCoord &p2y = p2[1];
  const VertCoord &p2z = p2[2];

  const VertCoord &p3x = p3[0];
  const VertCoord &p3y = p3[1];
  const VertCoord &p3z = p3[2];

  const VertCoord &px = p[0];
  const VertCoord &py = p[1];
  const VertCoord &pz = p[2];

  const VertCoord det = -p1z*p2y*p3x + p1y*p2z*p3x + p1z*p2x*p3y - p1x*p2z*p3y - p1y*p2x*p3z +
                         p1x*p2y*p3z + p1z*p2y*px -  p1y*p2z*px -  p1z*p3y*px +  p2z*p3y*px + 
                         p1y*p3z*px -  p2y*p3z*px -  p1z*p2x*py +  p1x*p2z*py +  p1z*p3x*py - 
                         p2z*p3x*py -  p1x*p3z*py +  p2x*p3z*py +  p1y*p2x*pz -  p1x*p2y*pz - 
                         p1y*p3x*pz +  p2y*p3x*pz +  p1x*p3y*pz -  p2x*p3y*pz;
  return sgn(det);
}

int MeshIntersectionGeometry::orientation(const InputVertex&p1, const InputVertex&p2,const InputVertex&p3, const InputVertex &v) const {
  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.orientation3DOO++;
  #endif

  return signDeterminant4(getCoordinates(p1),getCoordinates(p2),getCoordinates(p3),getCoordinates(v));
}

int MeshIntersectionGeometry::orientation(const InputVertex&p1, const InputVertex&p2,const InputVertex&p3, const VertexFromIntersection &v) const {
  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.orientation3DOI++;
  #endif
  return signDeterminant4(getCoordinates(p1),getCoordinates(p2),getCoordinates(p3),getCoordinates(v));
}

int MeshIntersectionGeometry::orientation(const InputTriangle&t, const InputVertex &v) const {
  return orientation(*(t.getInputVertex(0)),*(t.getInputVertex(1)),*(t.getInputVertex(2)),v);
}

int MeshIntersectionGeometry::orientation(const InputTriangle&t, const VertexFromIntersection &v) const {
  return orientation(*(t.getInputVertex(0)),*(t.getInputVertex(1)),*(t.getInputVertex(2)),v);
}  


//this is basically what we need for the "brute force" retesselation...
//this is essentially a 1D orientation!!!!
//TODO: implement this as 1D orientation...
//what is the signal of each coordinate the vector from orig to dest
//cannot be 0 (SoS)
int MeshIntersectionGeometry::signalVectorCoord(const Vertex &orig, const Vertex &dest, int coord) const {
  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.signVector++;
  #endif

  const Point &p0 =  getCoordinates(orig);
  const Point &p1 =  getCoordinates(dest);
  int ans = sgn(p1[coord]-p0[coord]);
  return  (ans==0)?1:ans;
}


/*
Predicates and functions using SoS
*/



//Is v1 closer to origV than v2 is?
bool MeshIntersectionGeometry::isCloserSoSImpl(const InputVertex &orig, const VertexFromIntersection &v1, const VertexFromIntersection &v2, TempVarsIsCloser &tempVars) const {
  //we know that v1 and v2 are formed by the intersection of v1.triangle and an edge with endpoint in orig ; and v2.triangle and an
  //edge with endpoint in orig. 
  //Thus, v1 will be closer to orig iff v1 and orig are on the same side of v2.triangle 
  return orientation(v2.triangle,orig) == orientation(v2.triangle,v1);
}


bool MeshIntersectionGeometry::isVertexInTriangleProjectionSoSImpl(const Vertex &v1,const Vertex &v2, const Vertex &v3, const Vertex &queryPoint,int whatPlaneProjectTrianglesTo,TempVarsIsVertexTriangleProjection &tempVars) const {
  int o1 = orientation(v1,v2,queryPoint,whatPlaneProjectTrianglesTo);
  int o2 = orientation(v2,v3,queryPoint,whatPlaneProjectTrianglesTo);
  int o3 = orientation(v3,v1,queryPoint,whatPlaneProjectTrianglesTo);
  //this will not work w/o SoS if point exactly on boundary...

  //TODO: change code after implementing SoS properly...
  /*
  if(o1!=o2) return false;
  
  return o2==o3;*/

  if(o2!=0 && o3!=0 && o2!=o3)
    return false;
  if(o1!=0 && o2!=0 && o1!=o2)
    return false;
  if(o1!=0 && o3!=0 && o1!=o3)
    return false;
  return true;
}


bool MeshIntersectionGeometry::isVertexConvexSoSImpl(const Vertex &v1,const Vertex &queryVertex, const Vertex &v3,int whatPlaneProjectTrianglesTo,TempVarsIsVertexConvex &tempVars) const {
  return orientation(v1,queryVertex,v3,whatPlaneProjectTrianglesTo)<0;
}


bool MeshIntersectionGeometry::isVertexInTriangleProjectionSoSImpl(const InputTriangle &t, const InputVertex &queryPoint,TempVarsIsVertexTriangleProjectionZ0 &tempVars) const {
  const Vertex &v1 = *(t.getInputVertex(0));
  const Vertex &v2 = *(t.getInputVertex(1));
  const Vertex &v3 = *(t.getInputVertex(2));

  const int whatPlaneProjectTrianglesTo = PLANE_Z0;

  //TODO: after implementing SoS properly... fix
  /*int o1 = orientation(v1,v2,queryPoint,whatPlaneProjectTrianglesTo);
  int o2 = orientation(v2,v3,queryPoint,whatPlaneProjectTrianglesTo);
  if(o1!=o2) return false;
  int o3 = orientation(v3,v1,queryPoint,whatPlaneProjectTrianglesTo);
  return o2==o3;
  */

  int o1 = orientation(v1,v2,queryPoint,whatPlaneProjectTrianglesTo);
  int o2 = orientation(v2,v3,queryPoint,whatPlaneProjectTrianglesTo);
  int o3 = orientation(v3,v1,queryPoint,whatPlaneProjectTrianglesTo);
  if(o2!=0 && o3!=0 && o2!=o3)
    return false;
  if(o1!=0 && o2!=0 && o1!=o2)
    return false;
  if(o1!=0 && o3!=0 && o1!=o3)
    return false;
  return true;
}


//Given a vertex p, is p below the triangle t ? (we know p projected to z=0 is on t projected to z=0...)
bool MeshIntersectionGeometry::isTriangleAbovePointSoSImpl(const InputTriangle &t, const InputVertex &p,TempVarIsTriangleAbovePointSoS &tempVars) const {
  int sideOfTriangle = orientation(t,p); //is p on the positive (1) or negative (-1) side of the triangle?

  TempVarIsTriangleNormalPointingPositiveZ temp;

  //if the point is on the positive side, it will be below the triangle if it points down
  if(sideOfTriangle==1) return !isTriangleNormalPointingPositiveZSoSImpl(t,temp);
  else return isTriangleNormalPointingPositiveZSoSImpl(t,temp);
  //if the point is on the negative sie, it will be below the triangle if it points up...
}



bool MeshIntersectionGeometry::isOrientationPositiveSoSImpl(const Vertex &origV, const Vertex &v1V, const Vertex &v2V, const int planeToProject, TempVarsIsAngleWith0Greater &tempVars) const {
  return orientation(origV,v1V,v2V,planeToProject)>0;
}

//we need to consider point on axis and the angle...
bool MeshIntersectionGeometry::isAngleWith0GreaterSoSImpl(const Vertex &origV, const Vertex &v1V, const Vertex &v2V, const int planeToProject, TempVarsIsAngleWith0Greater &tempVars) const {
  //TODO: check on what side of the "x" axis the vectors are (origV,v1V), (origV,v2V)...

  int xCoord = 0;
  int yCoord = 1;
  if(planeToProject==PLANE_Y0) {
    xCoord= 2;
    yCoord= 0;
  } else if(planeToProject==PLANE_X0) {
    xCoord= 1;
    yCoord= 2;    
  } 

  const int sgnV1y = signalVectorCoord(origV,v1V,yCoord);
  const int sgnV2y = signalVectorCoord(origV,v2V,yCoord);

  if(sgnV1y>0 && sgnV2y<0) return true; //check if the two vectors are in different sides of the x axis...
  if(sgnV1y<0 && sgnV2y>0) return false;

  //they are in the same side of the x axis...
  const int sgnV1x = signalVectorCoord(origV,v1V,xCoord);
  const int sgnV2x = signalVectorCoord(origV,v2V,xCoord);

  if(sgnV1y>0) { //are both at the positive y?
    //check if their x is different...
    if(sgnV1x > 0 && sgnV2x < 0) return true;
    if(sgnV1x < 0 && sgnV2x > 0) return false;
  } else { //are they both at the negative y?
    if(sgnV1x > 0 && sgnV2x < 0) return false;
    if(sgnV1x < 0 && sgnV2x > 0) return true;
  }

  
  //same side of the x-axis --> use vector orientation...
  return orientation(origV,v1V,v2V,planeToProject)>0;
}

//TODO: use pre-computed normals here...
//If all points from the same mesh are translated equally, we will not need SoS here
//because we will never chose an input triangle that is vertical 

//TO think: is it true that if we project the triangle to z=0, this function
//should return true if v3 is to the left of the vector v1-v2 ?
bool MeshIntersectionGeometry::isTriangleNormalPointingPositiveZSoSImpl(const InputTriangle &t, TempVarIsTriangleNormalPointingPositiveZ &tempVars) const {
    return orientation(*t.getInputVertex(0),*t.getInputVertex(1),*t.getInputVertex(2),PLANE_Z0)==1;
}




bool MeshIntersectionGeometry::intersectEdgeWithTriangleSoSImpl(const InputTriangle &triangle, const InputVertex &p1, const InputVertex &p2, TempVarsComputeIntersections &tempVars) const {
  int orientationP1Triangle = orientation(triangle,p1);
  int orientationP2Triangle = orientation(triangle,p2);
  if(orientationP1Triangle==orientationP2Triangle) return false; //both are on the same side of the triangle's plane...

  const InputVertex &a = *(triangle.getInputVertex(0));
  const InputVertex &b = *(triangle.getInputVertex(1));
  const InputVertex &c = *(triangle.getInputVertex(2));

  if(orientationP1Triangle>0) return (orientation(a,b,p2,p1)>0) && (orientation(b,c,p2,p1)>0) && (orientation(c,a,p2,p1)>0);
  else return (orientation(a,b,p1,p2)>0) && (orientation(b,c,p1,p2)>0) && (orientation(c,a,p1,p2)>0);
}



//If intersection happens at boundaries for example --> we really need SoS to be reliable...
//ex: maybe w/o SoS we would have edge-triangle, but with SoS the intersection is triangle-edge...
//think a little more to make sure we really cannot compute the coordinates..
bool MeshIntersectionGeometry::intersectTwoTrianglesSoSImpl(const InputTriangle &triMesh0,const InputTriangle &triMesh1,
             Point &coordsPt1,VertexFromIntersection &vertexThatCreatedPt1, Point &coordsPt2,
             VertexFromIntersection &vertexThatCreatedPt2, TempVarsComputeIntersections &tempVars) const {
  //because of SoS, we will have exactly two edge-triangle intersections
  //we need to test all 6 possibilities of edge-triangle intersections...

  int numIntersectionsFound = 0;

  for(int i=0;i<3;i++) {
    const InputVertex &v0 = *(triMesh1.getInputVertex(i));
    const InputVertex &v1 = *(triMesh1.getInputVertex((i+1)%3));
    if(intersectEdgeWithTriangleSoSImpl(triMesh0 ,v0,v1,tempVars)) {
      if(numIntersectionsFound==0) {
        vertexThatCreatedPt1.triangle = triMesh0;
        vertexThatCreatedPt1.edge[0] = v0;
        vertexThatCreatedPt1.edge[1] = v1;
      } else {
        vertexThatCreatedPt2.triangle = triMesh0;
        vertexThatCreatedPt2.edge[0] = v0;
        vertexThatCreatedPt2.edge[1] = v1;
      }
      numIntersectionsFound++;
    }
  }
  for(int i=0;i<3;i++) {
    const InputVertex &v0 = *(triMesh0.getInputVertex(i));
    const InputVertex &v1 = *(triMesh0.getInputVertex((i+1)%3));
    if(intersectEdgeWithTriangleSoSImpl(triMesh1 ,v0,v1,tempVars)) {
      if(numIntersectionsFound==0) {
        vertexThatCreatedPt1.triangle = triMesh1;
        vertexThatCreatedPt1.edge[0] = v0;
        vertexThatCreatedPt1.edge[1] = v1;
      } else {
        vertexThatCreatedPt2.triangle = triMesh1;
        vertexThatCreatedPt2.edge[0] = v0;
        vertexThatCreatedPt2.edge[1] = v1;
      }
      numIntersectionsFound++;
    }
  }



  assert(numIntersectionsFound==0 || numIntersectionsFound==2); //with SoS, we can have only either 2 or 0 intersections between edges of a triangle and another triangle...
  return numIntersectionsFound>0;
}






//Given two triangles above a point, where the height above point is equal for both triangles, decide which one is lower according after SoS
//Make sure the two triangles are not the same... (this could happen, but we avoid that by checkin in PinMesh... )
const InputTriangle * MeshIntersectionGeometry::getBestTrianglePointInObjectSoSImpl(const InputTriangle *candidateTriangle,const InputTriangle *bestTriangle, const InputVertex &p,TempVarGetBestTrianglePointInObjectSoS &tempVars) const {
  cerr << "TODO: SoS get best triangle point in object SoS" << endl;
  return candidateTriangle;
}

