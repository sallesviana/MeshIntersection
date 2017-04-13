
/*
Predicates and functions without SoS
Returns 0 if we have a coincidency...
*/

#include "tritri_isectline.c"


//Is v1 closer to origV than v2 is?
int MeshIntersectionGeometry::isCloserMainImpl(const InputVertex &origV, const VertexFromIntersection &v1V, const VertexFromIntersection &v2V, TempVarsIsCloser &tempVars) const {
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

  if(tempVertCoords[2]<tempVertCoords[3]) return 1;
  else if(tempVertCoords[2]>tempVertCoords[3]) return -1;
  else return 0;
}

//will return 0 if vector is on an axis or if they are coincident...
int MeshIntersectionGeometry::isAngleWith0GreaterMainImpl(const Vertex &origV, const Vertex &v1V, const Vertex &v2V, const int planeToProject, TempVarsIsAngleWith0Greater &tempVars) const {
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

  if(sgnV1y>0 && sgnV2y<0) return 1; //check if the two vectors are in different sides of the x axis...
  if(sgnV1y<0 && sgnV2y>0) return -1;

  /*
  const int sgnV1x = sgn(v1x);
  const int sgnV2x = sgn(v2x);  

  if(sgnV1x>=0 && sgnV1y==0) return 0; //vector on y=0 , x = 0+
  if(sgnV2x>=0 && sgnV2y==0) return 0; //vector on y=0 , x = 0+

  //are the points on different quadrants?
  if(sgnV1y!=0 && sgnV2y !=0 && sgnV1x!=0 && sgnV2x!=0) { //no vertex on an axis... trivial...

  }
  
  TODO:consider quadrant if slow......
  */

  
  if(sgnV1y==0 || sgnV2y==0) return 0; //coincidency... the point is on an axis..

  //TODO: maybe accelerate by checking if in different quadrants...

  //is the cross product positive?
  v1x *= v2y; //const VertCoord component1 = v1x*v2y; //v1.x * v2.y 
  v2x *= v1y; //const VertCoord component2 = v2x*v1y; //v2.x*v1.y

  //the cross product is = component1-component2
  //is it positive?
  if(v1x > v2x) return 1;//return (component1 > component2);
  else if(v1x < v2x) return -1;
  else return 0;
}


//this works only when no angle is 0 and no edge is degenerate
int MeshIntersectionGeometry::isAngleWith0GreaterNonZeroAngleMainImpl(const Vertex &origV, const Vertex &v1V, const Vertex &v2V, const int planeToProject, TempVarsIsAngleWith0Greater &tempVars) const {
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

  if(sgnV1y>0 && sgnV2y<0) return 1; //check if the two vectors are in different sides of the x axis...
  if(sgnV1y<0 && sgnV2y>0) return -1;


  //TODO: maybe accelerate by checking if in different quadrants...

  //is the cross product positive?
  v1x *= v2y; //const VertCoord component1 = v1x*v2y; //v1.x * v2.y 
  v2x *= v1y; //const VertCoord component2 = v2x*v1y; //v2.x*v1.y

  //the cross product is = component1-component2
  //is it positive?
  if(v1x > v2x) return 1;//return (component1 > component2);
  else if(v1x < v2x) return -1;
  else return 0;
}


int MeshIntersectionGeometry::isVertexInTriangleProjectionMainImpl(const Vertex &v1,const Vertex &v2, const Vertex &v3, const Vertex &queryPoint,int whatPlaneProjectTrianglesTo,TempVarsIsVertexTriangleProjection &tempVars) const {
  const Point &p = getCoordinates(queryPoint);
  const Point &p0 =  getCoordinates(v1);
  const Point &p1 =  getCoordinates(v2);
  const Point &p2 =  getCoordinates(v3);

  if ( p==p0 || p==p1 || p==p2) {
    return 0; //is the point directly on a vertex of the triangle?
  }

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
  if (denominator==0) { 
    return -1; //because of SoS, vertices will never be below vertical triangles...
  }
  VertCoord a = ((p1[coordY] - p2[coordY])*(p[coordX] - p2[coordX]) + (p2[coordX] - p1[coordX])*(p[coordY] - p2[coordY])) / denominator;
  if ( a<0 || a >1) return -1;
  
  VertCoord b = ((p2[coordY] - p0[coordY])*(p[coordX] - p2[coordX]) + (p0[coordX] - p2[coordX])*(p[coordY] - p2[coordY])) / denominator;

  if (b<0 || b>1) return -1;
  VertCoord c = 1 - a - b;
 
  //if ( (fabs(a) <= EPS && (fabs(b) <= EPS) || fabs(c) <= EPS)  || (fabs(b) <= EPS && fabs(c) <= EPS) ) return false; // the point is one of the 3 triangle vertices...
  //if ( (fabs(a) <= EPS && fabs(b) <= EPS) || (fabs(a) <= EPS && fabs(c) <= EPS)  || (fabs(b) <= EPS && fabs(c) <= EPS) ) return false; // the point is one of the 3 triangle vertices...
  //return 0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1; //maybe we should perform some tests using EPSILON to avoid floating point errors...
  
  if(a==0 || a==1 || b==0 || b==1 || c==0 || c==1) return 0;
  else if(0 < c && c < 1) return 1;
  else return -1; 
}


int MeshIntersectionGeometry::isVertexConvexMainImpl(const Vertex &v1,const Vertex &queryVertex, const Vertex &v3,int whatPlaneProjectTrianglesTo,TempVarsIsVertexConvex &tempVars) const {
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
  

  return -sgn(tempCoords[0]); //if the sign is negative --> the vertex is convex...
}


/*************** PinMesh Part... ****************/
int MeshIntersectionGeometry::isVertexInTriangleProjectionMainImpl(const InputTriangle &t, const InputVertex &queryPoint,TempVarsIsVertexTriangleProjectionZ0 &tempVars) const {
  const Point &p0 = getCoordinates(*t.getInputVertex(0));
  const Point &p1 = getCoordinates(*t.getInputVertex(1));
  const Point &p2 = getCoordinates(*t.getInputVertex(2));

  const Point &p = getCoordinates(queryPoint);

  if ( p==p0 || p==p1 || p==p2)  {
    return 0;
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
    //Because of SoS no point should be below a vertical triangle....
    return 0; //pointInTriangleProjSoS(p0, p1, p2, p) ; //SoS
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

  if (sgn(tempVertCoords[1])<0 || tempVertCoords[1]>1) return -1; //if ( a<0 || a >1) return false;
  
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

  if (sgn(tempVertCoords[2]) <0 || tempVertCoords[2]>1) return -1;//if (b<0 || b>1) return false;
  //VertCoord c = 1 - a - b;
  tempVertCoords[3] = 1;
  tempVertCoords[3] -= tempVertCoords[1];
  tempVertCoords[3] -= tempVertCoords[2];


  if (sgn(tempVertCoords[3]) < 0 || tempVertCoords[3] >1) return -1;

  //All coordinates are >=0 && <=1... so we have either a degenerate case or the point is in the triangle...
  if (sgn(tempVertCoords[3])==0 || sgn(tempVertCoords[2])==0 || sgn(tempVertCoords[1])==0) {
    return 0;
    //bool ans =  pointInTriangleProjSoS(p0, p1, p2, p) ; //SoS
    //return ans;
  }
  return 1;
}



//TODO: use pre-computed normals here...
//If all points from the same mesh are translated equally, we will not need SoS here
//because we will never chose an input triangle that is vertical 
int MeshIntersectionGeometry::isTriangleNormalPointingPositiveZMainImpl(const InputTriangle &t, TempVarIsTriangleNormalPointingPositiveZ &tempVars) const {
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

    /*if (C==0) {
      return 0;
    }
    return C>0; */
    return sgn(C);
}



