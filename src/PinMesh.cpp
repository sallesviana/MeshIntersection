#include "utils.h"
#include "3dGeometry.h"
#include "nested3DGrid.h"



#include "PinMesh.h"

const int LABEL_NON_EMPTY_CELL = -3;

/*

//Returns true iff the projection of the point p into a horizontal plane is inside of the projection of triangle
// (p0,p1,p2) into the horizontal plane...

//returns true if the point p is vertically directly above/below the point p...
//TODO special cases:
//- degenerate triangle?
//- point above/below vertex/edge?
//- vertical triangle?
bool PinMesh::pointInTriangleProj(const Point &p0, const Point &p1, const Point &p2, const Point &p)  {
  if ( p==p0 || p==p1 || p==p2) return false; //is the point directly above a vertex of the triangle?
  
  VertCoord denominator = ((p1[1] - p2[1])*(p0[0] - p2[0]) + (p2[0] - p1[0])*(p0[1] - p2[1]));
  if (denominator==0) { //TODO: check this.... degenerate triangles or vertical triangles (a segment never intersects a vertical triangle...)
    return false;
  }
  VertCoord a = ((p1[1] - p2[1])*(p[0] - p2[0]) + (p2[0] - p1[0])*(p[1] - p2[1])) / denominator;
  if ( a<0 || a >1) return false;
  
  VertCoord b = ((p2[1] - p0[1])*(p[0] - p2[0]) + (p0[0] - p2[0])*(p[1] - p2[1])) / denominator;

  if (b<0 || b>1) return false;
  VertCoord c = 1 - a - b;
 
  //if ( (fabs(a) <= EPS && (fabs(b) <= EPS) || fabs(c) <= EPS)  || (fabs(b) <= EPS && fabs(c) <= EPS) ) return false; // the point is one of the 3 triangle vertices...
  //if ( (fabs(a) <= EPS && fabs(b) <= EPS) || (fabs(a) <= EPS && fabs(c) <= EPS)  || (fabs(b) <= EPS && fabs(c) <= EPS) ) return false; // the point is one of the 3 triangle vertices...
  //return 0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1; //maybe we should perform some tests using EPSILON to avoid floating point errors...
  
  
  return 0 <= c && c <= 1; //maybe we should perform some tests using EPSILON to avoid floating point errors...
}




//Returns true iff the projection of the point p into a horizontal plane is inside of the projection of triangle
// (p0,p1,p2) into the horizontal plane...

//returns true if the point p is vertically directly above/below the point p...

// Uses SoS to treat the degenerate cases
// We pretend the point (p) is sligtly translated translated: translation (epsilon, epsilon^2, epsilon^3) 

//TODO special cases:
//- degenerate triangle?
//- point above/below vertex/edge?
//- vertical triangle?

// a_epsilon = a + (epsilon(p1y-p2y) + epsilon^2 (p2x-p1x) )/den
// b_epsilon = b + (epsilon(p2y-p0y) + epsilon^2 (p0x-p2x) )/den
// c_epsilon = 1 - a_epsilon - b_epsilon = 1 - a - (epsilon(p1y-p2y) + epsilon^2 (p2x-p1x) )/den - b - (epsilon(p2y-p0y) + epsilon^2 (p0x-p2x) )/den
// c_epsilon = 1 - a - b - ( epsilon(p1y-p0y) + epsilon^2(p0x-p1x)  )/den
bool PinMesh::pointInTriangleProjSoS(const Point &p0, const Point &p1, const Point &p2, const Point &p)  { 
  assert(p0!=p1 && p0!=p2 && p1!=p2);

  VertCoord denominator = ((p1[1] - p2[1])*(p0[0] - p2[0]) + (p2[0] - p1[0])*(p0[1] - p2[1]));
  if (sgn(denominator)==0) { //TODO: check this.... degenerate triangles or vertical triangles (a segment never intersects a vertical triangle...)
    return false;
  }
  VertCoord a = ((p1[1] - p2[1])*(p[0] - p2[0]) + (p2[0] - p1[0])*(p[1] - p2[1])) / denominator;
  VertCoord b = ((p2[1] - p0[1])*(p[0] - p2[0]) + (p0[0] - p2[0])*(p[1] - p2[1])) / denominator;
  VertCoord c = 1 - a - b;

  if (sgn(a)==0) {
    if (sgn(b)==0) { //a_epsilon,b_epsilon should be greater than 0 (thus, c_epsilon will be slightly smaller than 1...)
      bool aEpsilonGreater0 = (p1[1]>p2[1] || ( (p1[1]==p2[1]) && (p2[0]>p1[0]) )); //(p1[1]==p2[1]) && (p2[0]==p1[0]) will always be false! (otherwise the points p1 and p2 would be the same... what is invalid input..)
      if (sgn(denominator)<0) aEpsilonGreater0 = !aEpsilonGreater0;

      bool bEpsilonGreater0 = (p2[1]>p0[1] || ( (p2[1]==p0[1]) && (p0[0]>p2[0]) ));
      if (sgn(denominator)<0) bEpsilonGreater0 = !bEpsilonGreater0;

      return aEpsilonGreater0 && bEpsilonGreater0;
    }
    else if(b!=1) {  //a_epsilon should be slightly greater than 0....
      bool aEpsilonGreater0 = (p1[1]>p2[1] || ( (p1[1]==p2[1]) && (p2[0]>p1[0]) )); //(p1[1]==p2[1]) && (p2[0]==p1[0]) will always be false! (otherwise the points p1 and p2 would be the same... what is invalid input..)
      if (sgn(denominator)<0) aEpsilonGreater0 = !aEpsilonGreater0;
      return aEpsilonGreater0; //invert the return value if the denominator is negative (because, in this case, the sign of the expression changes...)
    } else { //b==1 --> b_epsilon should be slightly smaller than 1 and a_epsilon/c_epsilon should be slightly greater than 0 ...
      bool aEpsilonGreater0 = (p1[1]>p2[1] || ( (p1[1]==p2[1]) && (p2[0]>p1[0]) ));
      if (sgn(denominator)<0) aEpsilonGreater0 = !aEpsilonGreater0;


      bool bEpsilonSmaller1 = (p0[1]>p2[1]) || ( (p0[1]==p2[1]) && (p2[0]>p0[0]) );
      if (sgn(denominator)<0) bEpsilonSmaller1 = !bEpsilonSmaller1;

      // c_epsilon = 1 - a - b - ( epsilon(p1y-p0y) + epsilon^2(p0x-p1x)  )/den
      // (a+b = 1) --> c_epsilon = - ( epsilon(p1y-p0y) + epsilon^2(p0x-p1x)  )/den
      // --> c_epsilon > 0 <-> - ( epsilon(p1y-p0y) + epsilon^2(p0x-p1x)  )/den  > 0 <-> ( epsilon(p1y-p0y) + epsilon^2(p0x-p1x)  )/den  < 0
      bool cEpsilonGreater0 = (p0[1]>p1[1]) || ( (p0[1]==p1[1]) && (p1[0]>p0[0]) );
      if (sgn(denominator)<0) cEpsilonGreater0 = !cEpsilonGreater0;

      return aEpsilonGreater0 && bEpsilonSmaller1 && cEpsilonGreater0; 
    }
  } else if (sgn(b)==0) {
    if(a!=1) {  // 0<a<1 (otherwise the previous if statemente would be true...), 0<c<1
      bool bEpsilonGreater0 = (p2[1]>p0[1] || ( (p2[1]==p0[1]) && (p0[0]>p2[0]) ));
      if (sgn(denominator)<0) bEpsilonGreater0 = !bEpsilonGreater0;
      return bEpsilonGreater0; 
    } else { // a=1, b=0, c=0 --> aEpsilon should be slightly less than 1 and b epsilon slightly more than 0 (idem c_epsilon)...
      bool bEpsilonGreater0 = (p2[1]>p0[1] || ( (p2[1]==p0[1]) && (p0[0]>p2[0]) ));
      if (sgn(denominator)<0) bEpsilonGreater0 = !bEpsilonGreater0;


      bool aEpsilonSmaller1 = (p2[1]>p1[1]) || ( (p2[1]==p1[1]) && (p1[0]>p2[0]) );
      if (sgn(denominator)<0) aEpsilonSmaller1 = !aEpsilonSmaller1;

      bool cEpsilonGreater0 = (p0[1]>p1[1]) || ( (p0[1]==p1[1]) && (p1[0]>p0[0]) );
      if (sgn(denominator)<0) cEpsilonGreater0 = !cEpsilonGreater0;

      return bEpsilonGreater0 && aEpsilonSmaller1 && cEpsilonGreater0; 
    }
  } else { 
    assert(sgn(c)==0);
    //c == 0 . To be true, c_epsilon > 0 --> we need: a_epsilon+b_epsilon slightly less than 1.... (and a_epsilon>0, b_epsilon>0, what is already true because of the previous two if statements!)
    //So, we only need a_epsilon+b_epsilon slightly less than 1
    // a+b is = 1 --> we need the "epsilon" part of a_epsilon+b_epsilon to be < 0
    // That is, we need:  = (epsilon(p1y-p2y) + epsilon^2 (p2x-p1x) )/den + (epsilon(p2y-p0y) + epsilon^2 (p0x-p2x) )/den  < 0
    // ==  (epsilon*(p1y+p2y-p2y-p0y) + epsilon^2(p2x+p0x-p2x-p1x))/den <0
    // ==  (epsilon*(p1y-p0y) + epsilon^2(p0x-p1x))/den <0
    bool cGreater0 = ( (p0[1] > p1[1]) || ( (p0[1] == p1[1]) && (p1[0]>p0[0]) ) );
    if (sgn(denominator) <0) cGreater0 = !cGreater0;

    return cGreater0;
  }
}






//returns true if the point p is vertically directly above/below the point p...
//TODO special cases: --> DONE!
//- degenerate triangle?
//- point above/below vertex/edge?
//- vertical triangle?
// bool pointInTriangleProj(const Triangle &triangle,const Point &p) { 
//   VertCoord denominator = ((triangle[1][1] - triangle[2][1])*(triangle[0][0] - triangle[2][0]) + (triangle[2][0] - triangle[1][0])*(triangle[0][1] - triangle[2][1]));
//   VertCoord a = ((triangle[1][1] - triangle[2][1])*(p[0] - triangle[2][0]) + (triangle[2][0] - triangle[1][0])*(p[1] - triangle[2][1])) / denominator;
//   if ( a<0 || a >1) return false;
  
//   VertCoord b = ((triangle[2][1] - triangle[0][1])*(p[0] - triangle[2][0]) + (triangle[0][0] - triangle[2][0])*(p[1] - triangle[2][1])) / denominator;

//   if ( b<0 || b>1) return false;

//   VertCoord c = 1 - a - b;

//   return 0 <= c && c <= 1; 


// /*  point_coord denominator = ((p1.second - p2.second)*(p0.first - p2.first) + (p2.first - p1.first)*(p0.second - p2.second));
//   point_coord a = ((p1.second - p2.second)*(p.first - p2.first) + (p2.first - p1.first)*(p.second - p2.second)) / denominator;
//   if ( a<0 || a >1) return false;
  
//   point_coord b = ((p2.second - p0.second)*(p.first - p2.first) + (p0.first - p2.first)*(p.second - p2.second)) / denominator;

//   if ( b<0 || b>1) return false;

//   point_coord c = 1 - a - b;

//   return 0 <= c && c <= 1; 
// }

//p0,p1,p2 are the vertices of the triangle...
void PinMesh::getHeigthAbovePoint(VertCoord &heightAbovePoint,const Point &p0,const Point &p1,const Point &p2,const Point &p) {
  // lets compute the cross product of two vectors of the triangle in order to get the plane equation A(x-x0)+ B(y-y0)+C(z-z0)= 0
  //more info: http://tutorial.math.lamar.edu/Classes/CalcIII/EqnsOfPlanes.aspx
  VertCoord vec[2][3];
  vec[0][0] = p1[0]-p0[0];
  vec[0][1] = p1[1]-p0[1];
  vec[0][2] = p1[2]-p0[2];

  vec[1][0] = p2[0]-p0[0];
  vec[1][1] = p2[1]-p0[1];
  vec[1][2] = p2[2]-p0[2];

  //lets compute the values A,B,C basing on the two 3D vectors vec[0] and vec[1]
  //  | i             j           k      i             j       
  //  | vec[0][0]  vec[0][1]  vec[0][2]  vec[0][0]  vec[0][1]  
  //  | vec[1][0]  vec[1][1]  vec[1][2]  vec[1][0]  vec[1][1] 

  VertCoord A = vec[0][1]*vec[1][2] - vec[0][2]*vec[1][1];
  //VertCoord B = vec[0][0]*vec[1][2] - vec[0][2]*vec[1][0];
  VertCoord B = vec[0][2]*vec[1][0] - vec[0][0]*vec[1][2];
  VertCoord C = vec[0][0]*vec[1][1] - vec[0][1]*vec[1][0];

  assert(C!=0);

  // now we have an equation A(x-triangle[0][0]) +  B(y-triangle[0][1]) + C(z-triangle[0][2]) = 0
  // A(x-triangle[0][0]) +  B(y-triangle[0][1]) + C*z - C*triangle[0][2] = 0
  // C*z = C*triangle[0][2] - A(x-triangle[0][0]) - B(y-triangle[0][1])
  // C*z = A(triangle[0][0]-x) +  B(triangle[0][1]-y) + C*triangle[0][2] 

  heightAbovePoint = ( A*(p0[0]-p[0]) +  B*(p0[1]-p[1]) +  C*p0[2])/C; //TODO: check this...
}


//--------------------------------------------
//--------------------------------------------
//--------------------------------------------







//Returns true iff the projection of the point p into a horizontal plane is inside of the projection of triangle
// (p0,p1,p2) into the horizontal plane...
//Uses SoS to treat the special cases (degeneracies...)
//tempVertCoords shoudl have at least 5 elements (this variable is used to avoid allocating temporary memory).
bool PinMesh::pointInTriangleProj(const Point &p0, const Point &p1, const Point &p2, const Point &p, VertCoord *tempVertCoords)  {
  if ( p==p0 || p==p1 || p==p2)  {
    return pointInTriangleProjSoS(p0, p1, p2, p) ;//return 0; // degenerate case --> SoS
  }
  
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
    //cerr << "SoS: " << (sgn(tempVertCoords[3])==0)  << (sgn(tempVertCoords[2])==0) << (sgn(tempVertCoords[1])==0) << endl;
    bool ans =  pointInTriangleProjSoS(p0, p1, p2, p) ; //SoS
    //cerr << "ans: " << ans << endl;
    return ans;
  }
  return 1;
}

//returns true if the point p is vertically directly above/below the point p...
//TODO special cases:
//- degenerate triangle?
//- point above/below vertex/edge?
//- vertical triangle?
// bool pointInTriangleProj(const Triangle &triangle,const Point &p) { 
//   VertCoord denominator = ((triangle[1][1] - triangle[2][1])*(triangle[0][0] - triangle[2][0]) + (triangle[2][0] - triangle[1][0])*(triangle[0][1] - triangle[2][1]));
//   VertCoord a = ((triangle[1][1] - triangle[2][1])*(p[0] - triangle[2][0]) + (triangle[2][0] - triangle[1][0])*(p[1] - triangle[2][1])) / denominator;
//   if ( a<0 || a >1) return false;
  
//   VertCoord b = ((triangle[2][1] - triangle[0][1])*(p[0] - triangle[2][0]) + (triangle[0][0] - triangle[2][0])*(p[1] - triangle[2][1])) / denominator;

//   if ( b<0 || b>1) return false;

//   VertCoord c = 1 - a - b;

//   return 0 <= c && c <= 1; 


//   point_coord denominator = ((p1.second - p2.second)*(p0.first - p2.first) + (p2.first - p1.first)*(p0.second - p2.second));
//   point_coord a = ((p1.second - p2.second)*(p.first - p2.first) + (p2.first - p1.first)*(p.second - p2.second)) / denominator;
//   if ( a<0 || a >1) return false;
  
//   point_coord b = ((p2.second - p0.second)*(p.first - p2.first) + (p0.first - p2.first)*(p.second - p2.second)) / denominator;

//   if ( b<0 || b>1) return false;

//   point_coord c = 1 - a - b;

//   return 0 <= c && c <= 1; 
// }

//p0,p1,p2 are the vertices of the triangle...
// vec is a temporary matrix of coordinates
// tempVertCoords is a temporary array with size at least 4
void PinMesh::getHeigthAbovePoint(VertCoord &heightAbovePoint,const Point &p0,const Point &p1,const Point &p2,const Point &p, VertCoord vec[2][3],VertCoord *tempVertCoords) {
  // lets compute the cross product of two vectors of the triangle in order to get the plane equation A(x-x0)+ B(y-y0)+C(z-z0)= 0
  //more info: http://tutorial.math.lamar.edu/Classes/CalcIII/EqnsOfPlanes.aspx
  //VertCoord vec[2][3];
  vec[0][0] = p1[0]; vec[0][0] -= p0[0]; //vec[0][0] = p1[0]-p0[0];
  vec[0][1] = p1[1]; vec[0][1] -= p0[1]; //vec[0][1] = p1[1]-p0[1];
  vec[0][2] = p1[2]; vec[0][2] -= p0[2];  //vec[0][2] = p1[2]-p0[2];

  vec[1][0] = p2[0]; vec[1][0] -= p0[0]; //vec[1][0] = p2[0]-p0[0];
  vec[1][1] = p2[1]; vec[1][1] -= p0[1]; //vec[1][1] = p2[1]-p0[1];
  vec[1][2] = p2[2]; vec[1][2] -= p0[2]; //vec[1][2] = p2[2]-p0[2];

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
  heightAbovePoint = tempVertCoords[3];

  tempVertCoords[3] = p0[1];
  tempVertCoords[3] -= p[1];
  tempVertCoords[3] *= tempVertCoords[1];
  heightAbovePoint += tempVertCoords[3];

  tempVertCoords[3] = p0[2];
  tempVertCoords[3] *= tempVertCoords[2];
  heightAbovePoint += tempVertCoords[3];

  heightAbovePoint /= tempVertCoords[2];
}

void PinMesh::computeBarycentricCoordinates(const Point &p0,const Point &p1,const Point &p2, VertCoord &A, VertCoord &B, VertCoord &C) {
  VertCoord vec[2][3];

  vec[0][0] = p1[0]-p0[0];
  vec[0][1] = p1[1]-p0[1];
  vec[0][2] = p1[2]-p0[2];

  vec[1][0] = p2[0]-p0[0];
  vec[1][1] = p2[1]-p0[1];
  vec[1][2] = p2[2]-p0[2];

  //lets compute the values A,B,C basing on the two 3D vectors vec[0] and vec[1]
  //  | i             j           k      i             j       
  //  | vec[0][0]  vec[0][1]  vec[0][2]  vec[0][0]  vec[0][1]  
  //  | vec[1][0]  vec[1][1]  vec[1][2]  vec[1][0]  vec[1][1] 

  A = vec[0][1]*vec[1][2] - vec[0][2]*vec[1][1];
  //VertCoord B = vec[0][0]*vec[1][2] - vec[0][2]*vec[1][0];
  B = vec[0][2]*vec[1][0] - vec[0][0]*vec[1][2];
  C = vec[0][0]*vec[1][1] - vec[0][1]*vec[1][0];
  assert(C!=0);
}



const Triangle * PinMesh::getBestTrianglePointInObjectSoS(int meshId,const Triangle *newTriangle,const Triangle *bestTriangleSoFar,const Point &p) {
  //todo : consider meshId...
  VertCoord A,B,C;
  computeBarycentricCoordinates(vertices[meshId][(*bestTriangleSoFar)[0]],vertices[meshId][(*bestTriangleSoFar)[1]],vertices[meshId][(*bestTriangleSoFar)[2]], A,B,C);
  VertCoord Anew,Bnew,Cnew;
  computeBarycentricCoordinates(vertices[meshId][(*newTriangle)[0]],vertices[meshId][(*newTriangle)[1]],vertices[meshId][(*newTriangle)[2]], Anew,Bnew,Cnew);

  //since the two triangles have the same height --> they should share the query point p (in a common vertex/edge)
  //we want to check if..
  // (Anew*epsilon+Bnew*epsilon^2)/Cnew > (A*epsilon + B*epsilon^2)/C
  //They should never be equal (otherwise the triangles would be coplanar)

  if ( (Anew/Cnew) > (A/C)) return newTriangle;
  else if ( ((Anew/Cnew) == (A/C)) && (Bnew/Cnew > B/C)) return newTriangle;
  return bestTriangleSoFar;
}

bool PinMesh::triangleAbovePointSoS(int meshId,const Triangle &triangle,const Point &p) {
  VertCoord A,B,C;
  computeBarycentricCoordinates(vertices[meshId][triangle[0]],vertices[meshId][triangle[1]],vertices[meshId][triangle[2]], A,B,C);

  // now we have an equation A(x-triangle[0][0]) +  B(y-triangle[0][1]) + C(z-triangle[0][2]) = 0
  // A(x-triangle[0][0]) +  B(y-triangle[0][1]) + C*z - C*triangle[0][2] = 0
  // C*z = C*triangle[0][2] - A(x-triangle[0][0]) - B(y-triangle[0][1])
  // C*z = A(triangle[0][0]-x) +  B(triangle[0][1]-y) + C*triangle[0][2] 
  // z = (A(triangle[0][0]-x) +  B(triangle[0][1]-y) + C*triangle[0][2]) / C

  //invert answers if C<0... so replace true with (C>=0) and false with (C<0)
  //cerr << sgn(A) << " " << sgn(B) <<  " " << sgn(C) << endl;
  //cerr << A << " " << B << " " << C << endl;

  if (sgn(A)<0) return sgn(C)>=0; // true
  else if (sgn(A)==0) {
    if (sgn(B)<0) return sgn(C)>=0; // true
    else if (sgn(B)==0) return false;
    else return sgn(C)<0; // false 
  } else return sgn(C)<0; //false

}


*/

//lets suppose the uniform grid works with only one level...
//given a point, the coordinates of the grid cell where p is and a 1-level uniform grid, computes the object where p is
//if p is outside the objects, returns OUTSIDE_OBJECT
//meshId (0 or 1) represents the map that we want to verify 

//TODO:
//special cases.
//example of special case: vertical triangle, point below vertex of triangle, point below edge of triangle...
//temp_big_ints should have size at least 2
//tempVertCoords should have size at least 8
ObjectId PinMesh::computeObjectWherePointIsTwoLevel(const InputVertex &p,int globalGridCoordX,int globalGridCoordY, int globalGridCoordZ, int meshId, TempVarsComputeObjectWherePointIs &tempVars, bool &foundUsingGrid) {
  const int gridSizeLevel1 = uniformGrid->gridSizeLevel1;
  const int gridSizeLevel2 = uniformGrid->gridSizeLevel2;

  foundUsingGrid = false;

  int xGridFirstLevel = globalGridCoordX/gridSizeLevel2;
  int yGridFirstLevel = globalGridCoordY/gridSizeLevel2;
  int zGridFirstLevel = globalGridCoordZ/gridSizeLevel2;

  int xCellCoordSubGrid = globalGridCoordX%gridSizeLevel2;
  int yCellCoordSubGrid = globalGridCoordY%gridSizeLevel2;
  int zCellCoordSubGrid = globalGridCoordZ%gridSizeLevel2;

  if (xGridFirstLevel <0 || yGridFirstLevel<0 || zGridFirstLevel <0 || xGridFirstLevel>= gridSizeLevel1 || yGridFirstLevel>=gridSizeLevel1 || zGridFirstLevel>=gridSizeLevel1)
    return OUTSIDE_OBJECT; // if we are outside the bounding box --> we are outside any object.
  //cerr << "Computing.." << endl;
  int highestCellZToProcess = gridSizeLevel1;
  bool foundATriangleAboveP = false; //have we already found a triangle above p?

  

  const InputTriangle *bestTriangle = NULL; //the lowest triangle above p
  MeshIntersectionGeometry::HeightPointInTriangleProjection &heightAbovePoint = tempVars.heightAbovePoint;
  MeshIntersectionGeometry::HeightPointInTriangleProjection &heightAbovePointBestTriangle = tempVars.heightAbovePointBestTriangle;
  
  bool firstIteration = true;

  

  //const vector<vector<vector<Nested3DGridCell> > >  &gridCells = uniformGrid->grid.gridCells;
  const Nested3DGrid &gridFirstLevel = uniformGrid->grid;

  //lets process all cells above (and including...) the cell containing p in order to try to find the lowest triangle above p
  for(int cz = zGridFirstLevel;cz<gridSizeLevel1;cz++) {
    //To determine if a point p is above/below a triangle, we just need to project the triangle and the point in a horizontal
    //plane and check if the point is inside the projected triangle 

    //After checking if the triangle t is above or below the point p, we need to project p in t and get the z coordinate of this projection
    //we want to find the triangle with lowest z coordinate that is greater than p's z coordinate 
    //const Nested3DGridCell &gridCell = gridCells[xGridFirstLevel][yGridFirstLevel][cz];

    //const vector<Triangle *> &trianglesInCell = gridCell.triangles[meshId];
    InputTriangle * *trianglesInCell = gridFirstLevel.getPointerStartListTriangles(meshId,gridSizeLevel1,xGridFirstLevel,yGridFirstLevel,cz);
    //int numTrianglesInCell = trianglesInCell.size();
    int numTrianglesInCell = gridFirstLevel.numTrianglesInGridCell(meshId,gridSizeLevel1,xGridFirstLevel,yGridFirstLevel,cz);



    /*if (!foundATriangleAboveP && gridCellsLabels.labels[xGridFirstLevel][yGridFirstLevel][cz] >= 0) {
      foundUsingGrid = true;
      return gridCellsLabels.labels[xGridFirstLevel][yGridFirstLevel][cz];
    }*/

    if (!gridFirstLevel.hasSecondLevel(xGridFirstLevel,yGridFirstLevel,cz)) { 
      for (int iTriangle = 0;iTriangle < numTrianglesInCell;iTriangle++) {
        //cerr << "Testing triangle: " << iTriangle << " of " << numTrianglesInCell << endl;
        const InputTriangle &triangle = *(trianglesInCell[iTriangle]);

        if(&triangle==bestTriangle) continue; //the same triangle may be in different grid cells...
        assert(bestTriangle==NULL || (&triangle)->compare(*bestTriangle)!=0);

        //assert(pointInTriangleProj(vertices[meshId][triangle[0]],vertices[meshId][triangle[1]],vertices[meshId][triangle[2]],p, tempVertCoords + 3) == pointInTriangleProj(vertices[meshId][triangle[0]],vertices[meshId][triangle[1]],vertices[meshId][triangle[2]],p));
        //if (!pointInTriangleProj(vertices[meshId][triangle[0]],vertices[meshId][triangle[1]],vertices[meshId][triangle[2]],p, tempVertCoords + 3) ) continue; //the triangle is not above/below p
        if(!geometry->isVertexInTriangleProjection(triangle, p, tempVars.tempVarsIsVertexTriangleProjectionZ0)) 
            continue;

        //getHeigthAbovePoint(heightAbovePoint,vertices[meshId][triangle[0]], vertices[meshId][triangle[1]],vertices[meshId][triangle[2]],p , tempVertCoordMatrix, tempVertCoords + 3); //get the height of the triangle above p (considering that the point is projected vertically in the triangle)
        //for performance, we will first compute the height ignoring SoS... If the height is not equal --> SoS will not change the result!
        geometry->computeHeightAbovePointNoSoS(heightAbovePoint,triangle,p, tempVars.tempVarComputeHeightAbovePointNoSoS);

        //1 if heightAbovePoint > p (considering the z-coordinate)
        //0 if equal
        int heightComparison = geometry->compareHeightWithPointHeightNoSoS(p,heightAbovePoint);

        //p.print(); cerr << endl;
        //cerr << "Height comparison... " << endl;

        //if the height of the point "above" (actually, above or below) is greater than the height of p --> the point is above p
        //if the height comparison is 0, we need SoS to make sure the triangle is above the point (after SoS)
        if(heightComparison==1 || (heightComparison==0 && geometry->isTriangleAbovePointSoS(triangle,p,tempVars.tempVarIsTriangleAbovePointSoS)) ) {  
          int heightComparisonWithBestHeight = geometry->compareHeightWithBestHeightNoSoS(heightAbovePoint,heightAbovePointBestTriangle);   
          //1 if height is smaller than the best, 0 if equal...
          if (!foundATriangleAboveP || heightComparisonWithBestHeight==1) {
            //highestCellZToProcess  = uniformGrid->z_cell_from_coord_level1(heightAbovePoint, tempVertCoord,tempBigInts) ;
            highestCellZToProcess  = geometry->zCellLevel1FromProjectionOfPoint(heightAbovePoint, triangle, p, *uniformGrid, tempVars.tempVarZCellFromProjectionOfPoint) ;
            
            foundATriangleAboveP = true;
            bestTriangle = &triangle;
            heightAbovePointBestTriangle = heightAbovePoint;
          } else if (heightComparisonWithBestHeight == 0) { //the same triangle may be in different cells and considered twice 
            //We have one or more triangles with the same height...we need to use SoS to choose the correct one
            // we don't need to update the "highestCellZToProcess" since the height of the projection is the same as a previous one...
            // we only need to update the "best triangle"
            //cerr << "Same height" << endl;
            bestTriangle = geometry->getBestTrianglePointInObjectSoS(&triangle,bestTriangle,p, tempVars.tempVarGetBestTrianglePointInObjectSoS);
          }
        } 
      }
    } else { //we have a nested grid....
      const Nested3DGrid &childGrid = *gridFirstLevel.getChildGrid(xGridFirstLevel,yGridFirstLevel,cz);
      
      //in the first iteration we should not start in the coordinate 0 (because the vertex is inside the divided)
      int startZ = 0; 
      if (firstIteration) {        
        startZ = zCellCoordSubGrid;//childGrid.z_cell_from_coord(p[2],tempVertCoord,tempBigInts);      
      }
      
      int highestCellZToProcessNestedSubGrid = gridSizeLevel2;
      for(int zSubGrid=startZ;zSubGrid<gridSizeLevel2;zSubGrid++) { //iterate through the subgrid...
        /*const int ipol = gridCellsLabels.childGridLabels[xGridFirstLevel][yGridFirstLevel][cz]->labels[xCellCoordSubGrid][yCellCoordSubGrid][zSubGrid] ;//grid[1-imap][cx][cy].pol; 
        if (!foundATriangleAboveP && ipol >= 0) {
          foundUsingGrid = true;
          return ipol;
        }*/


        //const vector<Triangle *> &trianglesInNestedCell = childGrid.gridCells[xCellCoordSubGrid][yCellCoordSubGrid][zSubGrid].triangles[meshId];
        //int numTrianglesInNestedCell = trianglesInNestedCell.size();
        InputTriangle * *trianglesInNestedCell = childGrid.getPointerStartListTriangles(meshId,gridSizeLevel2,xCellCoordSubGrid,yCellCoordSubGrid,zSubGrid);
        int numTrianglesInNestedCell = childGrid.numTrianglesInGridCell(meshId,gridSizeLevel2,xCellCoordSubGrid,yCellCoordSubGrid,zSubGrid);


        for (int iTriangle = 0;iTriangle < numTrianglesInNestedCell;iTriangle++) {
          const InputTriangle &triangle = *(trianglesInNestedCell[iTriangle]);

          if(&triangle==bestTriangle) continue; //the same triangle may be in different grid cells...

         //cerr << "Testing: " << vertices[meshId][triangle[0]][0] << "," << vertices[meshId][triangle[0]][1] << "," << vertices[meshId][triangle[0]][2] << endl;
         // cerr << vertices[meshId][triangle[1]][0] << "," << vertices[meshId][triangle[1]][1] << "," << vertices[meshId][triangle[1]][2] << endl;
         // cerr << vertices[meshId][triangle[2]][0] << "," << vertices[meshId][triangle[2]][1] << "," << vertices[meshId][triangle[2]][2] << endl << endl;

          //if (!pointInTriangleProj(vertices[meshId][triangle[0]],vertices[meshId][triangle[1]],vertices[meshId][triangle[2]],p, tempVertCoords + 3) ) continue; //the triangle is not above/below p
          if(!geometry->isVertexInTriangleProjection(triangle, p, tempVars.tempVarsIsVertexTriangleProjectionZ0)) 
            continue;

          //getHeigthAbovePoint(heightAbovePoint,vertices[meshId][triangle[0]], vertices[meshId][triangle[1]],vertices[meshId][triangle[2]],p , tempVertCoordMatrix, tempVertCoords + 3); //get the height of the triangle above p (considering that the point is projected vertically in the triangle)
          //cerr << "Height: " << heightAbovePoint << endl;

          geometry->computeHeightAbovePointNoSoS(heightAbovePoint,triangle,p, tempVars.tempVarComputeHeightAbovePointNoSoS);

          //1 if heightAbovePoint > p (considering the z-coordinate)
          //0 if equal
          int heightComparison = geometry->compareHeightWithPointHeightNoSoS(p,heightAbovePoint);

          //cerr << endl;
          //triangle.print(); cerr << endl;
          //cerr << "Height comparison... " << heightComparison << endl;
          //cerr << "Is above point?  " << geometry->isTriangleAbovePointSoS(triangle,p,tempVars.tempVarIsTriangleAbovePointSoS) << endl;
         // cerr << "Positive or negative Z? " <<geometry->isTriangleNormalPointingPositiveZ(triangle, tempVars.tempVarIsTriangleNormalPointingPositiveZ) << endl;

          //if the height of the point "above" (actually, above or below) is greater than the height of p --> the point is above p
          if(heightComparison == 1 || (heightComparison == 0 && geometry->isTriangleAbovePointSoS(triangle,p,tempVars.tempVarIsTriangleAbovePointSoS))) {        
            //highestCellZToProcess =  min(highestCellZToProcess,uniformGrid->z_cell_from_coord(heightAbovePoint, tempVertCoord,tempBigInts) );
            int heightComparisonWithBestHeight = geometry->compareHeightWithBestHeightNoSoS(heightAbovePoint,heightAbovePointBestTriangle);   
            //1 if height is smaller than the best, 0 if equal...
            if (!foundATriangleAboveP || heightComparisonWithBestHeight==1) {
              //highestCellZToProcess  = uniformGrid->z_cell_from_coord_level1(heightAbovePoint, tempVertCoord,tempBigInts) ;
              int czGlobal = geometry->zCellGlobalFromProjectionOfPoint(heightAbovePoint, triangle, p, *uniformGrid,tempVars.tempVarZCellGlobalFromProjectionOfPoint) ;
              
              int cz2 = czGlobal/gridSizeLevel2;//uniformGrid->z_cell_from_coord_level1(heightAbovePoint, tempVertCoord,tempBigInts);

              highestCellZToProcess = cz2; //we dont need to process cells above cz2..
              //this triangle is the lowest one seen so far...
              //we do not need to process any grid cell above this point...

              foundATriangleAboveP = true;
              bestTriangle = &triangle;
              heightAbovePointBestTriangle = heightAbovePoint;

              if (cz2 <= cz) { //the lowest point so far above p is inside the nested grid... lets update the highest second level grid cell to process...
                highestCellZToProcessNestedSubGrid =  czGlobal%gridSizeLevel2;//childGrid.z_cell_from_coord(heightAbovePoint, tempVertCoord,tempBigInts);
              }
            } else if (heightComparisonWithBestHeight == 0) { //the same triangle may be in different cells and considered twice 
              //We have one or more triangles with the same height...we need to use SoS to choose the correct one
              // we don't need to update the "highestCellZToProcess" since the height of the projection is the same as a previous one...
              // we only need to update the "best triangle"
              //cerr << "Same height" << endl;
              bestTriangle = geometry->getBestTrianglePointInObjectSoS(&triangle,bestTriangle,p, tempVars.tempVarGetBestTrianglePointInObjectSoS);
            }         
          } 

        }
        if(foundATriangleAboveP && cz >= highestCellZToProcess  && zSubGrid >= highestCellZToProcessNestedSubGrid) {
           break;
        } 
      }
    }

    firstIteration = false;
    if (foundATriangleAboveP && cz >= highestCellZToProcess) break;
  }

  //cerr << "Best height: " << heightAbovePointBestTriangle << endl;
  /*cerr << "Found a triangle? " << foundATriangleAboveP << endl;
  if(foundATriangleAboveP) {
    bestTriangle->print(); cerr << endl;
    cerr << "Positive z? : " << geometry->isTriangleNormalPointingPositiveZ(*bestTriangle, tempVars.tempVarIsTriangleNormalPointingPositiveZ) << endl;
    cerr << "Above, below: " << bestTriangle->above << " , "<< bestTriangle->below << endl;
  }*/
  if (foundATriangleAboveP) {   
    if(geometry->isTriangleNormalPointingPositiveZ(*bestTriangle, tempVars.tempVarIsTriangleNormalPointingPositiveZ))  
      return bestTriangle->below;
    else 
      return bestTriangle->above;
  }
  else return OUTSIDE_OBJECT; //there is no triangle above p ... p must not be inside any object from the other map
}











//Locate a set of vertices in the objects in map 0...
//if meshIdToLocate=0, this means that vertices will be located in mesh 0
//#define PINMESH_VERBOSE
//the vectors gridXGlobalVertexToLocate (idem Y,Z) indicate in what grid cell (global coord) each vertex to locate is
void PinMesh::locateVerticesInObject(const vector<InputVertex> &verticesToLocate,  vector<ObjectId> &verticesIds,int meshIdToLocate) {
  timespec t0,t1,t01;

  int sum = 0; //we use this because the compiler may remove the call if we do not use the result...

  verticesIds.resize(verticesToLocate.size());

  int numVerticesToLocate = verticesToLocate.size();
  int gridSize = uniformGrid->gridSizeLevel1;
  int nestedGridSize = uniformGrid->gridSizeLevel2;
  
  #ifdef PINMESH_VERBOSE
  cerr << "Computing CCs of empty boxes..." << endl;
  clock_gettime(CLOCK_REALTIME, &t0);
  #endif

  //we will compute where the empty cells are (to accelerate the point in vol computations!)


/*  big_int tempVarsInt[2];
  VertCoord tempVertCoordMatrix[2][3];
  VertCoord tempVertCoords[8];
  VertCoord tempVar;

  Point box0 = uniformGrid->box[0];
  VertCoord cellWidth = uniformGrid->cellWidthLevel1 ;
  VertCoord nestedGridCellWidth  = uniformGrid->cellWidthLevel2 ;
  VertCoord halfCellWidth = cellWidth/2;
  VertCoord halfCellWidthNestedGrid = nestedGridCellWidth/2;*/
  

  //const vector<vector<vector<Nested3DGridCell> > >  &gridCells = uniformGrid->grid.gridCells;
  /*GridCellsLabels cellsLabels(gridSize); 

  const Nested3DGrid &gridFirstLevel = (uniformGrid->grid);
  for(int gx=0;gx<gridSize;gx++) 
    for(int gy=0;gy<gridSize;gy++)
      for(int gz=0;gz<gridSize;gz++) {
        int numTrianglesInCell = gridFirstLevel.numTrianglesInGridCell(meshIdToLocate,gridSize,gx,gy,gz);

        if(gridFirstLevel.hasSecondLevel(gx,gy,gz)) { //we have a nested grid...
            const Nested3DGrid &nestedGrid = *gridFirstLevel.childGrids[gx][gy][gz];//*gridCells[gx][gy][gz].childGrid;

            cellsLabels.childGridLabels[gx][gy][gz] = new GridCellsLabels(nestedGridSize);

            vector<vector<vector<int > > > &labels = cellsLabels.childGridLabels[gx][gy][gz]->labels;

            for(int nx=0;nx<nestedGridSize;nx++)  //lets initialize the nested grid labels...
              for(int ny=0;ny<nestedGridSize;ny++)
                for(int nz=0;nz<nestedGridSize;nz++)
                  if (nestedGrid.numTrianglesInGridCell(meshIdToLocate, nestedGridSize,nx,ny,nz) !=0)//if (nestedGrid.gridCells[nx][ny][nz].triangles[0].size()!=0)
                    labels[nx][ny][nz] = LABEL_NON_EMPTY_CELL;

        } else { //this grid cell does not have a nested grid... and there is at least one triangle in this grid cell...
          if(numTrianglesInCell !=0) //if (trianglesInCell.size()!=0)
            cellsLabels.labels[gx][gy][gz] = LABEL_NON_EMPTY_CELL;
        }        
      }
 // cerr << "End of filling labels.." << endl;
  clock_gettime(CLOCK_REALTIME, &t01);  */

  int numSearchsPerformedToFillGrid = 0;
  


 // #define DONT_USE_GRID_ACCELERATION
  /*
  #ifndef DONT_USE_GRID_ACCELERATION  


  assert( (gridSize%4) == 0);
  vector< array<array<int,3>,2> > startEndGridScan;
  array<int,3> start,end;

  int inc = gridSize/4;



  for(int gx=0;gx<gridSize;gx+=inc) 
      for(int gy=0;gy<gridSize;gy+=inc)
        for(int gz=0;gz<gridSize;gz+=inc) {
          start[0]=gx;
          start[1] = gy;
          start[2] = gz;

          end[0] = gx+inc;
          end[1] = gy+inc;
          end[2] = gz+inc;
          startEndGridScan.push_back({start,end});
        }

  int szBlocksProcess = startEndGridScan.size();


  

  #pragma omp parallel for schedule(dynamic,1)
  for(int i=0;i<szBlocksProcess;i++) {    
    array<int,3> start,end;
    start = startEndGridScan[i][0];
    end = startEndGridScan[i][1];

    big_int tempVarsInt[2];
    VertCoord tempVertCoordMatrix[2][3];
    VertCoord tempVertCoords[8];
    VertCoord tempVar;

    for(int gx=start[0];gx<end[0];gx++) 
      for(int gy=start[1];gy<end[1];gy++)
        for(int gz=start[2];gz<end[2];gz++) {
          if(!gridFirstLevel.hasSecondLevel(gx,gy,gz)) { //if(!gridCells[gx][gy][gz].hasChild()) {       
            
            if(cellsLabels.labels[gx][gy][gz]==DONT_KNOW_ID) {
              //cerr << "Computing id..." << endl;
              Point p;
              p[0] = gx*cellWidth + box0[0] + halfCellWidth;
              p[1] = gy*cellWidth + box0[1] + halfCellWidth;
              p[2] = gz*cellWidth + box0[2] + halfCellWidth;
              bool foundUsingGrid;
              //cerr << "computing id level 1......" << endl;
              ObjectId id = computeObjectWherePointIsTwoLevel(p,gx*nestedGridSize+nestedGridSize/2,gy*nestedGridSize+nestedGridSize/2,gz*nestedGridSize+nestedGridSize/2 ,meshIdToLocate, uniformGrid, tempVertCoords, tempVertCoordMatrix, tempVarsInt,cellsLabels,foundUsingGrid);
              //cerr << "setiing to " << id <<  endl;
              setObjectInWhereEmptyCellIs2(gx,gy,gz,-1,-1,-1,gridSize,nestedGridSize,cellsLabels,id,start[0],start[1],start[2],end[0],end[1],end[2]);
              //cellsLabels.labels[gx][gy][gz] = id;
             // cerr << "set level 1" << endl;

              numSearchsPerformedToFillGrid++;
              //cerr << "end" << endl;
            }   

          } else { //there is a nested uniform grid...
            //cerr << "The grid is nested... " << endl;
            //const Nested3DGrid &nestedGrid = *gridCells[gx][gy][gz].childGrid;
            vector<vector<vector<int > > > &labels = cellsLabels.childGridLabels[gx][gy][gz]->labels;

            //cerr << "Computing... " << endl;
            for(int nx=0;nx<nestedGridSize;nx++)  //lets initialize the nested grid labels...
              for(int ny=0;ny<nestedGridSize;ny++)
                for(int nz=0;nz<nestedGridSize;nz++)
                  if (labels[nx][ny][nz]==DONT_KNOW_ID) {
                    //cerr << "Aqui " << endl;
                    Point p;
                    p[0] = (gx*cellWidth + box0[0]) + nestedGridCellWidth*nx + halfCellWidthNestedGrid;
                    p[1] = (gy*cellWidth + box0[1]) + nestedGridCellWidth*ny + halfCellWidthNestedGrid;
                    p[2] = (gz*cellWidth + box0[2]) + nestedGridCellWidth*nz + halfCellWidthNestedGrid;
                    bool foundUsingGrid;
                   // cerr << "computing id level 2..." ;
                    ObjectId id = computeObjectWherePointIsTwoLevel(p,gx*nestedGridSize + nx,gy*nestedGridSize + ny,gz*nestedGridSize +nz,meshIdToLocate, uniformGrid, tempVertCoords, tempVertCoordMatrix, tempVarsInt,cellsLabels,foundUsingGrid);
                   // cerr << "setiing 2 to " << id <<  endl;
                    setObjectInWhereEmptyCellIs2(gx,gy,gz,nx,ny,nz,gridSize,nestedGridSize,cellsLabels,id,start[0],start[1],start[2],end[0],end[1],end[2]);
                   // cerr << "set level 2" << endl;
                  //  labels[nx][ny][nz] = id;

                    numSearchsPerformedToFillGrid++;
                  }
            //cerr << "End grid nested..." << endl;
          }                             
        } 
  }




  #endif 
  */

//  cerr << "Locating two level......" << endl;
  
  #ifdef PINMESH_VERBOSE
  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to fill grid and compute CCs: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;  
  cerr << "Time only to compute CCs: " << convertTimeMsecs(diff(t01,t1))/1000 << endl;
  #endif

  #ifdef PINMESH_VERBOSE   
  //computing some statistics...
  /*int numEmptyGridCells=0,numNonEmptyGridCells= 0;
   for(int gx=0;gx<gridSize;gx++) 
      for(int gy=0;gy<gridSize;gy++)
        for(int gz=0;gz<gridSize;gz++) {
          if(!gridFirstLevel.hasSecondLevel(gx,gy,gz)) { //if(!gridCells[gx][gy][gz].hasChild()) {       
            
            if(cellsLabels.labels[gx][gy][gz]==-3) {
              numNonEmptyGridCells++;
            }   else {
              numEmptyGridCells++;
            }

          } else { 
            vector<vector<vector<int > > > &labels = cellsLabels.childGridLabels[gx][gy][gz]->labels;

            //cerr << "Computing... " << endl;
            for(int nx=0;nx<nestedGridSize;nx++)  //lets initialize the nested grid labels...
              for(int ny=0;ny<nestedGridSize;ny++)
                for(int nz=0;nz<nestedGridSize;nz++)
                  if (labels[nx][ny][nz]==-3) {
                    numNonEmptyGridCells++;
                  }   else {
                    numEmptyGridCells++;
                  }
          }                             
        } 
   
  cerr << "Empty grid cells: " << numEmptyGridCells << endl;
  cerr << "Non empty grid cells: " << numNonEmptyGridCells << endl;
  cerr << "Percent empty cells: " << (100.0*numEmptyGridCells)/(numEmptyGridCells+numNonEmptyGridCells) << endl; */
  #endif


  #ifdef PINMESH_VERBOSE
  cerr << "Locating " << verticesToLocate.size() << " vertices..." << endl;
  clock_gettime(CLOCK_REALTIME, &t0);
  #endif
  

  int numVerticesFoundUsingGrid = 0;

  #pragma omp parallel
  {
    TempVarsComputeObjectWherePointIs tempVarsLocate;
    int numVerticesFoundUsingGridLocal = 0;
    bool foundUsingGrid;
    
    #pragma omp for schedule(dynamic,1000)
    for(int i=0;i<numVerticesToLocate;i++) {
      const InputVertex &p = verticesToLocate[i];        
    
      int gx = uniformGrid->get_global_x_coord_mesh_vertex(p.getMeshId(),p.getId());//gridXGlobalVertexToLocate[i];// gxVector[i];
      int gy = uniformGrid->get_global_y_coord_mesh_vertex(p.getMeshId(),p.getId());//gridYGlobalVertexToLocate[i];//gyVector[i];
      int gz = uniformGrid->get_global_z_coord_mesh_vertex(p.getMeshId(),p.getId());//gridZGlobalVertexToLocate[i];//gzVector[i];


      ObjectId id2 = computeObjectWherePointIsTwoLevel(p,gx,gy,gz,meshIdToLocate, tempVarsLocate ,foundUsingGrid);
      numVerticesFoundUsingGridLocal += foundUsingGrid;

      /*cerr << "Vertex from mesh: " << 1-meshIdToLocate << endl;
      cerr << "Point from mesh: " << p.getMeshId() << endl;
      cerr << "Locating... "; geometry->printPointForDebugging(p); cerr << endl;
      cerr << "Location: " << id2 << endl;
      cerr << "Ouside object is: "<< OUTSIDE_OBJECT << endl << endl;;*/

      verticesIds[i] = id2;
    }
    #pragma omp critical
    numVerticesFoundUsingGrid += numVerticesFoundUsingGridLocal;
  }
  

  #ifdef PINMESH_VERBOSE
  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to locate vertices: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  cerr << "Num vertices found using grid: " << numVerticesFoundUsingGrid << endl;
  cerr << "Percentage vertices found using grid: " << (100.0*numVerticesFoundUsingGrid)/numVerticesToLocate << "\n\n";
  #endif

  
}


