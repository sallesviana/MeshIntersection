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

#include <iostream>
#include <vector>
#include <array>
#include <cstdlib>
#include <time.h>
#include "rationals.h"




using namespace std;




// Utils for measuring time...

double convertTimeMsecs(const timespec &t){
  return (t.tv_sec*1000 + t.tv_nsec/1000000.0);
}

//source: http://www.guyrutenberg.com/2007/09/22/profiling-code-using-clock_gettime/
timespec diff(timespec start, timespec end)
{
        timespec temp;
        if ((end.tv_nsec-start.tv_nsec)<0) {
          temp.tv_sec = end.tv_sec-start.tv_sec-1;
          temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
        } else {
          temp.tv_sec = end.tv_sec-start.tv_sec;
          temp.tv_nsec = end.tv_nsec-start.tv_nsec;
        }
        return temp;
}

//===============================================================
// Constants...
//===============================================================



//===============================================================
#include "3d_objects.cpp"

#include "nested3DGrid.cpp"

#include "tritri_isectline.c"


//Each vector represents the vertices of a layer
vector<Point> vertices[2];

//Each vector represents a set of objects in the same layer
//The objects are represented by a set of triangles (defining their boundaries)
vector<Triangle> triangles[2]; //

Point bBox[2] ;//bounding box of each layer (each boundingbox has two vertices)
Point boundingBoxLayers[2]; //bounding box considering both layers togetter

Nested3DGridWrapper uniformGrid;



//=======================================================================================================================

vector< pair< array<VertCoord,3>,array<VertCoord,3> > > edges;

/*
VertCoord tempRationals[100];

bool intersect(const Triangle &a, const Triangle &b) {
  //return tri_tri_intersect(  vertices[0][a.p[0]].data(),vertices[0][a.p[1]].data(),vertices[0][a.p[2]].data()     ,     vertices[1][b.p[0]].data(),vertices[1][b.p[1]].data(),vertices[1][b.p[2]].data());
  int coplanar;
  VertCoord p0InterLine[3],p1InterLine[3];
  int ans = tri_tri_intersect_with_isectline(vertices[0][a.p[0]].data(),vertices[0][a.p[1]].data(),vertices[0][a.p[2]].data()     ,     
                                              vertices[1][b.p[0]].data(),vertices[1][b.p[1]].data(),vertices[1][b.p[2]].data(),
                                              &coplanar, p0InterLine ,p1InterLine,tempRationals);

  if (ans && !coplanar) {
    edges.push_back( pair<array<VertCoord,3>,array<VertCoord,3>>( {p0InterLine[0],p0InterLine[1],p0InterLine[2]} , {p1InterLine[0],p1InterLine[1],p1InterLine[2]}   ) );
    //cerr << p0InterLine[0] << "  " << p0InterLine[1] << " " << p0InterLine[2] << endl;
    //cerr << p1InterLine[0] << "  " << p1InterLine[1] << " " << p1InterLine[2] << endl;
  }

  return ans;
}
*/



//=======================================================================================================================








//Returns true iff the projection of the point p into a horizontal plane is inside of the projection of triangle
// (p0,p1,p2) into the horizontal plane...

//returns true if the point p is vertically directly above/below the point p...
//TODO special cases:
//- degenerate triangle?
//- point above/below vertex/edge?
//- vertical triangle?
bool pointInTriangleProj(const Point &p0, const Point &p1, const Point &p2, const Point &p)  {
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
bool pointInTriangleProjSoS(const Point &p0, const Point &p1, const Point &p2, const Point &p)  {	
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
//TODO special cases:
//- degenerate triangle?
//- point above/below vertex/edge?
//- vertical triangle?
/*bool pointInTriangleProj(const Triangle &triangle,const Point &p) { 
  VertCoord denominator = ((triangle[1][1] - triangle[2][1])*(triangle[0][0] - triangle[2][0]) + (triangle[2][0] - triangle[1][0])*(triangle[0][1] - triangle[2][1]));
  VertCoord a = ((triangle[1][1] - triangle[2][1])*(p[0] - triangle[2][0]) + (triangle[2][0] - triangle[1][0])*(p[1] - triangle[2][1])) / denominator;
  if ( a<0 || a >1) return false;
  
  VertCoord b = ((triangle[2][1] - triangle[0][1])*(p[0] - triangle[2][0]) + (triangle[0][0] - triangle[2][0])*(p[1] - triangle[2][1])) / denominator;

  if ( b<0 || b>1) return false;

  VertCoord c = 1 - a - b;

  return 0 <= c && c <= 1; 


/*  point_coord denominator = ((p1.second - p2.second)*(p0.first - p2.first) + (p2.first - p1.first)*(p0.second - p2.second));
  point_coord a = ((p1.second - p2.second)*(p.first - p2.first) + (p2.first - p1.first)*(p.second - p2.second)) / denominator;
  if ( a<0 || a >1) return false;
  
  point_coord b = ((p2.second - p0.second)*(p.first - p2.first) + (p0.first - p2.first)*(p.second - p2.second)) / denominator;

  if ( b<0 || b>1) return false;

  point_coord c = 1 - a - b;

  return 0 <= c && c <= 1; 
}*/

//p0,p1,p2 are the vertices of the triangle...
void getHeigthAbovePoint(VertCoord &heightAbovePoint,const Point &p0,const Point &p1,const Point &p2,const Point &p) {
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
bool pointInTriangleProj(const Point &p0, const Point &p1, const Point &p2, const Point &p, VertCoord *tempVertCoords)  {
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
/*bool pointInTriangleProj(const Triangle &triangle,const Point &p) { 
  VertCoord denominator = ((triangle[1][1] - triangle[2][1])*(triangle[0][0] - triangle[2][0]) + (triangle[2][0] - triangle[1][0])*(triangle[0][1] - triangle[2][1]));
  VertCoord a = ((triangle[1][1] - triangle[2][1])*(p[0] - triangle[2][0]) + (triangle[2][0] - triangle[1][0])*(p[1] - triangle[2][1])) / denominator;
  if ( a<0 || a >1) return false;
  
  VertCoord b = ((triangle[2][1] - triangle[0][1])*(p[0] - triangle[2][0]) + (triangle[0][0] - triangle[2][0])*(p[1] - triangle[2][1])) / denominator;

  if ( b<0 || b>1) return false;

  VertCoord c = 1 - a - b;

  return 0 <= c && c <= 1; 


/*  point_coord denominator = ((p1.second - p2.second)*(p0.first - p2.first) + (p2.first - p1.first)*(p0.second - p2.second));
  point_coord a = ((p1.second - p2.second)*(p.first - p2.first) + (p2.first - p1.first)*(p.second - p2.second)) / denominator;
  if ( a<0 || a >1) return false;
  
  point_coord b = ((p2.second - p0.second)*(p.first - p2.first) + (p0.first - p2.first)*(p.second - p2.second)) / denominator;

  if ( b<0 || b>1) return false;

  point_coord c = 1 - a - b;

  return 0 <= c && c <= 1; 
}*/

//p0,p1,p2 are the vertices of the triangle...
// vec is a temporary matrix of coordinates
// tempVertCoords is a temporary array with size at least 4
void getHeigthAbovePoint(VertCoord &heightAbovePoint,const Point &p0,const Point &p1,const Point &p2,const Point &p, VertCoord vec[2][3],VertCoord *tempVertCoords) {
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

void computeBarycentricCoordinates(const Point &p0,const Point &p1,const Point &p2, VertCoord &A, VertCoord &B, VertCoord &C) {
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



//check these computations...
const Triangle * getBestTrianglePointInObjectSoS(int imapObject,const Triangle *newTriangle,const Triangle *bestTriangleSoFar,const Point &p) {
	//todo : consider imapObject...
	VertCoord A,B,C;
	computeBarycentricCoordinates(vertices[imapObject][(*bestTriangleSoFar)[0]],vertices[imapObject][(*bestTriangleSoFar)[1]],vertices[imapObject][(*bestTriangleSoFar)[2]], A,B,C);
	VertCoord Anew,Bnew,Cnew;
	computeBarycentricCoordinates(vertices[imapObject][(*newTriangle)[0]],vertices[imapObject][(*newTriangle)[1]],vertices[imapObject][(*newTriangle)[2]], Anew,Bnew,Cnew);

	//since the two triangles have the same height --> they should share the query point p (in a common vertex/edge)
	//we want to check if..
	// (Anew*epsilon+Bnew*epsilon^2)/Cnew > (A*epsilon + B*epsilon^2)/C
	//They should never be equal (otherwise the triangles would be coplanar)

	if ( (Anew/Cnew) > (A/C)) return newTriangle;
	else if ( ((Anew/Cnew) == (A/C)) && (Bnew/Cnew > B/C)) return newTriangle;
	return bestTriangleSoFar;
}

bool triangleAbovePointSoS(int imapObject,const Triangle &triangle,const Point &p) {
	VertCoord A,B,C;
	computeBarycentricCoordinates(vertices[imapObject][triangle[0]],vertices[imapObject][triangle[1]],vertices[imapObject][triangle[2]], A,B,C);

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

//lets suppose the uniform grid works with only one level...
//given a point, the coordinates of the grid cell where p is and a 1-level uniform grid, computes the object where p is
//if p is outside the objects, returns OUTSIDE_OBJECT
//imapObject (0 or 1) represents the map that we want to verify 

//TODO:
//special cases.
//example of special case: vertical triangle, point below vertex of triangle, point below edge of triangle...
//temp_big_ints should have size at least 2
//tempVertCoords should have size at least 8
ObjectId computeObjectWherePointIsTwoLevel(const Point &p,int globalGridCoordX,int globalGridCoordY, int globalGridCoordZ, int imapObject, const Nested3DGridWrapper *uniformGrid, VertCoord tempVertCoords[],VertCoord tempVertCoordMatrix[2][3], big_int tempBigInts[], GridCellsLabels &gridCellsLabels, bool &foundUsingGrid, int &numGridsCellToCheck, int &numTrianglesChecked) {
  const int gridSizeLevel1 = uniformGrid->gridSizeLevel1;
  const int gridSizeLevel2 = uniformGrid->gridSizeLevel2;
  numGridsCellToCheck = 0;
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

  

  const Triangle *bestTriangle = NULL; //the lowest triangle above p
  VertCoord &heightAbovePointBestTriangle = tempVertCoords[0];
  VertCoord &heightAbovePoint = tempVertCoords[1]; 
  VertCoord &tempVertCoord = tempVertCoords[2];
  
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

    //const vector<Triangle *> &trianglesInCell = gridCell.triangles[imapObject];
    Triangle * *trianglesInCell = gridFirstLevel.getPointerStartListTriangles(imapObject,gridSizeLevel1,xGridFirstLevel,yGridFirstLevel,cz);
    //int numTrianglesInCell = trianglesInCell.size();
    int numTrianglesInCell = gridFirstLevel.numTrianglesInGridCell(imapObject,gridSizeLevel1,xGridFirstLevel,yGridFirstLevel,cz);



    if (!foundATriangleAboveP && gridCellsLabels.labels[xGridFirstLevel][yGridFirstLevel][cz] >= 0) {
    	numGridsCellToCheck++; 
      foundUsingGrid = true;
      return gridCellsLabels.labels[xGridFirstLevel][yGridFirstLevel][cz];
    }

    if (!gridFirstLevel.hasSecondLevel(xGridFirstLevel,yGridFirstLevel,cz)) {
    	numGridsCellToCheck++; 
      for (int iTriangle = 0;iTriangle < numTrianglesInCell;iTriangle++) {
      	numTrianglesChecked++;
        //cerr << "Testing triangle: " << iTriangle << " of " << numTrianglesInCell << endl;
        const Triangle &triangle = *(trianglesInCell[iTriangle]);

        //assert(pointInTriangleProj(vertices[imapObject][triangle[0]],vertices[imapObject][triangle[1]],vertices[imapObject][triangle[2]],p, tempVertCoords + 3) == pointInTriangleProj(vertices[imapObject][triangle[0]],vertices[imapObject][triangle[1]],vertices[imapObject][triangle[2]],p));
        if (!pointInTriangleProj(vertices[imapObject][triangle[0]],vertices[imapObject][triangle[1]],vertices[imapObject][triangle[2]],p, tempVertCoords + 3) ) continue; //the triangle is not above/below p
        getHeigthAbovePoint(heightAbovePoint,vertices[imapObject][triangle[0]], vertices[imapObject][triangle[1]],vertices[imapObject][triangle[2]],p , tempVertCoordMatrix, tempVertCoords + 3); //get the height of the triangle above p (considering that the point is projected vertically in the triangle)


        //if the height of the point "above" (actually, above or below) is greater than the height of p --> the point is above p
        if(heightAbovePoint > p[2]) {        
          if (!foundATriangleAboveP || heightAbovePoint < heightAbovePointBestTriangle) {
            highestCellZToProcess  = uniformGrid->z_cell_from_coord_level1(heightAbovePoint, tempVertCoord,tempBigInts) ;
            foundATriangleAboveP = true;
            bestTriangle = &triangle;
            heightAbovePointBestTriangle = heightAbovePoint;
          } else if (heightAbovePoint == heightAbovePointBestTriangle) {
          	//We have one or more triangles with the same height...we need to use SoS to choose the correct one
          	// we don't need to update the "highestCellZToProcess" since the height of the projection is the same as a previous one...
          	// we only need to update the "best triangle"
          	//cerr << "Same height" << endl;
          	bestTriangle = getBestTrianglePointInObjectSoS(imapObject,&triangle,bestTriangle,p);
          }
        } else if (heightAbovePoint == p[2]) { //if the height is equal (the point is on the triangle) and it is the first time we see this height --> SoS
        	//cerr << "On the triangle" << endl;
        	//Decide (with SoS) is we will consider that the point is above or below the triangle...
        	if (triangleAbovePointSoS(imapObject,triangle,p)) {
        		highestCellZToProcess  = uniformGrid->z_cell_from_coord_level1(heightAbovePoint, tempVertCoord,tempBigInts) ;
            foundATriangleAboveP = true;
            bestTriangle = &triangle;
            heightAbovePointBestTriangle = heightAbovePoint;
        	}
        }
      }
    } else { //we have a nested grid....
      const Nested3DGrid &childGrid = *gridFirstLevel.getChildGrid(xGridFirstLevel,yGridFirstLevel,cz);

      /*if (xCellCoordSubGrid<0) {
        xCellCoordSubGrid = childGrid.x_cell_from_coord(p[0],tempVertCoord,tempBigInts);
        yCellCoordSubGrid = childGrid.y_cell_from_coord(p[1],tempVertCoord,tempBigInts);
      }*/

      

      //in the first iteration we should not start in the coordinate 0 (because the vertex is inside the divided)
      int startZ = 0; 
      if (firstIteration) {        
        startZ = zCellCoordSubGrid;//childGrid.z_cell_from_coord(p[2],tempVertCoord,tempBigInts);      
      }
      
      int highestCellZToProcessNestedSubGrid = gridSizeLevel2;
      for(int zSubGrid=startZ;zSubGrid<gridSizeLevel2;zSubGrid++) { //iterate through the subgrid...
      	numGridsCellToCheck++;

        const int ipol = gridCellsLabels.childGridLabels[xGridFirstLevel][yGridFirstLevel][cz]->labels[xCellCoordSubGrid][yCellCoordSubGrid][zSubGrid] ;//grid[1-imap][cx][cy].pol; 
        if (!foundATriangleAboveP && ipol >= 0) {
          foundUsingGrid = true;
          return ipol;
        }


        //const vector<Triangle *> &trianglesInNestedCell = childGrid.gridCells[xCellCoordSubGrid][yCellCoordSubGrid][zSubGrid].triangles[imapObject];
        //int numTrianglesInNestedCell = trianglesInNestedCell.size();
        Triangle * *trianglesInNestedCell = childGrid.getPointerStartListTriangles(imapObject,gridSizeLevel2,xCellCoordSubGrid,yCellCoordSubGrid,zSubGrid);
        int numTrianglesInNestedCell = childGrid.numTrianglesInGridCell(imapObject,gridSizeLevel2,xCellCoordSubGrid,yCellCoordSubGrid,zSubGrid);


        for (int iTriangle = 0;iTriangle < numTrianglesInNestedCell;iTriangle++) {
        	numTrianglesChecked++;
          const Triangle &triangle = *(trianglesInNestedCell[iTriangle]);

         //cerr << "Testing: " << vertices[imapObject][triangle[0]][0] << "," << vertices[imapObject][triangle[0]][1] << "," << vertices[imapObject][triangle[0]][2] << endl;
         // cerr << vertices[imapObject][triangle[1]][0] << "," << vertices[imapObject][triangle[1]][1] << "," << vertices[imapObject][triangle[1]][2] << endl;
         // cerr << vertices[imapObject][triangle[2]][0] << "," << vertices[imapObject][triangle[2]][1] << "," << vertices[imapObject][triangle[2]][2] << endl << endl;

          if (!pointInTriangleProj(vertices[imapObject][triangle[0]],vertices[imapObject][triangle[1]],vertices[imapObject][triangle[2]],p, tempVertCoords + 3) ) continue; //the triangle is not above/below p
          getHeigthAbovePoint(heightAbovePoint,vertices[imapObject][triangle[0]], vertices[imapObject][triangle[1]],vertices[imapObject][triangle[2]],p , tempVertCoordMatrix, tempVertCoords + 3); //get the height of the triangle above p (considering that the point is projected vertically in the triangle)
          //cerr << "Height: " << heightAbovePoint << endl;

          //if the height of the point "above" (actually, above or below) is greater than the height of p --> the point is above p
          if(heightAbovePoint > p[2]) {        
            //highestCellZToProcess =  min(highestCellZToProcess,uniformGrid->z_cell_from_coord(heightAbovePoint, tempVertCoord,tempBigInts) );
            
            if (!foundATriangleAboveP || heightAbovePoint < heightAbovePointBestTriangle) {
              int czGlobal = uniformGrid->z_global_cell_from_coord(heightAbovePoint, tempVertCoord,tempBigInts);
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

            } else if (heightAbovePoint == heightAbovePointBestTriangle) {
            	//we need SoS to correctly break ties here...
            	//cerr << "Same height..." << endl;
            	bestTriangle = getBestTrianglePointInObjectSoS(imapObject,&triangle,bestTriangle,p);
            }            
          } else if(heightAbovePoint == p[2])  {
          	//the point is on the triangle... we need SoS
          		//cerr << "On the triangle" << endl;
          	if (triangleAbovePointSoS(imapObject,triangle,p)) {
          	//	cerr << "Updating..." << endl;
	        		highestCellZToProcess  = uniformGrid->z_cell_from_coord_level1(heightAbovePoint, tempVertCoord,tempBigInts) ;
	            foundATriangleAboveP = true;
	            bestTriangle = &triangle;
	            heightAbovePointBestTriangle = heightAbovePoint;
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

 // cerr << "Best height: " << heightAbovePointBestTriangle << endl;

  //cerr << "Found a triangle? " << foundATriangleAboveP << endl;
  if (foundATriangleAboveP) {
    //now we need to, basing on the orientation of the triangles vertices, find what is the object below the triangle
    //if the angle between the vector from p (point in the z+ direction) and the normal of the triangle is < 90 --> the point is below the triangle
    //if the angle is 90 --> error! the triangle is vertical (we should treat this... TODO)
    const Triangle &triangle = *bestTriangle; //this is the first triangle directly above p...

    const Point &p0 = vertices[imapObject][triangle[0]];
    const Point &p1 = vertices[imapObject][triangle[1]];
    const Point &p2 = vertices[imapObject][triangle[2]];

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
      cerr << "Error... vertical triangle (TODO: fix this)" << endl;
      exit(1);
    }
    if (C<0) {
      return triangle.above;
    } else {
      return triangle.below;
    }
  }
  else return OUTSIDE_OBJECT; //there is no triangle above p ... p must not be inside any object from the other map
}








#include "floodFillScanline.cpp"



//Locate a set of vertices in the objects in map 0...
void locateVerticesInObject(const Nested3DGridWrapper *uniformGrid, int gridSize,int nestedGridSize, const vector<Point> &verticesToLocate,std::vector<ObjectId> &verticesIds) {
  timespec t0,t1,t01;

  int sum = 0; //we use this because the compiler may remove the call if we do not use the result...


  int numVerticesToLocate = verticesToLocate.size();

  
  cerr << "Computing CCs of empty boxes..." << endl;
  clock_gettime(CLOCK_REALTIME, &t0);


  //we will compute where the empty cells are (to accelerate the point in vol computations!)

  GridCellsLabels cellsLabels(gridSize); 

  big_int tempVarsInt[2];
  VertCoord tempVertCoordMatrix[2][3];
  VertCoord tempVertCoords[8];
  VertCoord tempVar;

  Point box0 = uniformGrid->box[0];
  VertCoord cellWidth = uniformGrid->cellWidthLevel1 ;
  VertCoord nestedGridCellWidth  = uniformGrid->cellWidthLevel2 ;
  VertCoord halfCellWidth = cellWidth/2;
  VertCoord halfCellWidthNestedGrid = nestedGridCellWidth/2;
  

  //const vector<vector<vector<Nested3DGridCell> > >  &gridCells = uniformGrid->grid.gridCells;

  const Nested3DGrid &gridFirstLevel = (uniformGrid->grid);
  for(int gx=0;gx<gridSize;gx++) 
    for(int gy=0;gy<gridSize;gy++)
      for(int gz=0;gz<gridSize;gz++) {
        //const vector<Triangle *> &trianglesInCell = gridCells[gx][gy][gz].triangles[0];
        Triangle * *trianglesInCell = gridFirstLevel.getPointerStartListTriangles(0,gridSize,gx,gy,gz);
        int numTrianglesInCell = gridFirstLevel.numTrianglesInGridCell(0,gridSize,gx,gy,gz);

        if(gridFirstLevel.hasSecondLevel(gx,gy,gz)) { //we have a nested grid...
            const Nested3DGrid &nestedGrid = *gridFirstLevel.childGrids[gx][gy][gz];//*gridCells[gx][gy][gz].childGrid;

            cellsLabels.childGridLabels[gx][gy][gz] = new GridCellsLabels(nestedGridSize);

            vector<vector<vector<int > > > &labels = cellsLabels.childGridLabels[gx][gy][gz]->labels;

            for(int nx=0;nx<nestedGridSize;nx++)  //lets initialize the nested grid labels...
              for(int ny=0;ny<nestedGridSize;ny++)
                for(int nz=0;nz<nestedGridSize;nz++)
                  if (nestedGrid.numTrianglesInGridCell(0, nestedGridSize,nx,ny,nz) !=0)//if (nestedGrid.gridCells[nx][ny][nz].triangles[0].size()!=0)
                    labels[nx][ny][nz] = -3;

        } else { //this grid cell does not have a nested grid... and there is at least one triangle in this grid cell...
          if(numTrianglesInCell !=0) //if (trianglesInCell.size()!=0)
            cellsLabels.labels[gx][gy][gz] = -3;
        }        
      }
 // cerr << "End of filling labels.." << endl;
  clock_gettime(CLOCK_REALTIME, &t01);  



  int numSearchsPerformedToFillGrid = 0;
  

  //#define DONT_USE_GRID_ACCELERATION

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
              int temp;
              ObjectId id = computeObjectWherePointIsTwoLevel(p,gx*nestedGridSize,gy*nestedGridSize,gz*nestedGridSize ,0, uniformGrid, tempVertCoords, tempVertCoordMatrix, tempVarsInt,cellsLabels,foundUsingGrid,temp,temp);
              //cerr << "setiing to " << id <<  endl;
              setObjectInWhereEmptyCellIs2(gx,gy,gz,-1,-1,-1,gridSize,nestedGridSize,cellsLabels,id,start[0],start[1],start[2],end[0],end[1],end[2]);

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
                    int temp;
                    ObjectId id = computeObjectWherePointIsTwoLevel(p,gx*nestedGridSize + nx,gy*nestedGridSize + ny,gz*nestedGridSize +nz,0, uniformGrid, tempVertCoords, tempVertCoordMatrix, tempVarsInt,cellsLabels,foundUsingGrid,temp,temp);
                   // cerr << "setiing 2 to " << id <<  endl;
                    setObjectInWhereEmptyCellIs2(gx,gy,gz,nx,ny,nz,gridSize,nestedGridSize,cellsLabels,id,start[0],start[1],start[2],end[0],end[1],end[2]);
                   // cerr << "set level 2" << endl;
        

                    numSearchsPerformedToFillGrid++;
                  }
            //cerr << "End grid nested..." << endl;
          }                             
        } 
  }




  #endif 

//  cerr << "Locating two level......" << endl;
  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to fill grid and compute CCs: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  cerr << "Time only to compute CCs: " << convertTimeMsecs(diff(t01,t1))/1000 << endl;


  int numEmptyGridCells=0,numNonEmptyGridCells= 0;
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
  cerr << "Percent empty cells: " << (100.0*numEmptyGridCells)/(numEmptyGridCells+numNonEmptyGridCells) << endl; 

  /*clock_gettime(CLOCK_REALTIME, &t0);
  vector<array<int,3> > pointsGridCells(numVerticesInMap); //this will store in what grid cell each point is..
  #pragma omp parallel
  {
    big_int tempVarsInt[2];
    VertCoord tempVar;

    #pragma omp for schedule(static,32)
    for(int i=0;i<numVerticesInMap;i++) {
      const Point &p = verticesToLocate[i];  
      pointsGridCells[i][0] = uniformGrid->x_cell_from_coord(p[0], tempVar,tempVarsInt);
      pointsGridCells[i][1] = uniformGrid->y_cell_from_coord(p[1], tempVar,tempVarsInt);
      pointsGridCells[i][2] = uniformGrid->z_cell_from_coord(p[2], tempVar,tempVarsInt);
    }
  }
  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to compute in what grid cell each point is: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;*/


  vector<int> gxVector(numVerticesToLocate);
  vector<int> gyVector(numVerticesToLocate);
  vector<int> gzVector(numVerticesToLocate);

  cerr << "Computing vertices  " << verticesToLocate.size() << " grid coordinates..." << endl;
  clock_gettime(CLOCK_REALTIME, &t0);

  #pragma omp parallel
  {
    big_int tempVarsInt[2];
    VertCoord tempVar;
    
    #pragma omp for schedule(dynamic,1000)
    for(int i=0;i<numVerticesToLocate;i++) {
      const Point &p = verticesToLocate[i];        
    
      gxVector[i] = uniformGrid->x_global_cell_from_coord(p[0], tempVar,tempVarsInt);
      gyVector[i] = uniformGrid->y_global_cell_from_coord(p[1], tempVar,tempVarsInt);
      gzVector[i] = uniformGrid->z_global_cell_from_coord(p[2], tempVar,tempVarsInt);
      
    }
  }
  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to compute vertices grid...: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;



  cerr << "Locating " << verticesToLocate.size() << " vertices..." << endl;
  clock_gettime(CLOCK_REALTIME, &t0);

  int numVerticesFoundUsingGrid = 0;
  int totalNumGridsCellToCheck = 0;
  int numTrianglesChecked = 0;
  //#pragma omp parallel
  {
    big_int tempVarsInt[2];
    VertCoord tempVertCoordMatrix[2][3];
    VertCoord tempVertCoords[8];
    VertCoord tempVar;
    int numVerticesFoundUsingGridLocal = 0;
    bool foundUsingGrid;
    
    //#pragma omp for schedule(dynamic,1000)
    for(int i=0;i<numVerticesToLocate;i++) {
      const Point &p = verticesToLocate[i];        
    
      //int gx = uniformGrid->x_global_cell_from_coord(p[0], tempVar,tempVarsInt);
      //int gy = uniformGrid->y_global_cell_from_coord(p[1], tempVar,tempVarsInt);
      //int gz = uniformGrid->z_global_cell_from_coord(p[2], tempVar,tempVarsInt);
      int gx = gxVector[i];
      int gy = gyVector[i];
      int gz = gzVector[i];

      int numGridsCellToCheck;
      ObjectId id2 = computeObjectWherePointIsTwoLevel(p,gx,gy,gz,0, uniformGrid, tempVertCoords, tempVertCoordMatrix, tempVarsInt,cellsLabels,foundUsingGrid, numGridsCellToCheck, numTrianglesChecked);
      numVerticesFoundUsingGridLocal += foundUsingGrid;

     	totalNumGridsCellToCheck += numGridsCellToCheck;

      verticesIds[i] = id2;
    }
    //#pragma omp critical
    numVerticesFoundUsingGrid += numVerticesFoundUsingGridLocal;
  }
  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to locate vertices: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  cerr << "Num vertices found using grid: " << numVerticesFoundUsingGrid << endl;
  cerr << "Percentage vertices found using grid: " << (100.0*numVerticesFoundUsingGrid)/numVerticesToLocate << "\n\n";
  cerr << "Total grid cells checked: " << totalNumGridsCellToCheck << endl;
  cerr << "Total num triangles checked: " << numTrianglesChecked << endl;

}









//==================================================================================================================================

void saveEdges(string path) {
  ofstream fout(path.c_str());

  int numVert = 3*edges.size();
  int numEdges = 3*edges.size();
  int numFaces = edges.size();
  fout << numVert << " " << numEdges << " " << numFaces << "\n";
  for(pair< array<VertCoord,3>,array<VertCoord,3> > e:edges) {
    fout << e.first[0].get_d() << " " << e.first[1].get_d() << " " << e.first[2].get_d() << "\n";
    fout << e.second[0].get_d() << " " << e.second[1].get_d() << " " << e.second[2].get_d() << "\n";
    fout << e.second[0].get_d() << " " << e.second[1].get_d() << " " << e.second[2].get_d() << "\n";
  }
  fout << endl;

  int start = 1;
  for(int i=0;i<numFaces;i++) {
    fout << start << " " << start+1 << "\n";
    fout << start+1 << " " << start+2 << "\n";
    fout << start+2 << " " << start << "\n";
    start+= 3;
  }
  fout << endl;
  start = 1;
  for(int i=0;i<numFaces;i++) {
    fout << start << " " << start+1 << " " << start+2 <<  "\n";
    start+= 3;
  }


}






// TODO: if(a*b > 0 ) --> do not multiply!!! use test a>>0 && b>0 || a<0 && b<0

int main(int argc, char **argv) {
  if (argc!=6) {
      cerr << "Error... use ./point_in_volume inputMeshes gridSizeLevel1 gridSizeLevel2 triggerSecondLevel maxTreeDepth" << endl;
      cerr << "The mesh file may be in the gts format or in the lium format (multimaterial)" << endl;
      exit(1);
  }

  string inputObject = argv[1];
  int gridSizeLevel1 = atoi(argv[2]);
  int gridSizeLevel2 = atoi(argv[3]);
  int triggerSecondLevel = atoi(argv[4]);
  int maxTreeDepth = atoi(argv[5]);


  cerr << "Reading input points..." << endl;
  vector<Point> pointsToLocate;
  Point p;
  double x,y,z;
  while(cin >> x >> y >> z) {
    p[0] = x;
    p[1] = y;
    p[2] = z;
    pointsToLocate.push_back(p);
  }
  cerr << "Reading objects..." << endl;

  timespec t0,t1;
  clock_gettime(CLOCK_REALTIME, &t0); 

  assert(inputObject.size()>=4);
  string inputFileExtension = inputObject.substr(inputObject.size()-3,3);
  assert(inputFileExtension=="gts" || inputFileExtension=="ium");
  bool isGtsFile = inputFileExtension=="gts";
  if(isGtsFile)
    TIME(readGTSFile(inputObject,vertices[0],triangles[0],bBox));
  else 
    TIME(readLiumFile(inputObject,vertices[0],triangles[0],bBox));

  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to read: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;

    cerr << "Bounding box map 1: " << setprecision(18) << std::fixed << bBox[0][0].get_d() << " " << bBox[0][1].get_d() <<  " " << bBox[0][2].get_d() << " -- " << bBox[1][0].get_d() << " " << bBox[1][1].get_d() << " " << bBox[1][2].get_d() <<endl;


  clock_gettime(CLOCK_REALTIME, &t0); 
  cerr <<"Creating nested grid..." << endl;
  vector<Triangle *> trianglesPointers[2];

  int sz = triangles[0].size();
  trianglesPointers[0].resize(sz);
  for(int i=0;i<sz;i++) trianglesPointers[0][i]  = & (triangles[0][i]);


//  void init(const vector<Triangle *> trianglesInsert[2], const vector<Point> vertices[2], const int gridSizeLevel1, const int gridSizeLevel2, const Point &p0, const Point &p1,const long long prodThreshold);

  TIME(uniformGrid.init(trianglesPointers, vertices, gridSizeLevel1,gridSizeLevel2,bBox[0],bBox[1],triggerSecondLevel));


  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to create and refine grid: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;







  //locateVerticesOfOneObjectInTheOther(&uniformGrid,cellWidth, nestedGridCellWidth,gridSizeLevel1,gridSizeLevel2);

  vector<ObjectId> pointIds(pointsToLocate.size());
  locateVerticesInObject(&uniformGrid,gridSizeLevel1,gridSizeLevel2,pointsToLocate,pointIds);

  if(isGtsFile)
    for(const ObjectId& id:pointIds) {
      cout << id << "\n";
    }
  else
    for(const ObjectId& id:pointIds) {  //when we deal with LIUM files the -1 represents the outside object and objects start with id 0.... (we convert from -1 to 0 and from 0 to 1 in the input... so, we need to revert this change...s)
    
        cout << id-1 << "\n";
    } 

}