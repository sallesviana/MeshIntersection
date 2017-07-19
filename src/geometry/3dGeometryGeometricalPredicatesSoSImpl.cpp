
#include "sosPredicatesImpl.h"
#include "originalAlgFromMathematicaSosPredicatesImpl.h"


inline int getSignal(int i){
  if(i>0) return 1;
  if(i<0) return -1;
  return 0;
}

//#define SosPredicatesImpl OriginalAlgFromMathematicaSosPredicatesImpl

int MeshIntersectionGeometry::orientation(const InputVertex &v1, const InputVertex &v2, const InputVertex &queryPoint,int whatPlaneProjectTrianglesTo,TempVarsSoSPredicatesImpl &tempVars)  { 
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

  //int ans = sgn( (p1[coordX]-p0[coordX])*(p[coordY]-p0[coordY]) -  (p1[coordY]-p0[coordY])*(p[coordX]-p0[coordX]) );
  tempVars.tmp = p1[coordX];
  tempVars.tmp -= p0[coordX];
  tempVars.tmp2 = p[coordY];
  tempVars.tmp2 -= p0[coordY];
  tempVars.tmp *= tempVars.tmp2;
  tempVars.tmp2 = p1[coordY];
  tempVars.tmp2 -= p0[coordY];
  tempVars.tmp3 = p[coordX];
  tempVars.tmp3 -= p0[coordX];
  tempVars.tmp2 *= tempVars.tmp3;
  int ans = getSignal(cmp(tempVars.tmp,tempVars.tmp2));

  if(ans!=0) return ans;

  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.orientation2DOOO++;
  #endif  




  //a.x * b.y - a.y * b.x;
  int ans2 = SosPredicatesImpl(this,tempVars).orientation2D(v1,v2,queryPoint,whatPlaneProjectTrianglesTo);
  #ifdef DOUBLE_CHECK_SOS_PREDICATES_WITH_MATHEMATICA
  	assert(ans2==OriginalAlgFromMathematicaSosPredicatesImpl(this).orientation2D(v1,v2,queryPoint,whatPlaneProjectTrianglesTo) );
  #endif

  assert(ans==0 || ans==ans2);

  return ans2;
}

int MeshIntersectionGeometry::orientation(const InputVertex &v1,const InputVertex &v2, const VertexFromIntersection &queryPoint, int whatPlaneProjectTrianglesTo,TempVarsSoSPredicatesImpl &tempVars)  { 
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
  //int ans = sgn( (p1[coordX]-p0[coordX])*(p[coordY]-p0[coordY]) -  (p1[coordY]-p0[coordY])*(p[coordX]-p0[coordX]) );
  tempVars.tmp = p1[coordX];
  tempVars.tmp -= p0[coordX];
  tempVars.tmp2 = p[coordY];
  tempVars.tmp2 -= p0[coordY];
  tempVars.tmp *= tempVars.tmp2;
  tempVars.tmp2 = p1[coordY];
  tempVars.tmp2 -= p0[coordY];
  tempVars.tmp3 = p[coordX];
  tempVars.tmp3 -= p0[coordX];
  tempVars.tmp2 *= tempVars.tmp3;
  int ans = getSignal(cmp(tempVars.tmp,tempVars.tmp2));

  if(ans!=0) return ans;


  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.orientation2DOOI++;
  #endif  
  



  int ans2 = SosPredicatesImpl(this,tempVars).orientation2D(v1,v2,queryPoint,whatPlaneProjectTrianglesTo);
  #ifdef DOUBLE_CHECK_SOS_PREDICATES_WITH_MATHEMATICA
  	assert(ans2==OriginalAlgFromMathematicaSosPredicatesImpl(this).orientation2D(v1,v2,queryPoint,whatPlaneProjectTrianglesTo));
  #endif

  assert(ans==0 || ans==ans2);
    

  return ans2;
}

int MeshIntersectionGeometry::orientation(const InputVertex &v1, const VertexFromIntersection &v2, const VertexFromIntersection &queryPoint, int whatPlaneProjectTrianglesTo,TempVarsSoSPredicatesImpl &tempVars) { 
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

  //int ans = sgn( (p1[coordX]-p0[coordX])*(p[coordY]-p0[coordY]) -  (p1[coordY]-p0[coordY])*(p[coordX]-p0[coordX]) );
  tempVars.tmp = p1[coordX];
  tempVars.tmp -= p0[coordX];
  tempVars.tmp2 = p[coordY];
  tempVars.tmp2 -= p0[coordY];
  tempVars.tmp *= tempVars.tmp2;
  tempVars.tmp2 = p1[coordY];
  tempVars.tmp2 -= p0[coordY];
  tempVars.tmp3 = p[coordX];
  tempVars.tmp3 -= p0[coordX];
  tempVars.tmp2 *= tempVars.tmp3;
  int ans = getSignal(cmp(tempVars.tmp,tempVars.tmp2));

  if(ans!=0) return ans;


  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.orientation2DOII++;
  #endif  


  //a.x * b.y - a.y * b.x;  
  int ans2 = SosPredicatesImpl(this,tempVars).orientation2D(v1,v2,queryPoint,whatPlaneProjectTrianglesTo);
 	#ifdef DOUBLE_CHECK_SOS_PREDICATES_WITH_MATHEMATICA
  	assert(ans2==OriginalAlgFromMathematicaSosPredicatesImpl(this).orientation2D(v1,v2,queryPoint,whatPlaneProjectTrianglesTo));
  #endif

  assert(ans==0 || ans==ans2);


  return ans2;
}

int MeshIntersectionGeometry::orientation(const VertexFromIntersection &v1, const VertexFromIntersection &v2, const VertexFromIntersection &queryPoint, int whatPlaneProjectTrianglesTo,TempVarsSoSPredicatesImpl &tempVars) { 
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
  //int ans = sgn( (p1[coordX]-p0[coordX])*(p[coordY]-p0[coordY]) -  (p1[coordY]-p0[coordY])*(p[coordX]-p0[coordX]) );
  tempVars.tmp = p1[coordX];
  tempVars.tmp -= p0[coordX];
  tempVars.tmp2 = p[coordY];
  tempVars.tmp2 -= p0[coordY];
  tempVars.tmp *= tempVars.tmp2;
  tempVars.tmp2 = p1[coordY];
  tempVars.tmp2 -= p0[coordY];
  tempVars.tmp3 = p[coordX];
  tempVars.tmp3 -= p0[coordX];
  tempVars.tmp2 *= tempVars.tmp3;
  int ans = getSignal(cmp(tempVars.tmp,tempVars.tmp2));

  if(ans!=0) return ans;
  
  //I am not sure if this is true... double check..()
  /*if(v1.getMeshOfTriangleDefiningVertex() == v2.getMeshOfTriangleDefiningVertex() && v1.getMeshOfTriangleDefiningVertex()==queryPoint.getMeshOfTriangleDefiningVertex()) {
  	//if the points are generated from intersections with the same triangle, --> SoS will not modify the answer of the orientation!
  	if(v1.triangle.compare(v2.triangle)==0 && v2.triangle.compare(queryPoint.triangle)==0) {
  		#ifdef DOUBLE_CHECK_SOS_RESULTS
	  		int ans2 = SosPredicatesImpl(this).orientation2D(v1,v2,queryPoint,whatPlaneProjectTrianglesTo);
	  		assert(ans==ans2);
	  	#endif
	  	#ifdef DOUBLE_CHECK_SOS_PREDICATES_WITH_MATHEMATICA
		  	assert(SosPredicatesImpl(this).orientation2D(v1,v2,queryPoint,whatPlaneProjectTrianglesTo)==OriginalAlgFromMathematicaSosPredicatesImpl(this).orientation2D(v1,v2,queryPoint,whatPlaneProjectTrianglesTo));
		  #endif	

	  	return ans;
	  }
  }*/


  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.orientation2DIII++;
  #endif
  
  int ans2 = SosPredicatesImpl(this,tempVars).orientation2D(v1,v2,queryPoint,whatPlaneProjectTrianglesTo);
  #ifdef DOUBLE_CHECK_SOS_PREDICATES_WITH_MATHEMATICA
  	assert(ans2==OriginalAlgFromMathematicaSosPredicatesImpl(this).orientation2D(v1,v2,queryPoint,whatPlaneProjectTrianglesTo));
  #endif
      
  return ans2;
}

int cts[20] = {0};
//TODO: remove cts++ for performance...

//1 if make a left turn, -1 if make a right turn, 0 --> degeneracy (shouldn't happen...)
int MeshIntersectionGeometry::orientation(const Vertex &v1, const Vertex &v2, const Vertex &p, int whatPlaneProjectTrianglesTo,TempVarsSoSPredicatesImpl &tempVars)  { 
  bool isV1InputVertex = (&v1)->isInputVertex();
  bool isV2InputVertex = (&v2)->isInputVertex();  
  bool isPInputVertex = (&p)->isInputVertex();

  int numInputVertices = isV1InputVertex+isV2InputVertex+isPInputVertex;
  if(numInputVertices==3) {
    cts[0]++;
    return orientation(*static_cast<const InputVertex*>(&v1),*static_cast<const InputVertex*>(&v2),*static_cast<const InputVertex*>(&p),whatPlaneProjectTrianglesTo,tempVars);
  } else if(numInputVertices==0) {
    cts[1]++;
    return orientation(*static_cast<const VertexFromIntersection*>(&v1),*static_cast<const VertexFromIntersection*>(&v2),*static_cast<const VertexFromIntersection*>(&p),whatPlaneProjectTrianglesTo,tempVars);
  } else if(numInputVertices==1) {
    if(isV1InputVertex) {
      cts[2]++;
      return orientation(*static_cast<const InputVertex*>(&v1),*static_cast<const VertexFromIntersection*>(&v2),*static_cast<const VertexFromIntersection*>(&p),whatPlaneProjectTrianglesTo,tempVars);
    } else if(isV2InputVertex) {
      cts[3]++;
      return -orientation(*static_cast<const InputVertex*>(&v2),*static_cast<const VertexFromIntersection*>(&v1),*static_cast<const VertexFromIntersection*>(&p),whatPlaneProjectTrianglesTo,tempVars);
    } else { //p is the input vertex...
      cts[4]++;
      return -orientation(*static_cast<const InputVertex*>(&p),*static_cast<const VertexFromIntersection*>(&v2),*static_cast<const VertexFromIntersection*>(&v1),whatPlaneProjectTrianglesTo,tempVars);
    }
  } else { //numInputVertices is 2...
    if(!isV1InputVertex) { // v2 and p are input vertices
      cts[5]++;
      return orientation(*static_cast<const InputVertex*>(&v2),*static_cast<const InputVertex*>(&p),*static_cast<const VertexFromIntersection*>(&v1),whatPlaneProjectTrianglesTo,tempVars);
    } else if(!isV2InputVertex) { // v1 and p are input vertices
      cts[6]++;
      return orientation(*static_cast<const InputVertex*>(&p),*static_cast<const InputVertex*>(&v1),*static_cast<const VertexFromIntersection*>(&v2),whatPlaneProjectTrianglesTo,tempVars);
    } else { // v1 and v2 are input vertices...
      cts[7]++;
      return orientation(*static_cast<const InputVertex*>(&v1),*static_cast<const InputVertex*>(&v2),*static_cast<const VertexFromIntersection*>(&p),whatPlaneProjectTrianglesTo,tempVars);
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


/*
int signDeterminant4(const Point &p1,const Point &p2,const Point &p3,const Point &p4) {
  const VertCoord &p1x = p1[0];
  const VertCoord &p1y = p1[1];
  const VertCoord &p1z = p1[2];

  const VertCoord &p2x = p2[0];
  const VertCoord &p2y = p2[1];
  const VertCoord &p2z = p2[2];

  const VertCoord &p3x = p3[0];
  const VertCoord &p3y = p3[1];
  const VertCoord &p3z = p3[2];

  const VertCoord &p4x = p4[0];
  const VertCoord &p4y = p4[1];
  const VertCoord &p4z = p4[2];

  const VertCoord det = -p1z*p2y*p3x + p1y*p2z*p3x + p1z*p2x*p3y - p1x*p2z*p3y - p1y*p2x*p3z +
                         p1x*p2y*p3z + p1z*p2y*p4x - p1y*p2z*p4x - p1z*p3y*p4x + p2z*p3y*p4x + 
                         p1y*p3z*p4x - p2y*p3z*p4x - p1z*p2x*p4y + p1x*p2z*p4y + p1z*p3x*p4y - 
                         p2z*p3x*p4y - p1x*p3z*p4y + p2x*p3z*p4y + p1y*p2x*p4z - p1x*p2y*p4z - 
                         p1y*p3x*p4z + p2y*p3x*p4z + p1x*p3y*p4z - p2x*p3y*p4z;
  return sgn(det);
}*/


//3d orientation...
int MeshIntersectionGeometry::orientation(const InputVertex&p1, const InputVertex&p2,const InputVertex&p3, const InputVertex &v,TempVarsSoSPredicatesImpl &tempVars) {
  int ans = 0;//signDeterminant4(getCoordinates(p1),getCoordinates(p2),getCoordinates(p3),getCoordinates(v));
  if(ans!=0) return ans;


  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.orientation3DOO++;
  #endif

  int ans2 = SosPredicatesImpl(this,tempVars).orientation3D(p1,p2,p3,v);    
  #ifdef DOUBLE_CHECK_SOS_PREDICATES_WITH_MATHEMATICA
  	assert(ans2==OriginalAlgFromMathematicaSosPredicatesImpl(this).orientation3D(p1,p2,p3,v));
  #endif

  assert(ans==0 || ans==ans2);

  return ans2;
}

int MeshIntersectionGeometry::orientation(const InputVertex&p1, const InputVertex&p2,const InputVertex&p3, const VertexFromIntersection &v,TempVarsSoSPredicatesImpl &tempVars)  {
  int ans = 0;// signDeterminant4(getCoordinates(p1),getCoordinates(p2),getCoordinates(p3),getCoordinates(v));
  if(ans!=0) return ans;


  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.orientation3DOI++;
  #endif  

  int ans2 = SosPredicatesImpl(this,tempVars).orientation3D(p1,p2,p3,v);
  #ifdef DOUBLE_CHECK_SOS_PREDICATES_WITH_MATHEMATICA
  	assert(ans2==OriginalAlgFromMathematicaSosPredicatesImpl(this).orientation3D(p1,p2,p3,v));
  #endif
  
  assert(ans==0 || ans==ans2);

  return ans2;
}


//3d orientation.
int MeshIntersectionGeometry::orientation(const InputTriangle&t, const InputVertex &v,TempVarsSoSPredicatesImpl &tempVars)  {
  return orientation(*(t.getInputVertex(0)),*(t.getInputVertex(1)),*(t.getInputVertex(2)),v,tempVars);
}

int MeshIntersectionGeometry::orientation(const InputTriangle&t, const VertexFromIntersection &v,TempVarsSoSPredicatesImpl &tempVars)  {
  return orientation(*(t.getInputVertex(0)),*(t.getInputVertex(1)),*(t.getInputVertex(2)),v,tempVars);
}  








int MeshIntersectionGeometry::signalVectorCoordOnlyCallWhenCoincident(const InputVertex &orig, const InputVertex &dest, int coord,TempVarsSoSPredicatesImpl &tempVars)  {
  //signal vector coord(orig,dest) = -orientation(orig,dest)
  //
  //
  /*
  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.orientation1DOO++;
  #endif*/

  int ans = 0;
  int meshId0 = orig.getMeshId();
  int meshId1 = dest.getMeshId();
  if(meshId0!=meshId1) {
    if(meshId0==0) ans = 1; // orig is in mesh 0 --> non perturbed , dest is in 1--> positive perturb, dest-orig = positive perturb.
    else ans = -1;
  }

  #ifdef DOUBLE_CHECK_SOS_RESULTS    
  	int ans2 = -SosPredicatesImpl(this,tempVars).orientation1D(orig,dest,coord);
    assert(ans==ans2);
  #endif
  #ifdef DOUBLE_CHECK_SOS_PREDICATES_WITH_MATHEMATICA
  	assert(SosPredicatesImpl(this,tempVars).orientation1D(orig,dest,coord)==OriginalAlgFromMathematicaSosPredicatesImpl(this).orientation1D(orig,dest,coord));
  #endif

  return ans;
 
}


int MeshIntersectionGeometry::signalVectorCoordOnlyCallWhenCoincident(const InputVertex &orig, const VertexFromIntersection &dest, int coord,TempVarsSoSPredicatesImpl &tempVars)  {
  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.orientation1DOI++;
  #endif

  #ifdef DOUBLE_CHECK_SOS_PREDICATES_WITH_MATHEMATICA
  	assert(SosPredicatesImpl(this,tempVars).orientation1D(orig,dest,coord)==OriginalAlgFromMathematicaSosPredicatesImpl(this).orientation1D(orig,dest,coord));
  #endif  

  return -SosPredicatesImpl(this,tempVars).orientation1D(orig,dest,coord);
}

int MeshIntersectionGeometry::signalVectorCoordOnlyCallWhenCoincident(const VertexFromIntersection &orig, const VertexFromIntersection &dest, int coord,TempVarsSoSPredicatesImpl &tempVars)  {
  int ans = 0;

  /*
	//this does not seem to be true!
	//imagine the tip of a triangle touching another triangle... after perturbation the signal can become non zero, right?
  if(orig.getMeshOfTriangleDefiningVertex() == dest.getMeshOfTriangleDefiningVertex() ) { //if they are generated from triangles from same mesh --> if there is a coincidence it is really a coincidence
  	#ifdef DOUBLE_CHECK_SOS_RESULTS
  		int ans2 = -SosPredicatesImpl(this).orientation1D(orig,dest,coord);
  		assert(ans2==ans);
  	#endif

  	return ans;
  }*/

  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.orientation1DII++;
  #endif  

  int ans2 = -SosPredicatesImpl(this,tempVars).orientation1D(orig,dest,coord);

  #ifdef DOUBLE_CHECK_SOS_PREDICATES_WITH_MATHEMATICA
  	assert(ans2==-OriginalAlgFromMathematicaSosPredicatesImpl(this).orientation1D(orig,dest,coord));
  #endif

  return ans2; 
}


//this is basically what we need for the "brute force" retesselation...
//this is essentially a 1D orientation!!!!
//TODO: implement this as 1D orientation...
//what is the signal of each coordinate the vector from orig to dest
//can be 0 (SoS)
//I think this could be 0... suppose the two vertices are from same mesh, for example...
//TODO: review...
int MeshIntersectionGeometry::signalVectorCoord(const Vertex &orig, const Vertex &dest, int coord,TempVarsSoSPredicatesImpl &tempVars)  {
  const Point &p0 =  getCoordinates(orig);
  const Point &p1 =  getCoordinates(dest);
  int ans = getSignal(cmp(p1[coord],p0[coord])); // sgn(p1[coord]-p0[coord]);
  //if(p1[coord] > p0[coord]) ans = 1;
  //if(p1[coord] < p0[coord]) ans = -1;

  if(ans!=0) return ans;

  
  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.signVector++;
  #endif

  bool isV1InputVertex = (&orig)->isInputVertex();
  bool isV2InputVertex = (&dest)->isInputVertex();  

  if(isV1InputVertex) {
    if(isV2InputVertex) //input,input
      return signalVectorCoordOnlyCallWhenCoincident(*static_cast<const InputVertex*>(&orig),*static_cast<const InputVertex*>(&dest),coord,tempVars); 
    else //input, intersection
      return signalVectorCoordOnlyCallWhenCoincident(*static_cast<const InputVertex*>(&orig),*static_cast<const VertexFromIntersection*>(&dest),coord,tempVars); 
  }
  else { 
    if(isV2InputVertex)
      return -signalVectorCoordOnlyCallWhenCoincident(*static_cast<const InputVertex*>(&dest),*static_cast<const VertexFromIntersection*>(&orig),coord,tempVars);
    else
      return signalVectorCoordOnlyCallWhenCoincident(*static_cast<const VertexFromIntersection*>(&orig),*static_cast<const VertexFromIntersection*>(&dest),coord,tempVars); 
  }
}

/*
//TODO: review this...
int MeshIntersectionGeometry::signalVectorCoordCanBe0(const Vertex &orig, const Vertex &dest, int coord) const {
  #ifdef COLLECT_GEOMETRY_STATISTICS
    #pragma omp atomic
    geometryStatisticsDegenerateCases.signVector++;
  #endif

  const Point &p0 =  getCoordinates(orig);
  const Point &p1 =  getCoordinates(dest);
  int ans = sgn(p1[coord]-p0[coord]);
  return  ans;
}
*/

/*
Predicates and functions using SoS
*/



//Is v1 closer to origV than v2 is?
bool MeshIntersectionGeometry::isCloserSoSImpl(const InputVertex &orig, const VertexFromIntersection &v1, const VertexFromIntersection &v2, TempVarsIsCloser &tempVars)  {
  //we know that v1 and v2 are formed by the intersection of v1.triangle and an edge with endpoint in orig ; and v2.triangle and an
  //edge with endpoint in orig. 
  //Thus, v1 will be closer to orig iff v1 and orig are on the same side of v2.triangle 
  return orientation(v2.triangle,orig,tempVars.tempVarsSoSPredicatesImpl) == orientation(v2.triangle,v1,tempVars.tempVarsSoSPredicatesImpl);
}


bool MeshIntersectionGeometry::isVertexInTriangleProjectionSoSImpl(const Vertex &v1,const Vertex &v2, const Vertex &v3, const Vertex &queryPoint,int whatPlaneProjectTrianglesTo,TempVarsIsVertexTriangleProjection &tempVars)  {
  int o1 = orientation(v1,v2,queryPoint,whatPlaneProjectTrianglesTo,tempVars.tempVarsSoSPredicatesImpl);
  int o2 = orientation(v2,v3,queryPoint,whatPlaneProjectTrianglesTo,tempVars.tempVarsSoSPredicatesImpl);
  int o3 = orientation(v3,v1,queryPoint,whatPlaneProjectTrianglesTo,tempVars.tempVarsSoSPredicatesImpl);
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


bool MeshIntersectionGeometry::isVertexConvexSoSImpl(const Vertex &v1,const Vertex &queryVertex, const Vertex &v3,int whatPlaneProjectTrianglesTo,TempVarsIsVertexConvex &tempVars)  {
  return orientation(v1,queryVertex,v3,whatPlaneProjectTrianglesTo,tempVars.tempVarsSoSPredicatesImpl)<0;
}


bool MeshIntersectionGeometry::isVertexInTriangleProjectionSoSImpl(const InputTriangle &t, const InputVertex &queryPoint,TempVarsIsVertexTriangleProjectionZ0 &tempVars)  {
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

  int o1 = orientation(v1,v2,queryPoint,whatPlaneProjectTrianglesTo,tempVars.tempVarsSoSPredicatesImpl);
  int o2 = orientation(v2,v3,queryPoint,whatPlaneProjectTrianglesTo,tempVars.tempVarsSoSPredicatesImpl);
  int o3 = orientation(v3,v1,queryPoint,whatPlaneProjectTrianglesTo,tempVars.tempVarsSoSPredicatesImpl);
  if(o2!=0 && o3!=0 && o2!=o3)
    return false;
  if(o1!=0 && o2!=0 && o1!=o2)
    return false;
  if(o1!=0 && o3!=0 && o1!=o3)
    return false;
  return true;
}

//The input triangle should not be vertical (PinMesh does not use vertical triangles because of the perturbation!)
//Given a vertex p, is p below the triangle t ? (we know p projected to z=0 is on t projected to z=0...)
bool MeshIntersectionGeometry::isTriangleAbovePointSoSImpl(const InputTriangle &t, const InputVertex &p,TempVarIsTriangleAbovePointSoS &tempVars)  {
  //The orientation(t,p) < 0 iff the point is on the side pointed by the normal
  int sideOfTriangle = -orientation(t,p,tempVars.tempVarsSoSPredicatesImpl); //is p on the positive (1) or negative (-1) side of the triangle?

 // TempVarIsTriangleNormalPointingPositiveZ temp;

  //cerr << "Which side? " << sideOfTriangle << endl;
  //if the point is on the positive side, it will be below the triangle if it points down
  if(sideOfTriangle==1) return !isTriangleNormalPointingPositiveZSoSImpl(t,tempVars.tempVarsTriangleNormal);
  else return isTriangleNormalPointingPositiveZSoSImpl(t,tempVars.tempVarsTriangleNormal);
  //if the point is on the negative sie, it will be below the triangle if it points up...
}



bool MeshIntersectionGeometry::isOrientationPositiveSoSImpl(const Vertex &origV, const Vertex &v1V, const Vertex &v2V, const int planeToProject, TempVarsIsAngleWith0Greater &tempVars)  {
  return orientation(origV,v1V,v2V,planeToProject,tempVars.tempVarsSoSPredicatesImpl)>0;
}



int getCodeVertexInQuadrant(int sgnX,int sgnY) {
  //Let's suppose: 0 is x+ axis, 1 is x+,y+ , 2 is y+ axis, 3 is x-,y+, 4 is x-, 5 is ...
  //      2
  //    3 | 1
  //      |
  //  4-------0
  //      |
  //    5 | 7
  //      6

  assert(!(sgnX==0 && sgnY==0));
  //Given the signals of a vector, return its code... Indefinite result for 0,0 (shouldn't happen...)
  if(sgnY==0) {
    if(sgnX>0) return 0;
    else return 4;
  } else if(sgnY>0) { //above x axis...
    if(sgnX==0) return 2;
    else if(sgnX>0) return 1;
    else return 3;
  } else { //sgnY is < 0
    if(sgnX==0) return 6;
    else if(sgnX>0) return 7;
    else return 5;
  }
  
}


//TODO
//we need to consider point on axis and the angle...
//this function have to be complete... we will call it, for example, when we have degenerate edges.
bool MeshIntersectionGeometry::isAngleWith0GreaterSoSImpl(const Vertex &origV, const Vertex &v1V, const Vertex &v2V, const int planeToProject, TempVarsIsAngleWith0Greater &tempVars)  {
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

  const int sgnV1y = signalVectorCoord(origV,v1V,yCoord,tempVars.tempVarsSoSPredicatesImpl);
  const int sgnV2y = signalVectorCoord(origV,v2V,yCoord,tempVars.tempVarsSoSPredicatesImpl);

  //TODO: consider the case where one of the vectors is exactly on 0... this can happen!
  
  if(sgnV1y>0 && sgnV2y<0) return true; //check if the two vectors are in different sides of the x axis...
  if(sgnV1y<0 && sgnV2y>0) return false;

  //they are in the same side of the x axis... (or exactly on the x axis...)
  // ---> No!!!! ONE may be on the x axis and the other one in any side...
  const int sgnV1x = signalVectorCoord(origV,v1V,xCoord,tempVars.tempVarsSoSPredicatesImpl);
  const int sgnV2x = signalVectorCoord(origV,v2V,xCoord,tempVars.tempVarsSoSPredicatesImpl);

  //Let's suppose: 0 is x+ axis, 1 is x+,y+ , 2 is y+ axis, 3 is x-,y+, 4 is x-, 5 is ...
  //      2
  //    3 | 1
  //      |
  //  4-------0
  //      |
  //    5 | 7
  //      6
  int codeV1 = getCodeVertexInQuadrant(sgnV1x,sgnV1y);
  int codeV2 = getCodeVertexInQuadrant(sgnV2x,sgnV2y);
  if(codeV1!=codeV2) return codeV1 < codeV2;
 /*
 Buggy code... int posV1 = 0;

  if(sgnV1y>0) { //are both at the non-negative y?
    //check if their x is different...
    if(sgnV1x > 0 && sgnV2x < 0) return true;
    if(sgnV1x < 0 && sgnV2x > 0) return false;
  } else if(sgnV1y<0) { //are they both at the negative y?
    if(sgnV1x > 0 && sgnV2x < 0) return false;
    if(sgnV1x < 0 && sgnV2x > 0) return true;
  }

  //is one of them on the x+ axis?
  if(sgnV1y==0 && sgnV1x>0) return true; // v1 is on x+ --> angle is 0
  if(sgnV2y==0 && sgnV2x>0) return false; //v2 is on x+ --> angle of v2 is 0 --> v2 < v1..*/

  
  //same side of the x-axis --> use vector orientation...
  return orientation(origV,v1V,v2V,planeToProject,tempVars.tempVarsSoSPredicatesImpl)>0;
}

//TODO: use pre-computed normals here...
//If all points from the same mesh are translated equally, we will not need SoS here
//because we will never chose an input triangle that is vertical 

//TO think: is it true that if we project the triangle to z=0, this function
//should return true if v3 is to the left of the vector v1-v2 ?
bool MeshIntersectionGeometry::isTriangleNormalPointingPositiveZSoSImpl(const InputTriangle &t, TempVarIsTriangleNormalPointingPositiveZ &tempVars)  {
    return orientation(*t.getInputVertex(0),*t.getInputVertex(1),*t.getInputVertex(2),PLANE_Z0,tempVars.tempVarsSoSPredicatesImpl)==1;
}




bool MeshIntersectionGeometry::intersectEdgeWithTriangleSoSImpl(const InputTriangle &triangle, const InputVertex &p1, const InputVertex &p2, TempVarsComputeIntersections &tempVars)  {
  int orientationP1Triangle = orientation(triangle,p1,tempVars.tempVarsSoSPredicatesImpl);
  int orientationP2Triangle = orientation(triangle,p2,tempVars.tempVarsSoSPredicatesImpl);
  if(orientationP1Triangle==orientationP2Triangle) return false; //both are on the same side of the triangle's plane...

  const InputVertex &a = *(triangle.getInputVertex(0));
  const InputVertex &b = *(triangle.getInputVertex(1));
  const InputVertex &c = *(triangle.getInputVertex(2));

  if(orientationP1Triangle>0) return (orientation(a,b,p2,p1,tempVars.tempVarsSoSPredicatesImpl)>0) && (orientation(b,c,p2,p1,tempVars.tempVarsSoSPredicatesImpl)>0) && (orientation(c,a,p2,p1,tempVars.tempVarsSoSPredicatesImpl)>0);
  else return (orientation(a,b,p1,p2,tempVars.tempVarsSoSPredicatesImpl)>0) && (orientation(b,c,p1,p2,tempVars.tempVarsSoSPredicatesImpl)>0) && (orientation(c,a,p1,p2,tempVars.tempVarsSoSPredicatesImpl)>0);
}



//If intersection happens at boundaries for example --> we really need SoS to be reliable...
//ex: maybe w/o SoS we would have edge-triangle, but with SoS the intersection is triangle-edge...
//think a little more to make sure we really cannot compute the coordinates..
bool MeshIntersectionGeometry::intersectTwoTrianglesSoSImpl(const InputTriangle &triMesh0,const InputTriangle &triMesh1,
            VertexFromIntersection &vertexThatCreatedPt1, VertexFromIntersection &vertexThatCreatedPt2, TempVarsComputeIntersections &tempVars)  {
  //because of SoS, we will have exactly two edge-triangle intersections
  //we need to test all 6 possibilities of edge-triangle intersections...

  int numIntersectionsFound = 0;

  for(int i=0;i<3;i++) {
    const InputVertex &v0 = *(triMesh1.getInputVertex(i));
    const InputVertex &v1 = *(triMesh1.getInputVertex((i+1)%3));
    if(intersectEdgeWithTriangleSoSImpl(triMesh0 ,v0,v1,tempVars)) {
      if(numIntersectionsFound==0) {
        vertexThatCreatedPt1.triangle = triMesh0;
        vertexThatCreatedPt1.setEdges(v0,v1);
        //vertexThatCreatedPt1.edge[0] = v0;
        //vertexThatCreatedPt1.edge[1] = v1;
      } else {
        vertexThatCreatedPt2.triangle = triMesh0;
        vertexThatCreatedPt2.setEdges(v0,v1);
        //vertexThatCreatedPt2.edge[0] = v0;
        //vertexThatCreatedPt2.edge[1] = v1;
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
        vertexThatCreatedPt1.setEdges(v0,v1);
        //vertexThatCreatedPt1.edge[0] = v0;
        //vertexThatCreatedPt1.edge[1] = v1;
      } else {
        vertexThatCreatedPt2.triangle = triMesh1;
        vertexThatCreatedPt2.setEdges(v0,v1);
        //vertexThatCreatedPt2.edge[0] = v0;
        //vertexThatCreatedPt2.edge[1] = v1;
      }
      numIntersectionsFound++;
    }
  }



  assert(numIntersectionsFound==0 || numIntersectionsFound==2); //with SoS, we can have only either 2 or 0 intersections between edges of a triangle and another triangle...
  return numIntersectionsFound>0;
}






//Given two triangles above a point, where the height above point is equal for both triangles, decide which one is lower according after SoS
//Make sure the two triangles are not the same... (this could happen, but we avoid that by checkin in PinMesh... )
/*
Getting the best triangle by perturbing the triangles by epsV is equivalent to get the best one by perturbing the query point by -epsV

Suppose there are two triangles t1, t2 with plane equations: a1 x + b1 y + c1 z = d1 and a2 x + b2 y + c2 z = d2

Then, z1 = (d1 - (a1 x + b1 y))/c1  (c1 cannot be 0 otherwise the triangle would be vertical)
And z2 =  (d2 - (a2 x + b2 y))/c2

To check if t1 < t2, we have to check if: 
z1  < z2 
(d1 - (a1 x + b1 y))/c1 <  (d2 - (a2 x + b2 y))/c2 
c2 d1 - c2(a1 x + b1 y) < c1d2 - c1 (a2 x + b2 y)

Considering the query point q is perturbed (using a negative or positive epsilon), we have its coordinates equal to (x+ (2mesh-1)eps,  y +  (2mesh-1)eps^2).
c2 d1 - c2(a1 x + a1 eps (2mesh-1) + b1 y + b1 eps^2  (2mesh-1)) < c1 d2 - c1(a2 x + a2 eps (2mesh-1) + b2 y + b2 eps^2  (2mesh-1)) 

c2 d1 - c2(a1 x + b1 y) - c2 (a1 eps (2mesh-1) + b1 eps^2  (2mesh-1))  < c1 d2 - c1(a2 x + b2 y ) - c1(a2 eps (2mesh-1) + b2 eps^2  (2mesh-1))

Since the heights are the same at (x,y), we can elliminate the non-epsilon terms (they sum to 0):
 - c2 (a1 eps (2mesh-1) + b1 eps^2  (2mesh-1))  < - c1(a2 eps (2mesh-1) + b2 eps^2  (2mesh-1))
 c2 (a1 eps (2mesh-1) + b1 eps^2  (2mesh-1))  > c1(a2 eps (2mesh-1) + b2 eps^2  (2mesh-1))
 */
const void computePlaneEquationsInputTriangle(const Point &t0,const Point &t1,const Point &t2, VertCoord &a, VertCoord &b, VertCoord &c) ;
const InputTriangle * MeshIntersectionGeometry::getBestTrianglePointInObjectSoSImpl(const InputTriangle *candidateTriangle,const InputTriangle *bestTriangle, const InputVertex &p,TempVarGetBestTrianglePointInObjectSoS &tempVars)  {
  //Let candidate triangle be triangle 1 and bestTriangle be triangle2...

  //cerr << "Getting best triangle" << endl;
  //c2 (a1 eps (2mesh-1) + b1 eps^2  (2mesh-1))  > c1(a2 eps (2mesh-1) + b2 eps^2  (2mesh-1))
  VertCoord a1,b1,c1; //plane equation of triangle 1 (triangle 1: a1 x + b1 y + c1 z = d1 -- we do not need the "D")
  VertCoord a2,b2,c2; 

  //Mesh 0 is unperturbed, mesh 1 is perturbed by (eps,eps^2, eps^3)
  //if the query point is in mesh 1 --> the computation should be similar to PinMesh's
  //otherwise, it is like translating the query point by -(eps,eps^2, eps^3) and keeping the triangles at the same place

  int meshIdQueryPoint = p.getMeshId();
  int epsSignal = (2*meshIdQueryPoint - 1); //1 if meshId is 1, -1 if 0


  computePlaneEquationsInputTriangle(getCoordinates(*(candidateTriangle->getInputVertex(0))),getCoordinates(*(candidateTriangle->getInputVertex(1))),getCoordinates(*(candidateTriangle->getInputVertex(2))),a1,b1,c1);
  computePlaneEquationsInputTriangle(getCoordinates(*(bestTriangle->getInputVertex(0))),getCoordinates(*(bestTriangle->getInputVertex(1))),getCoordinates(*(bestTriangle->getInputVertex(2))),a2,b2,c2);


  //Because of SoS, let's first check the first epsilon coefficient (eps):
  //c2 (a1 eps (2mesh-1))  > c1(a2 eps (2mesh-1) ) --> c2 (a1 eps (2mesh-1))  - c1(a2 eps (2mesh-1) ) > 0 --> (c2 (a1)  - c1(a2))*eps (2mesh-1) > 0
  int signEps1Term = sgn( a1/c1 - a2/c2 )*epsSignal;
  if(signEps1Term==1) return candidateTriangle;
  else if(signEps1Term==0) {
    int signEps2Term = sgn( b1/c1 - b2/c2 )*epsSignal; ////c2 (b1 ) - c1(b2) * (eps^2  (2mesh-1))
    assert(signEps2Term!=0); //the triangle shouldn't be vertical...  both signals can't be 0...
    if(signEps2Term==1) return candidateTriangle;
    else return bestTriangle;
  }
  else return bestTriangle;
}

//we do not need the "d" of the plane equation
const void computePlaneEquationsInputTriangle(const Point &t0,const Point &t1,const Point &t2, VertCoord &a, VertCoord &b, VertCoord &c)  {
  Point vec[2];
  for(int i=0;i<3;i++) {vec[0][i] = t1[i]-t0[i];} //t1-t0
  for(int i=0;i<3;i++) {vec[1][i] = t2[i]-t0[i];} //t2-t0


  //lets compute the values A,B,C basing on the two 3D vectors vec[0] and vec[1]
  //  | i             j           k      i             j       
  //  | vec[0][0]  vec[0][1]  vec[0][2]  vec[0][0]  vec[0][1]  
  //  | vec[1][0]  vec[1][1]  vec[1][2]  vec[1][0]  vec[1][1] 

  a = vec[0][1]*vec[1][2] - vec[0][2]*vec[1][1];
  b = vec[0][2]*vec[1][0] - vec[0][0]*vec[1][2];
  c = vec[0][0]*vec[1][1] - vec[0][1]*vec[1][0];    
}