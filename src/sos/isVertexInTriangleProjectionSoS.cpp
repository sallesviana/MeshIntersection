
int MeshIntersectionGeometry::orientation(const InputVertex &v1,
																					 const InputVertex &v2, 
																					 const InputVertex &queryPoint, 
																					 int whatPlaneProjectTrianglesTo) {	

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

int MeshIntersectionGeometry::orientation(const InputVertex &v1,
																					 const InputVertex &v2, 
																					 const VertexFromIntersection &queryPoint, 
																					 int whatPlaneProjectTrianglesTo) {	

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

int MeshIntersectionGeometry::orientation(const InputVertex &v1,
																					 const VertexFromIntersection &v2, 
																					 const VertexFromIntersection &queryPoint, 
																					 int whatPlaneProjectTrianglesTo) {	

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

int MeshIntersectionGeometry::orientation(const VertexFromIntersection &v1,
																					 const VertexFromIntersection &v2, 
																					 const VertexFromIntersection &queryPoint, 
																					 int whatPlaneProjectTrianglesTo) {	

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
//1 if make a left turn, -1 if make a right turn, 0 --> degeneracy (shouldn't happen...)
int MeshIntersectionGeometry::orientation(const Vertex &v1,
																					 const Vertex &v2, 
																					 const Vertex &p, 
																					 int whatPlaneProjectTrianglesTo) {	
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
	} else { //numImputVertices is 2...
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

bool MeshIntersectionGeometry::isVertexInTriangleProjectionSoS(const Vertex &v1,const Vertex &v2, const Vertex &v3, const Vertex &queryPoint,int whatPlaneProjectTrianglesTo) {
	#pragma omp critical
	for(int i=0;i<8;i++) cerr << i << " --> " << cts[i] << endl;

	int o1 = orientation(v1,v2,queryPoint,whatPlaneProjectTrianglesTo);
	int o2 = orientation(v2,v3,queryPoint,whatPlaneProjectTrianglesTo);
	if(o1!=o2) return false;
	int o3 = orientation(v3,v1,queryPoint,whatPlaneProjectTrianglesTo);
	return o2==o3;
}