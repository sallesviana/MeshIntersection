#include "3dGeometry.h"





#include "nested3DGrid.h"


#include "geometry/3dGeometryGeometricalPredicatesSoSImpl.cpp"
#include "geometry/3dGeometryGeometricalPredicatesMainImpl.cpp"
#include "geometry/3dGeometryGeometricalPredicatesMainImplOrig.cpp"
#include "geometry/3dGeometryGeometricalPredicates.cpp"







#include "boundaryPolygon.cpp"


const MeshIntersectionGeometry::PlaneEquation &MeshIntersectionGeometry::getPlaneEquationInputTriangle(int meshId, int triId,TempVarsComputePlaneEquation &tempVars) {
	initPlaneEquationInputTriangle(meshId,triId,tempVars);
	return planeEquationsInputTriangles[meshId][triId];
}

void MeshIntersectionGeometry::initPlaneEquationInputTriangle(int meshId, int triId,TempVarsComputePlaneEquation &tempVars) {
  if(!isPlaneEquationInputTrianglesInitialized[meshId][triId]) {
    const InputTriangle &t = inputTriangles[meshId][triId];
    computePlaneEquation(planeEquationsInputTriangles[meshId][triId], getCoordinates(*t.getInputVertex(0)), getCoordinates(*t.getInputVertex(1)),getCoordinates(*t.getInputVertex(2)), tempVars);
    isPlaneEquationInputTrianglesInitialized[meshId][triId] = true;
  }
}



void MeshIntersectionGeometry::storeIntersectionVerticesCoordinatesAndUpdateVerticesIds(vector< pair<VertexFromIntersection, VertexFromIntersection> >  &edgesFromIntersection,const vector< pair<Point, Point> > &coordsVerticesOfEdges) {
  timespec t0,t1;
  clock_gettime(CLOCK_REALTIME, &t0);
  
  int numEdges = edgesFromIntersection.size();
  clog << "Number of edges: " << numEdges << "\n";

  //get pointers to the vertices in edges...
  vector< const Point *> vertPtrs(numEdges*2);

  #pragma omp parallel for 
  for(int i=0;i<numEdges;i++) {
    vertPtrs[2*i] = &(coordsVerticesOfEdges[i].first);
    vertPtrs[2*i+1] = &(coordsVerticesOfEdges[i].second);
  }

  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Vertices pointersinitialized..." <<  convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 
  clock_gettime(CLOCK_REALTIME, &t0);

  //sort the pointers basing on the vertices (indirect sorting)
  auto VertPtrCompare = [&](const Point *a, const Point *b) {
    if( (*a)[0] != (*b)[0] ) return (*a)[0] < (*b)[0];
    if( (*a)[1] != (*b)[1] ) return (*a)[1] < (*b)[1];
    return (*a)[2] < (*b)[2];
  };
  __gnu_parallel::sort(vertPtrs.begin(),vertPtrs.end(),VertPtrCompare);

  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Sort vertices: " << vertPtrs.size() << " time: "  << convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 
  clock_gettime(CLOCK_REALTIME, &t0);

  //get unique vertices (indirect)
  auto VertPtrEqCompare = [&](const Point *a, const Point *b) {
    if( (*a)[0] != (*b)[0] ) return false;
    if( (*a)[1] != (*b)[1] ) return false;
    if( (*a)[2] != (*b)[2] ) return false;
    return true;
  };
  auto it = unique (vertPtrs.begin(), vertPtrs.end(), VertPtrEqCompare);                                              
  vertPtrs.resize( distance(vertPtrs.begin(),it) );

  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Unique vertices: " << vertPtrs.size() << " time: "  << convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 
  clock_gettime(CLOCK_REALTIME, &t0);
  //copy unique vertices to vertices[2]

  const int numUniqueVerts = vertPtrs.size();
  verticesCoordinates[2].resize(numUniqueVerts);
  #pragma omp parallel for
  for(int i=0;i<numUniqueVerts;i++)
    verticesCoordinates[2][i] = *vertPtrs[i];

  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Vertices copied" << endl;
  clock_gettime(CLOCK_REALTIME, &t0);

  //use binary search to find ids? or hash? or map?
  auto getVertexId = [&](const array<VertCoord,3> *elem) {
    auto it = lower_bound(vertPtrs.begin(),vertPtrs.end(),elem,VertPtrCompare);
    #ifdef SANITY_CHECKS
      assert(it!=vertPtrs.end());
      assert(VertPtrEqCompare(*it,elem)); //it should point to the same vertex as elem points...
    #endif
    return it-vertPtrs.begin(); //the element should be in the vector...
  };

  #pragma omp parallel for
  for(int i=0;i<numEdges;i++) {
    //cerr << "i " << i << " " << numEdges << endl;
    int idA = getVertexId(&(coordsVerticesOfEdges[i].first));    
    int idB = getVertexId(&(coordsVerticesOfEdges[i].second));
    
    edgesFromIntersection[i].first.id = idA;  
    edgesFromIntersection[i].second.id = idB;

    edgesFromIntersection[i].first.meshId = 2;  
    edgesFromIntersection[i].second.meshId = 2;
  }

  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "End!" <<  convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 
}

void MeshIntersectionGeometry::computeIntersections(const vector<pair<InputTriangle *,InputTriangle *> > &inputTrianglesToConsider, 
                          vector< pair<InputTriangle *,InputTriangle *> >  &intersectingTrianglesThatGeneratedEdges, 
                          vector< pair<VertexFromIntersection, VertexFromIntersection> >  &edgesFromIntersection, 
                          unsigned long long &numIntersectionTests){ 

  const long long numPairsTrianglesToProcess = inputTrianglesToConsider.size();
  numIntersectionTests = numPairsTrianglesToProcess;

  
  cerr << "Initializing plane equations...\n";
  
  for(int meshId=0;meshId<2;meshId++) {
    //cerr << "Determining triangles to compute equations\n";
    /*vector<int> trianglesToProcess;
    if(meshId==0)
      for(int i=0;i<numPairsTrianglesToProcess;i++) 
        trianglesToProcess.push_back(inputTrianglesToConsider[i].first-&inputTriangles[0][0]);
    else
      for(int i=0;i<numPairsTrianglesToProcess;i++) 
        trianglesToProcess.push_back(inputTrianglesToConsider[i].second-&inputTriangles[1][0]);

    sort(trianglesToProcess.begin(),trianglesToProcess.end());

    vector<int>::iterator it = std::unique (trianglesToProcess.begin(), trianglesToProcess.end());
    trianglesToProcess.resize( std::distance(trianglesToProcess.begin(),it) );*/

    int numTri = inputTriangles[meshId].size();
    cerr << "Triangles to compute equations determined\n";

    #pragma omp parallel
    {
      TempVarsComputePlaneEquation tempVars;

      #pragma omp for
      for(int i=0;i<numTri;i++)
        initPlaneEquationInputTriangle(meshId, i,tempVars);
    }
  }
  cerr << "Computing the intersections...\n"; 


  vector<pair<Point,Point> > coordsVerticesOfEdges;

  #pragma omp parallel
  {
    TempVarsComputeIntersections tempVars;

    //vector<Point> coordsVerticesOfEdgesTemp[2];
    vector<pair<Point,Point> > coordsVerticesOfEdgesTemp;
    vector< pair<VertexFromIntersection, VertexFromIntersection> > edgesFromIntersectionTemp;
    vector< pair<InputTriangle *,InputTriangle *> > intersectingTrianglesThatGeneratedEdgesTemp;

    VertexFromIntersection tempVertexFromIntersection[2]; //TODO: avoid copying...
    Point tempCoordsVerticesFromIntersection[2];

    #pragma omp for
    for(int i=0;i<numPairsTrianglesToProcess;i++) {
      InputTriangle *t1 = inputTrianglesToConsider[i].first;
      InputTriangle *t2 = inputTrianglesToConsider[i].second;

      
      int coplanar;
      int intersect = intersectTwoTriangles(*t1,*t2,
             tempCoordsVerticesFromIntersection[0], tempVertexFromIntersection[0], tempCoordsVerticesFromIntersection[1],
             tempVertexFromIntersection[1],tempVars);

      if(intersect) {
	      //coordsVerticesOfEdgesTemp[0].push_back(tempCoordsVerticesFromIntersection[0]);
	      //coordsVerticesOfEdgesTemp[1].push_back(tempCoordsVerticesFromIntersection[1]);
        coordsVerticesOfEdgesTemp.push_back(make_pair(tempCoordsVerticesFromIntersection[0],tempCoordsVerticesFromIntersection[1]));
        edgesFromIntersectionTemp.push_back(make_pair(tempVertexFromIntersection[0],tempVertexFromIntersection[1]));
	      //verticesOfEdgesTemp[0].push_back(tempVertexFromIntersection[0]);
	      //verticesOfEdgesTemp[1].push_back(tempVertexFromIntersection[1]);
        intersectingTrianglesThatGeneratedEdgesTemp.push_back(inputTrianglesToConsider[i]);
	    }

    }

    #pragma omp critical
    {
        coordsVerticesOfEdges.insert(coordsVerticesOfEdges.end(),coordsVerticesOfEdgesTemp.begin(),coordsVerticesOfEdgesTemp.end());
      	//coordsVerticesOfEdges[i].insert(coordsVerticesOfEdges[i].end(),coordsVerticesOfEdgesTemp[i].begin(),coordsVerticesOfEdgesTemp[i].end());
      	edgesFromIntersection.insert(edgesFromIntersection.end(),edgesFromIntersectionTemp.begin(),edgesFromIntersectionTemp.end());
        //verticesOfEdges[i].insert(verticesOfEdges[i].end(),verticesOfEdgesTemp[i].begin(),verticesOfEdgesTemp[i].end());
        intersectingTrianglesThatGeneratedEdges.insert(intersectingTrianglesThatGeneratedEdges.end(),intersectingTrianglesThatGeneratedEdgesTemp.begin(),intersectingTrianglesThatGeneratedEdgesTemp.end());
    }
  }

  //After the intersections are computed, we need to "register" the coordinates of the new vertices and store their ids...
  storeIntersectionVerticesCoordinatesAndUpdateVerticesIds(edgesFromIntersection,coordsVerticesOfEdges);


  //#define DEBUGGING_MODE

  #ifdef DEBUGGING_MODE  
  for(int i=0;i<2;i++) {
  	cerr << "Testing i = " << i << endl;

  	int numV = edgesFromIntersection.size();
  	cerr << "Number of vertices: " << numV << endl;

  	int percent = 0;
  	for(int j=0;j<numV;j++) {
  		
  		
  		VertexFromIntersection v = edgesFromIntersection[j].first;
      if(i==1) v = edgesFromIntersection[j].second;

  		Point computed = computePointFromIntersectionVertex(v);

      const Point &stored = getCoordinates(v);
  		assert(computed==stored);

  		if(((100*j)/numV)!=percent) {
  			percent = (100*j)/numV;
  			cerr << "Percent: " << percent << endl;
  			cerr << computed[1] << "\n" << stored[1] << endl;
  		}
  	}
  }
  #endif




  cerr << "Number of edges from intersection: " << edgesFromIntersection.size() << endl;



}





#include <iomanip>


// IO and utils... //
/********MeshIntersectionGeometry************/

void MeshIntersectionGeometry::printBoundingBoxes() {
	for(int meshId=0;meshId<2;meshId++)
  	cerr << "Bounding box mesh " << meshId << ": " << setprecision(18) << std::fixed << meshBoundingBoxes[meshId][0][0].get_d() << " " << meshBoundingBoxes[meshId][0][1].get_d() <<  " " << meshBoundingBoxes[meshId][0][2].get_d() << " -- " << meshBoundingBoxes[meshId][1][0].get_d() << " " << meshBoundingBoxes[meshId][1][1].get_d() << " " << meshBoundingBoxes[meshId][1][2].get_d() <<endl;
  cerr << "Bounding box two meshes togetter " << ": " << setprecision(18) << std::fixed << boundingBoxTwoMeshesTogetter[0][0].get_d() << " " << boundingBoxTwoMeshesTogetter[0][1].get_d() <<  " " << boundingBoxTwoMeshesTogetter[0][2].get_d() << " -- " << boundingBoxTwoMeshesTogetter[1][0].get_d() << " " << boundingBoxTwoMeshesTogetter[1][1].get_d() << " " << boundingBoxTwoMeshesTogetter[1][2].get_d() <<endl;
}


MeshIntersectionGeometry::MeshIntersectionGeometry(const string &pathMesh0, const string &pathMesh1) {
	loadInputMesh(0,pathMesh0);
	loadInputMesh(1,pathMesh1);
 
	//initializes the bounding boxes...
	boundingBoxTwoMeshesTogetter[0] = meshBoundingBoxes[0][0];
  boundingBoxTwoMeshesTogetter[1] = meshBoundingBoxes[0][1];

  for(int i=0;i<3;i++) {
  	if(meshBoundingBoxes[1][0][i] < boundingBoxTwoMeshesTogetter[0][i])
  		boundingBoxTwoMeshesTogetter[0][i] = meshBoundingBoxes[1][0][i];
  	if(meshBoundingBoxes[1][1][i] > boundingBoxTwoMeshesTogetter[1][i])
  		boundingBoxTwoMeshesTogetter[1][i] = meshBoundingBoxes[1][1][i];
  }
}

void MeshIntersectionGeometry::loadInputMesh(int meshId,const string &path) {
	assert(path.size()>=4);
  string inputFileExtension = path.substr(path.size()-3,3);
  assert(inputFileExtension=="gts" || inputFileExtension=="ium");
  bool isGtsFile = inputFileExtension=="gts";
  if(isGtsFile)
    readGTSFile(path,verticesCoordinates[meshId],inputTriangles[meshId],meshBoundingBoxes[meshId],meshId);
  else 
    readLiumFile(path,verticesCoordinates[meshId],inputTriangles[meshId],meshBoundingBoxes[meshId],meshId);

  cerr << "Initializing bounding-boxes of triangles\n"; 
  for(int meshId=0;meshId<2;meshId++) {
  	const int numInputTrianglesThisMesh = inputTriangles[meshId].size();

  	inputTrianglesBoundingBox[meshId].resize(numInputTrianglesThisMesh);
  	#pragma omp parallel
  	for(int tid=0;tid<numInputTrianglesThisMesh;tid++) {
  		initTriangleBoundingBox(meshId,tid);
  	}

  	planeEquationsInputTriangles[meshId].resize(numInputTrianglesThisMesh);
  	isPlaneEquationInputTrianglesInitialized[meshId] = vector<bool>(numInputTrianglesThisMesh,false);
  }
  cerr << "Bunding-boxes of triangles initialized\n"; 
}

void MeshIntersectionGeometry::initTriangleBoundingBox(int meshId, int triangleId) {
	array<array<int,3>,2> &boundingBox = inputTrianglesBoundingBox[meshId][triangleId];
	int p0 = (inputTriangles[meshId][triangleId].getInputVertex(0))->getId();
	int p1 = (inputTriangles[meshId][triangleId].getInputVertex(1))->getId();
	int p2 = (inputTriangles[meshId][triangleId].getInputVertex(2))->getId();
	const vector<Point> &vertices = verticesCoordinates[meshId];
	for(int i=0;i<3;i++) {
			boundingBox[0][i] = p0;
			if (vertices[p1][i] < vertices[p0][i]) boundingBox[0][i] = 	p1;
			if (vertices[p2][i] < vertices[boundingBox[0][i] ][i]) boundingBox[0][i] = 	p2;
	}
	for(int i=0;i<3;i++) {
		boundingBox[1][i] = p0;
		if (vertices[p1][i] > vertices[p0][i]) boundingBox[1][i] = 	p1;
		if (vertices[p2][i] > vertices[ boundingBox[1][i] ][i]) boundingBox[1][i] = 	p2;
	}
}

//Reads a GTS file, fills the boundingBox with the boundingBox of the triangles read
void readGTSFile(string fileName, vector<Point> &vertices,vector<InputTriangle> &triangles, Point boundingBox[2],const int meshId) {
	FILE *inputFile = fopen(fileName.c_str(),"r");
	if (inputFile==NULL ) {
    cerr << "ERROR: failed to open file " << fileName << endl;
    exit(1);
  }
  cerr << "Reading file " << fileName << endl;
  int numVertices, numEdges, numTriangles;
	assert(fscanf(inputFile,"%d %d %d",&numVertices,&numEdges,&numTriangles)==3);
	vertices.resize(numVertices);
	triangles.resize(numTriangles);
	vector< array<VertexId,2> > edges(numEdges);

	//rational xr,yr,zr;
	for(int i=0;i<numVertices;i++) {
		double x,y,z;
		assert(fscanf(inputFile,"%lf %lf %lf",&x,&y,&z)==3);

		vertices[i][0] = x;
		vertices[i][1] = y;
		vertices[i][2] = z;

		if (i==0) { //initialize the bounding box in the first iteration...
			boundingBox[0] = vertices[0];
			boundingBox[1] = vertices[0];
		}

		accum_min(boundingBox[0][0],vertices[i][0]);
		accum_min(boundingBox[0][1],vertices[i][1]);
		accum_min(boundingBox[0][2],vertices[i][2]);

		accum_max(boundingBox[1][0],vertices[i][0]);
		accum_max(boundingBox[1][1],vertices[i][1]);
		accum_max(boundingBox[1][2],vertices[i][2]);
	}
	for(int i=0;i<numEdges;i++) {
		VertexId a,b;
		assert(fscanf(inputFile,"%d %d",&a,&b)==2);
		edges[i][0] = a-1; //GTS starts counting from 1... we will start the ids from 0
		edges[i][1] = b-1;
	}

	for(int i=0;i<numTriangles;i++) {
		int a,b,c;
		assert(fscanf(inputFile,"%d %d %d",&a,&b,&c)==3);
		a--;b--;c--; //we count from 0, not from 1...

		if (edges[a][0] == edges[b][0] || edges[a][0] == edges[b][1]) swap(edges[a][0], edges[a][1]); // we will ensure that the edges are in the format (a,b)-(b,c)-(d,e) or (a,b)-(c,b)-(d,e)
		if (edges[a][1] == edges[b][1]) swap(edges[b][0],edges[b][1]); //ensure that edges are in the format (a,b)-(b,c)-(d,e)
		if (edges[b][1] == edges[c][1]) swap(edges[c][0],edges[c][1]); //ensure that edges are in the format (a,b)-(b,c)-(c,a)

		assert(edges[a][1]==edges[b][0]);
		assert(edges[b][1]==edges[c][0]);
		assert(edges[c][1]==edges[a][0]);

		triangles[i] = InputTriangle( InputVertex(meshId,edges[a][0]), InputVertex(meshId,edges[a][1]),InputVertex(meshId,edges[b][1]),OUTSIDE_OBJECT,1);
	}
}


//Reads a Lium file, fills the boundingBox with the boundingBox of the triangles read
void readLiumFile(string fileName, vector<Point> &vertices,vector<InputTriangle> &triangles, Point boundingBox[2],const int meshId) {
	FILE *inputFile = fopen(fileName.c_str(),"r");
	if (inputFile==NULL ) {
    cerr << "ERROR: failed to open file " << fileName << endl;
    exit(1);
  }
  cerr << "Reading file " << fileName << endl;
  int numVertices, numTriangles,numObjects;
	assert(fscanf(inputFile,"%d %d %d",&numVertices,&numTriangles,&numObjects)==3);
	cerr << "Vertices, triangles, objects: " << numVertices << " " << numTriangles << " " << numObjects << endl;
 	vertices.resize(numVertices);
	triangles.resize(numTriangles);

	//rational xr,yr,zr;
	for(int i=0;i<numVertices;i++) {
		double x,y,z;
		assert(fscanf(inputFile,"%lf %lf %lf",&x,&y,&z)==3);

		vertices[i][0] = x;
		vertices[i][1] = y;
		vertices[i][2] = z;

		if (i==0) { //initialize the bounding box in the first iteration...
			boundingBox[0] = vertices[0];
			boundingBox[1] = vertices[0];
		}

		accum_min(boundingBox[0][0],vertices[i][0]);
		accum_min(boundingBox[0][1],vertices[i][1]);
		accum_min(boundingBox[0][2],vertices[i][2]);

		accum_max(boundingBox[1][0],vertices[i][0]);
		accum_max(boundingBox[1][1],vertices[i][1]);
		accum_max(boundingBox[1][2],vertices[i][2]);
	}
	
	for(int i=0;i<numTriangles;i++) {
		int a,b,c;
		int objContrNormal,objNormal;
		assert(fscanf(inputFile,"%d %d %d %d %d",&a,&b,&c,&objNormal,&objContrNormal)==5);
		//a--;b--;c--; //we count from 0, not from 1...
		
		objNormal++;
		objContrNormal++;
		//cerr << objNormal << " " << objContrNormal << endl;
		triangles[i] = InputTriangle(InputVertex(meshId,a),InputVertex(meshId,b),InputVertex(meshId,c),objNormal,objContrNormal);
	}
}

void MeshIntersectionGeometry::storeAllVertices(ostream &out) {
  for(int meshId=0;meshId<3;meshId++)
    for(const Point &p:verticesCoordinates[meshId]) 
      out << p[0].get_d() << " " << p[1].get_d() << " " << p[2].get_d() << "\n";    
}

//For debugging purposes
Point MeshIntersectionGeometry::computePointFromIntersectionVertex(VertexFromIntersection &vertexFromIntersection) {
	Point   u, v, n;              // triangle vectors
  Point    dir, w0, w;           // ray vectors
  VertCoord   r, a, b;              // params to calc ray-plane intersect

  const Point &V0 = getCoordinates(*(vertexFromIntersection.triangle.getInputVertex(0)));
  const Point &V1 = getCoordinates(*(vertexFromIntersection.triangle.getInputVertex(1)));
  const Point &V2 = getCoordinates(*(vertexFromIntersection.triangle.getInputVertex(2)));

  const Point &P0 = getCoordinates((vertexFromIntersection.edge[0]));
  const Point &P1 = getCoordinates((vertexFromIntersection.edge[1]));


  VertCoord tmp;


  // get triangle edge vectors and plane normal
  SUB(u,V1,V0); //u = T.V1 - T.V0;

  SUB(v,V2,V0);//v = T.V2 - T.V0;
  CROSS(n,u,v,tmp);;//n = u * v;              // cross product

 
  SUB(dir,P1,P0);//dir = R.P1 - R.P0;              // ray direction vector
  SUB(w0,P0,V0);//w0 = R.P0 - T.V0;
  DOT(a,n,w0,tmp);;//a = -dot(n,w0);
  a = -a;


  DOT(b,n,dir,tmp);//b = dot(n,dir);


  // get intersect point of ray with triangle plane
  r = a / b;

  Point ans;
  ans[0] = P0[0] + r*dir[0];
  ans[1] = P0[1] + r*dir[1];
  ans[2] = P0[2] + r*dir[2];
  //*I = R.P0 + r * dir;

  return ans;

}