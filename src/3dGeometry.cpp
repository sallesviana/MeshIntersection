#include "3dGeometry.h"







#include "3dGeometryGeometricalPart.cpp"










#include "tritri_isectline.c"
void MeshIntersectionGeometry::computeIntersections(const vector<pair<InputTriangle *,InputTriangle *> > &inputTrianglesToConsider, 
                          vector< pair<InputTriangle *,InputTriangle *> >  &intersectingTrianglesThatGeneratedEdges, 
                          vector< pair<VertexFromIntersection, VertexFromIntersection> >  &edgesFromIntersection, 
                          unsigned long long &numIntersectionTests){ 

  const long long numPairsTrianglesToProcess = inputTrianglesToConsider.size();
  numIntersectionTests = numPairsTrianglesToProcess;


  vector<Point> coordsVerticesOfEdges[2];
  vector<VertexFromIntersection> verticesOfEdges[2];
  #pragma omp parallel
  {
    TempVarsComputeIntersections tempVars;

    vector<Point> coordsVerticesOfEdgesTemp[2];
    vector<VertexFromIntersection> verticesOfEdgesTemp[2];

    VertexFromIntersection tempVertexFromIntersection[2]; //TODO: avoid copying...
    Point tempCoordsVerticesFromIntersection[2];

    #pragma omp for
    for(int i=0;i<numPairsTrianglesToProcess;i++) {
      InputTriangle *t1 = inputTrianglesToConsider[i].first;
      InputTriangle *t2 = inputTrianglesToConsider[i].second;

      
      int coplanar;
      int intersect = intersectTwoTriangles(*t1,*t2,
             tempCoordsVerticesFromIntersection[0], tempVertexFromIntersection[0], tempCoordsVerticesFromIntersection[1],
             tempVertexFromIntersection[1],tempVars.tempRationals);

      if(intersect) {
	      coordsVerticesOfEdgesTemp[0].push_back(tempCoordsVerticesFromIntersection[0]);
	      coordsVerticesOfEdgesTemp[1].push_back(tempCoordsVerticesFromIntersection[1]);

	      verticesOfEdgesTemp[0].push_back(tempVertexFromIntersection[0]);
	      verticesOfEdgesTemp[1].push_back(tempVertexFromIntersection[1]);
	    }

    }

    #pragma omp critical
    {
    	for(int i=0;i<2;i++) {
      	coordsVerticesOfEdges[i].insert(coordsVerticesOfEdges[i].end(),coordsVerticesOfEdgesTemp[i].begin(),coordsVerticesOfEdgesTemp[i].end());
      	verticesOfEdges[i].insert(verticesOfEdges[i].end(),verticesOfEdgesTemp[i].begin(),verticesOfEdgesTemp[i].end());
    	}
    }

  }

  cerr << "Number of edges from intersection: " << verticesOfEdges[0].size() << endl;

  //TODO: create edges, vertices, etc.
  //TODO: improve tri-tri intersection: reduce copies, avoid re-computing the normals every time

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
