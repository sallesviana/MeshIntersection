#include <queue>

#include "triangleClassification.h"
#include "PinMesh.h"

bool isVertexSharedBetweenMeshes(int vertexId);
Point * getPointFromVertexId(int vertexId, int meshId);

int labelConnectedComponentsEachVertex(const vector<list<int> > &adjList,vector<int> &connectedComponentEachVertex,vector<int> &sampleVertexIdFromEachConnectedComponent) {
  int numComponentsSeenSoFar = 0;
  queue< int >   vertexToLabel;
  int numVertices = adjList.size();

  vector<bool> shouldLabel(adjList.size(),true);
  for(int i=0;i<numVertices;i++) 
    if(shouldLabel[i]) {
      vertexToLabel.push(i);
      connectedComponentEachVertex[i] = numComponentsSeenSoFar;
      shouldLabel[i] = false;
      sampleVertexIdFromEachConnectedComponent.push_back(i);

      while(!vertexToLabel.empty()) {
        int v = vertexToLabel.front();
        vertexToLabel.pop(); 
        

        for(int neighbor:adjList[v]) {
          if(shouldLabel[neighbor]) {
            shouldLabel[neighbor] = false;
            connectedComponentEachVertex[neighbor] = numComponentsSeenSoFar;
            vertexToLabel.push(neighbor);
          }
        }
      }

      numComponentsSeenSoFar++;
    }
  

  
  

  return numComponentsSeenSoFar;
}


int labelConnectedComponentsEachVertex(const vector<vector<int> > &adjList,vector<int> &connectedComponentEachVertex,vector<int> &sampleVertexIdFromEachConnectedComponent) {
  int numComponentsSeenSoFar = 0;
  queue< int >   vertexToLabel;
  int numVertices = adjList.size();

  vector<bool> shouldLabel(adjList.size(),true);
  for(int i=0;i<numVertices;i++) 
    if(shouldLabel[i]) {
      vertexToLabel.push(i);
      connectedComponentEachVertex[i] = numComponentsSeenSoFar;
      shouldLabel[i] = false;
      sampleVertexIdFromEachConnectedComponent.push_back(i);

      while(!vertexToLabel.empty()) {
        int v = vertexToLabel.front();
        vertexToLabel.pop(); 
        

        for(int neighbor:adjList[v]) {
          if(shouldLabel[neighbor]) {
            shouldLabel[neighbor] = false;
            connectedComponentEachVertex[neighbor] = numComponentsSeenSoFar;
            vertexToLabel.push(neighbor);
          }
        }
      }

      numComponentsSeenSoFar++;
    }
   

  return numComponentsSeenSoFar;
}


//Checks if at least one of the vertices in the interval is not shared between the two meshes
//Returns a pointer for the first vertex (pointer to the first vertex id) that is non shared
//(or NULL if there is no such a vertex)
//[begin,end) is the interval of vertices representing the polygon
const int* getNonSharedVertextFromPolygon(const int *begin,const int *end) {
  while(begin!=end) {
    if(!isVertexSharedBetweenMeshes(*begin)) return begin;
    begin++;
  }
  return NULL;  
}


void locateTrianglesAndPolygonsInOtherMesh(const Nested3DGridWrapper *uniformGrid, 
                                            const unordered_set<const Triangle *> trianglesThatIntersect[2],
                                            vector<BoundaryPolygon> polygonsFromRetesselation[2],
                                            int thisMeshId,
                                            vector<ObjectId> &locationOfEachNonIntersectingTrianglesInOtherMesh,
                                            vector<ObjectId> &locationOfPolygonsFromRetesselationInOtherMesh,
                                            const vector<Point> vertices[3],
                                            const vector<Triangle> triangles[2]) {  
    const int DONT_KNOW_FLAG = -999999999;
    timespec t0,t1,t0Function;

    
    clock_gettime(CLOCK_REALTIME, &t0);
    t0Function  = t0;

    vector<const Point *> verticesToLocateInOtherMesh;

    clock_gettime(CLOCK_REALTIME, &t0);
    
    vector<ObjectId> locationEachVertex(vertices[thisMeshId].size(),DONT_KNOW_FLAG);
    vector<int> connectedComponentEachVertex(vertices[thisMeshId].size(),DONT_KNOW_FLAG);
    vector<int> sampleVertexIdFromEachConnectedComponent;
    int numComponents;

{
    vector<vector<int> > adjList(vertices[thisMeshId].size());
    for(const Triangle&t:triangles[thisMeshId]) {
      if(trianglesThatIntersect[thisMeshId].count(&t)==0) { //this triangle does not intersect the other mesh...
        adjList[t.p[0]].push_back(t.p[1]);
        adjList[t.p[1]].push_back(t.p[2]);
        adjList[t.p[2]].push_back(t.p[0]);
        
        adjList[t.p[1]].push_back(t.p[0]);
        adjList[t.p[2]].push_back(t.p[1]);
        adjList[t.p[0]].push_back(t.p[2]);
      } 
    }

    clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "Total time to create adj. list: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
    clock_gettime(CLOCK_REALTIME, &t0);

    

    cerr << "Labeling connected components\n";
    
    numComponents = labelConnectedComponentsEachVertex(adjList,connectedComponentEachVertex,sampleVertexIdFromEachConnectedComponent);
    assert(numComponents==sampleVertexIdFromEachConnectedComponent.size());

    clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "Total time to compute CCs: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
    clock_gettime(CLOCK_REALTIME, &t0);
}


    cerr << "Num connected components to locate: " << numComponents << "\n";

    clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "Total time to free adj. list memory: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;


    clock_gettime(CLOCK_REALTIME, &t0);

    vector<int> global_x_coord_vertex_to_locate,global_y_coord_vertex_to_locate,global_z_coord_vertex_to_locate;
    for(int i=0;i<numComponents;i++) {      
      int vertexToLocate = sampleVertexIdFromEachConnectedComponent[i];

      verticesToLocateInOtherMesh.push_back(&(vertices[thisMeshId][vertexToLocate]));     
      global_x_coord_vertex_to_locate.push_back(uniformGrid->get_global_x_coord_mesh_vertex(thisMeshId,vertexToLocate));  
      global_y_coord_vertex_to_locate.push_back(uniformGrid->get_global_y_coord_mesh_vertex(thisMeshId,vertexToLocate));  
      global_z_coord_vertex_to_locate.push_back(uniformGrid->get_global_z_coord_mesh_vertex(thisMeshId,vertexToLocate));         
    }


    int posStartVerticesOfIntersectingTrianglesInThisMesh = verticesToLocateInOtherMesh.size();
    const int numPolygonsFromRetesselation = polygonsFromRetesselation[thisMeshId].size();
    cerr << "Mesh " << thisMeshId << " Num polygons from retesselation: " << numPolygonsFromRetesselation << endl;
    
    

 
    
    clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "Total time to select vertices to locate: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
    clock_gettime(CLOCK_REALTIME, &t0);




    //Now, let's locate the triangles/polygons from retesselation in the other mesh
    int numberPolygonsFromRetesselationWithNonSharedVertices =0;
    int numberPolygonsFromRetesselationWithOnlySharedVertices =0;


    vector<int> polygonFromRetesselationWithOnlySharedVertices;
    for(int tid=0;tid<numPolygonsFromRetesselation;tid++) {
        const BoundaryPolygon& polygon = polygonsFromRetesselation[thisMeshId][tid];
        if(getNonSharedVertextFromPolygon(&(*polygon.vertexSequence.begin()),&(*polygon.vertexSequence.end()))!=NULL)  {
          numberPolygonsFromRetesselationWithNonSharedVertices++;
        }
        else {          
          polygonFromRetesselationWithOnlySharedVertices.push_back(tid);     
          numberPolygonsFromRetesselationWithOnlySharedVertices++;   
        }
    }

    //Resize all vectors to fit the data related to the intersecting triangles...
    global_x_coord_vertex_to_locate.resize(posStartVerticesOfIntersectingTrianglesInThisMesh+numberPolygonsFromRetesselationWithOnlySharedVertices);
    global_y_coord_vertex_to_locate.resize(posStartVerticesOfIntersectingTrianglesInThisMesh+numberPolygonsFromRetesselationWithOnlySharedVertices);
    global_z_coord_vertex_to_locate.resize(posStartVerticesOfIntersectingTrianglesInThisMesh+numberPolygonsFromRetesselationWithOnlySharedVertices);
    verticesToLocateInOtherMesh.resize(posStartVerticesOfIntersectingTrianglesInThisMesh+numberPolygonsFromRetesselationWithOnlySharedVertices);
    vector<Point> centerOfIntersectingTriangles(numberPolygonsFromRetesselationWithOnlySharedVertices);

    //Polygons with only shared vertices cannot be located by just locating one of its vertices (except SoS... TODO)
    //Thus, we need to locate a vertex in its interior... TODO
    cerr << "Number of polygons from retesselation with all vertices shared: " << polygonFromRetesselationWithOnlySharedVertices.size() << "\n";
    cerr << "Number of polygons from retesselation with non shared vertices: " << numberPolygonsFromRetesselationWithNonSharedVertices << "\n";
   
    /*
    #pragma omp parallel 
    {
      VertCoord tempVar;
      big_int tempVarsInt[3];

      #pragma omp for 
      for(int t =0;t<polygonFromRetesselationWithOnlySharedVertices.size();t++) {
        int tid = polygonFromRetesselationWithOnlySharedVertices[t];
        const TriangleNoBB& triangle = trianglesFromRetesselation[thisMeshId][tid];

        Point *a = getPointFromVertexId(triangle.p[0], thisMeshId);
        Point *b = getPointFromVertexId(triangle.p[1], thisMeshId);
        Point *c = getPointFromVertexId(triangle.p[2], thisMeshId);
        
        //for(int i=0;i<3;i++) centerOfIntersectingTriangles[tid][i] = 0;
        for(int i=0;i<3;i++) {
          //centerOfIntersectingTriangles[tid][i] += (*a)[i];
          //centerOfIntersectingTriangles[tid][i] += (*b)[i];
          //centerOfIntersectingTriangles[tid][i] += (*c)[i];
          //centerOfIntersectingTriangles[tid][i] /= 3;
          tempVar = (*a)[i];
          tempVar += (*b)[i];
          tempVar += (*c)[i];
          tempVar /= 3;
          centerOfIntersectingTriangles[t][i] = tempVar;
        }  
        verticesToLocateInOtherMesh[t +posStartVerticesOfIntersectingTrianglesInThisMesh] = (&centerOfIntersectingTriangles[t]);
        global_x_coord_vertex_to_locate[t +posStartVerticesOfIntersectingTrianglesInThisMesh] = (uniformGrid->x_global_cell_from_coord(centerOfIntersectingTriangles[t][0], tempVar,tempVarsInt)); 
        global_y_coord_vertex_to_locate[t +posStartVerticesOfIntersectingTrianglesInThisMesh] = (uniformGrid->y_global_cell_from_coord(centerOfIntersectingTriangles[t][1], tempVar,tempVarsInt)); 
        global_z_coord_vertex_to_locate[t +posStartVerticesOfIntersectingTrianglesInThisMesh] = (uniformGrid->z_global_cell_from_coord(centerOfIntersectingTriangles[t][2], tempVar,tempVarsInt));      
                
      }
    }
    */
    
    //TODO: locate unique vertices???
    vector<ObjectId> locationOfEachVertexInOtherMesh(verticesToLocateInOtherMesh.size());
    //vertices of mesh "meshId" will be located in mesh "1-meshId"

    clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "Total time to compute center and get grid cell of center of intersecting triangles: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;


    //Locate the vertices...    

    PinMesh pointLocationAlgorithm(uniformGrid,vertices,triangles);
    clock_gettime(CLOCK_REALTIME, &t0);
    pointLocationAlgorithm.locateVerticesInObject(verticesToLocateInOtherMesh,global_x_coord_vertex_to_locate,global_y_coord_vertex_to_locate,global_z_coord_vertex_to_locate,locationOfEachVertexInOtherMesh,1-thisMeshId);
    clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "Total time to locate: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
    //now, we know in what object of the other mesh each triangle that does not intersect other triangles is...
    clock_gettime(CLOCK_REALTIME, &t0);




    
    locationOfEachNonIntersectingTrianglesInOtherMesh = vector<ObjectId>(triangles[thisMeshId].size(),DONT_KNOW_FLAG);
    locationOfPolygonsFromRetesselationInOtherMesh = vector<ObjectId>(numPolygonsFromRetesselation,DONT_KNOW_FLAG);


    int ct =0;
    for(int tid = 0; tid < triangles[thisMeshId].size();tid++) {
      const Triangle &t = triangles[thisMeshId][tid];
      if(trianglesThatIntersect[thisMeshId].count(&t)==0) {
        int connectedComponentOfThisVertex = connectedComponentEachVertex[t.p[0]];
        locationOfEachNonIntersectingTrianglesInOtherMesh[tid] = locationOfEachVertexInOtherMesh[connectedComponentOfThisVertex];
      }
    }

    /*

    #pragma omp parallel for 
    for(int tid=0;tid<numPolygonsFromRetesselation;tid++) {
      const TriangleNoBB& t = trianglesFromRetesselation[thisMeshId][tid];
      if(!isVertexSharedBetweenMeshes(t.p[0]) || !isVertexSharedBetweenMeshes(t.p[1]) || !isVertexSharedBetweenMeshes(t.p[2]) ) {
        int vertexFromOriginalTriangle = -1;
        if(!isVertexSharedBetweenMeshes(t.p[0]))
          vertexFromOriginalTriangle = t.p[0];
        else if(!isVertexSharedBetweenMeshes(t.p[1]))
          vertexFromOriginalTriangle = t.p[1];
        else if(!isVertexSharedBetweenMeshes(t.p[2]))
          vertexFromOriginalTriangle = t.p[2];

        int connectedComponentOfThisVertex = connectedComponentEachVertex[vertexFromOriginalTriangle];
        locationOfTrianglesFromRetesselationInTheOtherMesh[tid] = locationOfEachVertexInOtherMesh[connectedComponentOfThisVertex];
      } 
    }

         

      #pragma omp parallel for 
      for(int t =0;t<triangleIdToProcessNonIntersectingTriangleWithoutSharedVertices.size();t++) {
        int tid = triangleIdToProcessNonIntersectingTriangleWithoutSharedVertices[t];
        
        assert(locationOfTrianglesFromRetesselationInTheOtherMesh[tid] ==DONT_KNOW_FLAG);
    
          //triangleIdToProcessNonIntersectingTriangleWithoutSharedVertices.push_back(numberIntersectingTrianglesWithoutNonSharedVertices);     
          //numberIntersectingTrianglesWithoutNonSharedVertices++;   
        locationOfTrianglesFromRetesselationInTheOtherMesh[tid] = locationOfEachVertexInOtherMesh[t +posStartVerticesOfIntersectingTrianglesInThisMesh];
      }
     // cerr << "Ct set: " << ctSet << " " << locationOfTrianglesFromRetesselationInTheOtherMesh.size() << "\n";
      /*
      if(!isVertexSharedBetweenMeshes(t.p[0]) || !isVertexSharedBetweenMeshes(t.p[1]) || !isVertexSharedBetweenMeshes(t.p[2]) ) {
          int vertexToLocate = -1;
          if(!isVertexSharedBetweenMeshes(t.p[0]))
            vertexToLocate = t.p[0];
          else if(!isVertexSharedBetweenMeshes(t.p[1]))
            vertexToLocate = t.p[1];
          else if(!isVertexSharedBetweenMeshes(t.p[2]))
            vertexToLocate = t.p[2];

          verticesToLocateInOtherMesh[tid+posStartVerticesOfIntersectingTrianglesInThisMesh] = (&(vertices[thisMeshId][vertexToLocate]));     
          global_x_coord_vertex_to_locate[tid+posStartVerticesOfIntersectingTrianglesInThisMesh] = (uniformGrid->get_global_x_coord_mesh_vertex(thisMeshId,vertexToLocate));  
          global_y_coord_vertex_to_locate[tid+posStartVerticesOfIntersectingTrianglesInThisMesh] = (uniformGrid->get_global_y_coord_mesh_vertex(thisMeshId,vertexToLocate));  
          global_z_coord_vertex_to_locate[tid+posStartVerticesOfIntersectingTrianglesInThisMesh] = (uniformGrid->get_global_z_coord_mesh_vertex(thisMeshId,vertexToLocate));     
      } 

      locationOfTrianglesFromRetesselationInTheOtherMesh[tid] = locationOfEachVertexInOtherMesh[tid+posStartVerticesOfIntersectingTrianglesInThisMesh];
      */
    clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "Total copy triangle labels: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;


    
    clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "Total entire location functions: " << convertTimeMsecs(diff(t0Function,t1))/1000 << endl;
}



void triangulatePolygonsFromRetesselation(vector<BoundaryPolygon> &polygonsFromRetesselation, const vector<Point> vertices[3],int meshIdWherePolygonIs) {
  int numPolygons = polygonsFromRetesselation.size();

  int onePercentNumPolygons = numPolygons/100;
  if(onePercentNumPolygons==0) onePercentNumPolygons = 1;

  #pragma omp parallel for
  for(int i=0;i<numPolygons;i++) {
    //cerr << "Triangulating " << i << " of " << numPolygons << " Percent= " << i*100/numPolygons << endl;
    if((i%onePercentNumPolygons)==0) clog << "Triangulating " << i << " of " << numPolygons << " Percent= " << i*100/numPolygons << "\n";
    polygonsFromRetesselation[i].triangulatePolygon(vertices,meshIdWherePolygonIs);
  }
  //cerr << "End of triangulations..." << endl;
}



//TODO: if all the tree vertices of the original triangle are in the same grid cell --> the new triangles will also be!

//----------------------------------------------------------------------------
// TODO: not all vertices from the retesselated triangles will be in the output (?)... we do not need to store them!
/*//Each vector represents the vertices of a layer
vector<Point> vertices[2];

//Each vector represents a set of objects in the same layer
//The objects are represented by a set of triangles (defining their boundaries)
vector<Triangle> triangles[2]; //
*/
void classifyTrianglesAndGenerateOutput(const Nested3DGridWrapper *uniformGrid, const unordered_set<const Triangle *> trianglesThatIntersect[2],vector<BoundaryPolygon> polygonsFromRetesselation[2], const vector<Point> vertices[3],const vector<Triangle> triangles[2],ostream &outputStream) {
	timespec t0,t0ThisFunction,t1;
	clock_gettime(CLOCK_REALTIME, &t0);
	t0ThisFunction = t0;

	//int ctIntersectingTrianglesTotal =0;
	vector<Triangle> outputTriangles[2]; //output triangles generated from triangles of each mesh
	vector<TriangleNoBB> outputTrianglesFromRetesselation[2];
	for(int meshId=0;meshId<2;meshId++){
    timespec t0,t1;

    //before locating the polygons from retesselation we need to triangulate them...
    triangulatePolygonsFromRetesselation(polygonsFromRetesselation[0], vertices,0);
    triangulatePolygonsFromRetesselation(polygonsFromRetesselation[1], vertices,1);


    vector<ObjectId> locationOfEachNonIntersectingTrianglesInOtherMesh,locationOfTrianglesFromRetesselationInTheOtherMesh;

    clock_gettime(CLOCK_REALTIME, &t0);
    locateTrianglesAndPolygonsInOtherMesh(uniformGrid,trianglesThatIntersect,polygonsFromRetesselation,meshId,locationOfEachNonIntersectingTrianglesInOtherMesh,locationOfTrianglesFromRetesselationInTheOtherMesh,vertices,triangles);
     clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "# Time to locate vertices in other mesh: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
    clock_gettime(CLOCK_REALTIME, &t0);

		int numTrianglesThisMesh = triangles[meshId].size();
		int ctTrianglesProcessed = 0;


		for(int i=0;i<numTrianglesThisMesh;i++) {
			const Triangle &t=triangles[meshId][i];
			if(trianglesThatIntersect[meshId].count(&t)==0) { //this triangle does not intersect the other mesh...
				//this will (probably) be an output triangle...
				ObjectId objWhereTriangleIs = locationOfEachNonIntersectingTrianglesInOtherMesh[i];//locationOfEachVertexInOtherMesh[ctTrianglesProcessed++];		
				//cerr << "obj: " << objWhereTriangleIs << endl;		
				if (objWhereTriangleIs!=OUTSIDE_OBJECT) {
					//if the triangle is not outside the other mesh, it will be in the output (we still need to update the left/right objects correctly...)
					outputTriangles[meshId].push_back(t);
				}				
			}			
		}


    //TODO
    /*
		int numTrianglesFromIntersectionThisMesh = trianglesFromRetesselation[meshId].size();
		
		for(int i=0;i<numTrianglesFromIntersectionThisMesh;i++) {
			const TriangleNoBB&t = trianglesFromRetesselation[meshId][i];
			ObjectId objWhereTriangleIs = locationOfTrianglesFromRetesselationInTheOtherMesh[i];//locationOfEachVertexInOtherMesh[ctTrianglesProcessed++];					
			
			if (objWhereTriangleIs!=OUTSIDE_OBJECT) {				

					//if the triangle is not outside the other mesh, it will be in the output (we still need to update the left/right objects correctly...)
					outputTrianglesFromRetesselation[meshId].push_back(t);
			}								
  	}*/  

    //for each polygon   


    clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "# Time to classify triangles vertices in other mesh: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
    clock_gettime(CLOCK_REALTIME, &t0);
	}
	clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to locate vertices and classify triangles: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  
  clock_gettime(CLOCK_REALTIME, &t0); 
  
	vector<bool> verticesToWriteOutputMesh0(vertices[0].size(),false);
	vector<int> newIdVerticesMesh0(vertices[0].size(),-1);
	for(const Triangle &t : outputTriangles[0]) {
		verticesToWriteOutputMesh0[t.p[0]] = verticesToWriteOutputMesh0[t.p[1]] = verticesToWriteOutputMesh0[t.p[2]] = true;
	}
	for(TriangleNoBB &t : outputTrianglesFromRetesselation[0]) {
		//cerr << t.p[0] << " " << vertices[0].size() << endl;
		assert(t.p[0] < (int)vertices[0].size());
		assert(t.p[1] < (int)vertices[0].size());
		assert(t.p[2] < (int)vertices[0].size());

		if(t.p[0]>=0) verticesToWriteOutputMesh0[t.p[0]] = true;
		if(t.p[1]>=0) verticesToWriteOutputMesh0[t.p[1]] = true;
		if(t.p[2]>=0) verticesToWriteOutputMesh0[t.p[2]] = true;
	}
	int numVerticesMesh0InOutput = 0;
	int sz = vertices[0].size();
	for(int i=0;i<sz;i++) 
		if(verticesToWriteOutputMesh0[i]) {			
			newIdVerticesMesh0[i] = numVerticesMesh0InOutput;
			numVerticesMesh0InOutput++;
		}

	vector<bool> verticesToWriteOutputMesh1(vertices[1].size(),false);
	vector<int> newIdVerticesMesh1(vertices[1].size(),-1);
	for(const Triangle &t : outputTriangles[1]) {
		verticesToWriteOutputMesh1[t.p[0]] = verticesToWriteOutputMesh1[t.p[1]] = verticesToWriteOutputMesh1[t.p[2]] = true;
	}
	for(TriangleNoBB &t : outputTrianglesFromRetesselation[1]) {
		assert(t.p[0] < (int)vertices[1].size());
		assert(t.p[1] < (int)vertices[1].size());
		assert(t.p[2] < (int)vertices[1].size());
		if(t.p[0]>=0) verticesToWriteOutputMesh1[t.p[0]] = true;
		if(t.p[1]>=0) verticesToWriteOutputMesh1[t.p[1]] = true;
		if(t.p[2]>=0) verticesToWriteOutputMesh1[t.p[2]] = true;
	}

	int numVerticesMesh1InOutput = 0;
	sz = vertices[1].size();
	for(int i=0;i<sz;i++) 
		if(verticesToWriteOutputMesh1[i]) {			
			newIdVerticesMesh1[i] = numVerticesMesh0InOutput + numVerticesMesh1InOutput;
			numVerticesMesh1InOutput++;
		}

	//compute the new ids of the shared vertices...

	for(Triangle &t : outputTriangles[0]) {
		t.p[0] = newIdVerticesMesh0[t.p[0]];
		t.p[1] = newIdVerticesMesh0[t.p[1]];
		t.p[2] = newIdVerticesMesh0[t.p[2]];
	}
	for(Triangle &t : outputTriangles[1]) {
		t.p[0] = newIdVerticesMesh1[t.p[0]];
		t.p[1] = newIdVerticesMesh1[t.p[1]];
		t.p[2] = newIdVerticesMesh1[t.p[2]];
	}

	
	//the new ids of the shared vertices will be i+numVerticesMesh0InOutput + numVerticesMesh1InOutput, for i = 0...num shared vertices-1.
	const int baseIdsSharedVertices = numVerticesMesh0InOutput + numVerticesMesh1InOutput -1;

	//cerr << baseIdsSharedVertices << endl;
	//cerr << newIdVerticesMesh0.size() << endl;
	//cerr << outputTrianglesFromRetesselation[0][351].p[2] << endl;
	for(TriangleNoBB &t : outputTrianglesFromRetesselation[0]) {
		t.p[0] = (t.p[0]<0)?baseIdsSharedVertices-t.p[0]:newIdVerticesMesh0[t.p[0]];
		t.p[1] = (t.p[1]<0)?baseIdsSharedVertices-t.p[1]:newIdVerticesMesh0[t.p[1]];
		t.p[2] = (t.p[2]<0)?baseIdsSharedVertices-t.p[2]:newIdVerticesMesh0[t.p[2]];
	}
	//cerr << outputTrianglesFromRetesselation[0][351].p[2] << endl;
	//
	for(TriangleNoBB &t : outputTrianglesFromRetesselation[1]) {
		t.p[0] = (t.p[0]<0)?baseIdsSharedVertices-t.p[0]:newIdVerticesMesh1[t.p[0]];
		t.p[1] = (t.p[1]<0)?baseIdsSharedVertices-t.p[1]:newIdVerticesMesh1[t.p[1]];
		t.p[2] = (t.p[2]<0)?baseIdsSharedVertices-t.p[2]:newIdVerticesMesh1[t.p[2]];
	}
	

	clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to update ids of new vertices: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  clock_gettime(CLOCK_REALTIME, &t0); 

  int totalNumberOutputVertices = numVerticesMesh0InOutput+numVerticesMesh1InOutput + vertices[2].size();
  int totalNumberOutputVerticesFromNonIntersectingTriangles = numVerticesMesh0InOutput+numVerticesMesh1InOutput;
  int totalNumberOutputTriangles = outputTriangles[0].size() + outputTriangles[1].size() + outputTrianglesFromRetesselation[0].size() + outputTrianglesFromRetesselation[1].size();

  
  unordered_map<pair<int,int>,int> edgesIds; //maybe use unordered_map for performance (if necessary...)
  vector<pair<int,int> > outputEdges;
 
	//vector<pair<int,int> > outputEdgesWithRepetition[2];

	for(int meshId=0;meshId<2;meshId++) {
		int numNewEdgesToAdd =0;
		#pragma omp parallel
		{			
			const int sz = 	outputTriangles[meshId].size();
			vector<pair<int,int> > myEdgesFound;	
			

			#pragma omp for
		  for(int i=0;i<sz;i++) {
		  	const Triangle &t = outputTriangles[meshId][i];
				int a = t.p[0];
				int b = t.p[1];
				int c = t.p[2];

				pair<int,int> e;
				if (a<b) {e.first = a; e.second = b;}
				else     {e.first = b; e.second = a;}
				myEdgesFound.push_back(e);

				if (b<c) {e.first = b; e.second = c;}
				else     {e.first = c; e.second = b;}
				myEdgesFound.push_back(e);

				if (c<a) {e.first = c; e.second = a;}
				else     {e.first = a; e.second = c;}
				myEdgesFound.push_back(e);
			}

			const int szOutputTrianglesFromIntersection = outputTrianglesFromRetesselation[meshId].size();
			#pragma omp for
		  for(int i=0;i<szOutputTrianglesFromIntersection;i++) {
		  	const TriangleNoBB &t = outputTrianglesFromRetesselation[meshId][i];
				int a = t.p[0];
				int b = t.p[1];
				int c = t.p[2];
			//	cerr << meshId << " " << i << " "  << szOutputTrianglesFromIntersection << " " << a << " " << b << " " << c << endl;
				assert(a>=0);
				assert(b>=0);
				assert(c>=0);
				
				pair<int,int> e;
				if (a<b) {e.first = a; e.second = b;}
				else     {e.first = b; e.second = a;}
				myEdgesFound.push_back(e);

				if (b<c) {e.first = b; e.second = c;}
				else     {e.first = c; e.second = b;}
				myEdgesFound.push_back(e);

				if (c<a) {e.first = c; e.second = a;}
				else     {e.first = a; e.second = c;}
				myEdgesFound.push_back(e);
			}


			sort(myEdgesFound.begin(),myEdgesFound.end());
			auto newEndItr = unique(myEdgesFound.begin(),myEdgesFound.end());
			myEdgesFound.resize(newEndItr- myEdgesFound.begin());

			

			#pragma omp critical
			{
				outputEdges.insert(outputEdges.end(),myEdgesFound.begin(),myEdgesFound.end());			
			}
		}
	}
	sort(outputEdges.begin(),outputEdges.end());
	auto newEndItr = unique(outputEdges.begin(),outputEdges.end());
	outputEdges.resize(newEndItr- outputEdges.begin());

	clock_gettime(CLOCK_REALTIME, &t1);
	cerr << "T so far: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;

	
	int totalNumberOutputEdges = outputEdges.size();
	vector<map<int,int> > mapEdgesIds2(totalNumberOutputVertices); //maps each edge (pair of vertices) to ids id
	for (int i=0;i<totalNumberOutputEdges;i++) {
		const pair<int,int> &e = outputEdges[i];
		//edgesIds[e] = i;

		assert(e.first>=0 );
		assert(e.second>=0);
		assert(e.first<totalNumberOutputVertices);
		//cerr << e.first << " " << e.second << " " << totalNumberOutputEdges << endl;
		assert(e.second<totalNumberOutputVertices);		
		mapEdgesIds2[e.first][e.second] = i;
	}

	clock_gettime(CLOCK_REALTIME, &t1);
	//timeClassifyTriangles = convertTimeMsecs(diff(t0ThisFunction,t1))/1000;

  cerr << "Time to create edges: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
 // cerr << "Total time (excluding I/O) so far: " << convertTimeMsecs(diff(t0AfterDatasetRead,t1))/1000 << endl;
  clock_gettime(CLOCK_REALTIME, &t0); 

  

  
 

  //now, let's write everything in the output!
  outputStream << totalNumberOutputVertices << " " << totalNumberOutputEdges << " " << totalNumberOutputTriangles << '\n';

  //print the coordinates of the vertices...
  sz = verticesToWriteOutputMesh0.size();
  for(int i=0;i<sz;i++) 
		if(verticesToWriteOutputMesh0[i]) 
			outputStream << vertices[0][i][0].get_d() << " " << vertices[0][i][1].get_d() << " " << vertices[0][i][2].get_d() << "\n";
	sz = verticesToWriteOutputMesh1.size();
  for(int i=0;i<sz;i++) 
		if(verticesToWriteOutputMesh1[i]) 
			outputStream << vertices[1][i][0].get_d() << " " << vertices[1][i][1].get_d() << " " << vertices[1][i][2].get_d() << "\n";
	sz = vertices[2].size();
  for(int i=0;i<sz;i++) 
			outputStream << vertices[2][i][0].get_d() << " " << vertices[2][i][1].get_d() << " " << vertices[2][i][2].get_d() << "\n";		

	//print edges...
	for(const pair<int,int> &p:outputEdges) {
		outputStream << p.first+1 << " " << p.second+1 << "\n"; //in a GTS file we start counting from 1...
	}

	//print triangles...
	for(int meshId=0;meshId<2;meshId++) 
		for(Triangle &t : outputTriangles[meshId]) {
			int a = t.p[0];
			int b = t.p[1];
			int c = t.p[2];


			if(t.above != OUTSIDE_OBJECT) {
				//according to the right hand rule, the ordering of the vertices should be (c,b,a)
				swap(a,c);
			} 
			pair<int,int> e;
			if (a<b) {e.first = a; e.second = b;}
			else     {e.first = b; e.second = a;}
			//outputStream << edgesIds[e]+1 << " "; //we start counting from 1 in GTS files...
			outputStream << mapEdgesIds2[e.first][e.second]+1 << " ";

			if (b<c) {e.first = b; e.second = c;}
			else     {e.first = c; e.second = b;}
			//outputStream << edgesIds[e]+1 << " ";
			outputStream << mapEdgesIds2[e.first][e.second]+1 << " ";

			if (c<a) {e.first = c; e.second = a;}
			else     {e.first = a; e.second = c;}
			//outputStream << edgesIds[e]+1 << "\n";
			outputStream << mapEdgesIds2[e.first][e.second]+1 << "\n";
		}

		for(int meshId=0;meshId<2;meshId++) 
			for(TriangleNoBB &t : outputTrianglesFromRetesselation[meshId]) {
				int a = t.p[0];
				int b = t.p[1];
				int c = t.p[2];

				//cerr << a << " " << b << " " << c << endl;
				if(a<0) { //if the vertex refers to a shared vertex (created from the intersection)
					a = -a -1 + totalNumberOutputVerticesFromNonIntersectingTriangles;
				}
				if(b<0) {
					b = -b -1 + totalNumberOutputVerticesFromNonIntersectingTriangles;
				}
				if(c<0) {
					c = -c -1 + totalNumberOutputVerticesFromNonIntersectingTriangles;
				}

				if(t.above != OUTSIDE_OBJECT) {
					//according to the right hand rule, the ordering of the vertices should be (c,b,a)
					swap(a,c);
				} 
				pair<int,int> e;
				if (a<b) {e.first = a; e.second = b;}
				else     {e.first = b; e.second = a;}
				//outputStream << edgesIds[e]+1 << " "; //we start counting from 1 in GTS files...
				outputStream << mapEdgesIds2[e.first][e.second]+1 << " ";

				if (b<c) {e.first = b; e.second = c;}
				else     {e.first = c; e.second = b;}
				//outputStream << edgesIds[e]+1 << " ";
				outputStream << mapEdgesIds2[e.first][e.second]+1 << " ";

				if (c<a) {e.first = c; e.second = a;}
				else     {e.first = a; e.second = c;}
				//outputStream << edgesIds[e]+1 << "\n";
				outputStream << mapEdgesIds2[e.first][e.second]+1 << "\n";
			}
	
		cerr << "Output vertices         : " << totalNumberOutputVertices << endl;
		cerr << "Output edges            : " << totalNumberOutputEdges << endl;
		cerr << "Output triangles non int: " << totalNumberOutputTriangles << endl;
		//cerr << "Intersecting triangles  : " << ctIntersectingTrianglesTotal << endl;

}