#include <queue>

#include "triangleClassification.h"


const int DONT_KNOW_FLAG = -999999999;


#include "PinMesh.h"




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
const Vertex* getNonSharedVertextFromPolygon(const Vertex * const*begin,const Vertex * const*end) {
  while(begin!=end) {
    if(!(*begin)->isFromIntersection()) return *begin;
    begin++;
  }
  return NULL;  
}



void locatePolygonsOtherMeshUsingDFS(const MeshIntersectionGeometry &geometry,
                                    pair<const InputTriangle *,vector<BoundaryPolygon>> &tri) {
  stack< BoundaryPolygon * > bpToProcess;

  for(BoundaryPolygon&bp:tri.second) {
    if(bp.getPolyhedronWherePolygonIs()!=DONT_KNOW_ID) { //do we know the location of this polygon?
      bpToProcess.push(&bp);
    }    
  }

  assert(!bpToProcess.empty()); //at least one polygon should be known (at least the polygons containing the input vertices of the triangle...)

  while(!bpToProcess.empty()) {
    BoundaryPolygon *bp = bpToProcess.top();
    bpToProcess.pop();

    //cerr << "Found one with id... " << bp->getPolyhedronWherePolygonIs() << endl;
    //for each edge, let's add the neighbor (if it was not processed yet...)
    const int numEdgesBp = bp->boundaryPolygonOtherSideEdge.size();
    for(int i=0;i<numEdgesBp;i++) 
      if(bp->boundaryPolygonOtherSideEdge[i]!=NULL) { //is there anything on the other side? (anything belonging to this triangle?)
        ObjectId objectThisMesh = bp->getPolyhedronWherePolygonIs();
        pair<ObjectId,ObjectId> &objectsBoundedByThisEdge = bp->objectsOtherMeshBoundedByThisEdge[i];

        ObjectId objectOtherPolygon;
        if(objectsBoundedByThisEdge.first==DONT_KNOW_ID) {
          objectOtherPolygon = objectThisMesh; //the edge is not an edge from intersection --> both polygons are in the same object...
        } else {
          objectOtherPolygon = (objectThisMesh==objectsBoundedByThisEdge.first)?objectsBoundedByThisEdge.second:objectsBoundedByThisEdge.first;
          assert(objectOtherPolygon!=objectThisMesh); //the edge should separate different objects...
        }
        BoundaryPolygon *bpOtherSide = bp->boundaryPolygonOtherSideEdge[i];
        if(bpOtherSide->getPolyhedronWherePolygonIs()==DONT_KNOW_ID) { //not inserted into the stack yet...
          bpOtherSide->setPolyhedronWherePolygonIs(objectOtherPolygon);
          bpToProcess.push(bpOtherSide); //let's process it latter (DFS...)
        } else {
          //already labeled... should be already in the stack (or have already been processed...)
          ObjectId objAlreadyOtherSide = bpOtherSide->getPolyhedronWherePolygonIs();
          if(objAlreadyOtherSide!=objectOtherPolygon) {
            cerr << "Processing triangle from meshid: " << tri.first->getMeshId() << endl;
            cerr << "Objects bounded by this edge: " << objectsBoundedByThisEdge.first << " " << objectsBoundedByThisEdge.second << endl;
            cerr << "Already other side, should be, this side: " << objAlreadyOtherSide << " " << objectOtherPolygon << " " << objectThisMesh << endl << endl;
          }
          assert(bpOtherSide->getPolyhedronWherePolygonIs()==objectOtherPolygon);
        }
        
      }
    
  }

  //cerr << "Checking if everything ok... " << endl;
  int ctDontKnow = 0;
  for(BoundaryPolygon&bp:tri.second) {
    if(bp.getPolyhedronWherePolygonIs()==DONT_KNOW_ID || bp.getPolyhedronWherePolygonIs()==DONT_KNOW_FLAG) { //do we know the location of this polygon?
      ctDontKnow++;
    }    
  }
  assert(ctDontKnow==0); //the DFS should find all polygons and label them...

}

void locateTrianglesAndPolygonsInOtherMesh(const Nested3DGridWrapper *uniformGrid, 
                                            MeshIntersectionGeometry &geometry, 
                                            const unordered_set<const InputTriangle *> trianglesThatIntersect[2],
                                            vector< pair<const InputTriangle *,vector<BoundaryPolygon>> > polygonsFromRetesselationOfEachTriangle[2],
                                            int meshId,
                                            vector<ObjectId> &locationOfEachNonIntersectingTrianglesInOtherMesh) {  
    
    vector<InputTriangle> *inputTriangles = geometry.inputTriangles;

    timespec t0,t1,t0Function;

    
    clock_gettime(CLOCK_REALTIME, &t0);
    t0Function  = t0;

    

    clock_gettime(CLOCK_REALTIME, &t0);
    
    int numInputVerticesCoordinatesThisMesh = geometry.getNumVertices(meshId);

    vector<ObjectId> locationEachVertex(numInputVerticesCoordinatesThisMesh,DONT_KNOW_FLAG);
    vector<int> connectedComponentEachVertex(numInputVerticesCoordinatesThisMesh,DONT_KNOW_FLAG);
    vector<int> sampleVertexIdFromEachConnectedComponent;
    int numComponents;

{
    vector<vector<int> > adjList(numInputVerticesCoordinatesThisMesh);

    
    for(const InputTriangle&t:inputTriangles[meshId]) {
      if(trianglesThatIntersect[meshId].count(&t)==0) { //this triangle does not intersect the other mesh...
        adjList[t.getInputVertex(0)->getId()].push_back(t.getInputVertex(1)->getId());
        adjList[t.getInputVertex(0)->getId()].push_back(t.getInputVertex(2)->getId()); 

        adjList[t.getInputVertex(1)->getId()].push_back(t.getInputVertex(0)->getId());
        adjList[t.getInputVertex(1)->getId()].push_back(t.getInputVertex(2)->getId());

        adjList[t.getInputVertex(2)->getId()].push_back(t.getInputVertex(0)->getId());
        adjList[t.getInputVertex(2)->getId()].push_back(t.getInputVertex(1)->getId());        
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


    vector<InputVertex> verticesToLocateInOtherMesh;
    verticesToLocateInOtherMesh.reserve(numComponents);

    for(int i=0;i<numComponents;i++) {      
      int vertexToLocateId = sampleVertexIdFromEachConnectedComponent[i];

      verticesToLocateInOtherMesh.push_back(InputVertex(meshId,vertexToLocateId));           
    }
      
    int posStartVerticesOfIntersectingTrianglesInThisMesh = verticesToLocateInOtherMesh.size();

    int numPolygonsFromRetesselation = 0;
    int numTriFromRetesselation = 0;
    for(const auto &tri:polygonsFromRetesselationOfEachTriangle[meshId]) {
      const auto &boundaryPolygons = tri.second;
      numPolygonsFromRetesselation += boundaryPolygons.size();
      for(const BoundaryPolygon &b:boundaryPolygons)
        numTriFromRetesselation += b.triangulatedPolygon.size();
    }

    cerr << "Mesh " << meshId << " Num boundary polygons : " << numPolygonsFromRetesselation << endl;         
    cerr << "Mesh " << meshId << " Num tri from retesselation: " << numTriFromRetesselation << endl;
    

 
    
    clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "Total time to select vertices to locate: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
    clock_gettime(CLOCK_REALTIME, &t0);


    //First let's locate the triangles/polygons that contain at least one non-shared vertex.
    
    
    vector<ObjectId> locationOfEachVertexInOtherMesh(verticesToLocateInOtherMesh.size());
    //vertices of mesh "meshId" will be located in mesh "1-meshId"

    clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "Total time to compute center and get grid cell of center of intersecting triangles: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;


    //Locate the vertices...    

    PinMesh pointLocationAlgorithm(uniformGrid,&geometry);
    clock_gettime(CLOCK_REALTIME, &t0);
    pointLocationAlgorithm.locateVerticesInObject(verticesToLocateInOtherMesh,locationOfEachVertexInOtherMesh,1-meshId);
    clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "Total time to locate: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;


    //----------------------------------------------------------------------------------------------------------
    //now, we know in what object of the other mesh each triangle that does not intersect other triangles is...

    clock_gettime(CLOCK_REALTIME, &t0);


    const int numInputTrianglesThisMesh = inputTriangles[meshId].size();

    locationOfEachNonIntersectingTrianglesInOtherMesh = vector<ObjectId>(numInputTrianglesThisMesh,DONT_KNOW_FLAG);

    int ct =0;
    for(int tid = 0; tid < numInputTrianglesThisMesh;tid++) {
      const InputTriangle &t = inputTriangles[meshId][tid];
      if(trianglesThatIntersect[meshId].count(&t)==0) {
        int connectedComponentOfThisVertex = connectedComponentEachVertex[t.getInputVertex(0)->getId()];
        locationOfEachNonIntersectingTrianglesInOtherMesh[tid] = locationOfEachVertexInOtherMesh[connectedComponentOfThisVertex];
      }
    }

    //next step: locate triangles intersecting other triangles...

    int ctOnlySharedVertices = 0;
    //locationOfPolygonsFromRetesselationInOtherMesh = vector<ObjectId>(numPolygonsFromRetesselation,DONT_KNOW_FLAG);


    clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "Total copy triangle labels: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
    clock_gettime(CLOCK_REALTIME, &t0);

    const int numPolygonsFromRetesselationToProcess = polygonsFromRetesselationOfEachTriangle[meshId].size();
    #pragma omp parallel for
    for(int i =0;i<numPolygonsFromRetesselationToProcess;i++) {
      auto &tri= polygonsFromRetesselationOfEachTriangle[meshId][i];
      auto &boundaryPolygons = tri.second;

      bool locatedAllPolygons = true;
      for(BoundaryPolygon &polygon:boundaryPolygons) {
        const Vertex *nonSharedVertex = getNonSharedVertextFromPolygon(&(*polygon.vertexSequence.begin()),&(*polygon.vertexSequence.end()));
        if(nonSharedVertex!=NULL) {
          int connectedComponentOfThisVertex = connectedComponentEachVertex[nonSharedVertex->getId()];
          polygon.setPolyhedronWherePolygonIs(locationOfEachVertexInOtherMesh[connectedComponentOfThisVertex]);          
        } else {
          locatedAllPolygons = false;
          ctOnlySharedVertices++;
        }
      }
      if(!locatedAllPolygons) {
        locatePolygonsOtherMeshUsingDFS(geometry,tri);
      }
    }

    clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "Total locate polygons using DFS: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
    
    cerr << "Num polygons with only shared vertices: " << ctOnlySharedVertices << endl;
    cerr << "Num polygons with input vertices: " << numPolygonsFromRetesselation-ctOnlySharedVertices << endl;
     
    


    
    clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "Total entire location functions: " << convertTimeMsecs(diff(t0Function,t1))/1000 << endl;


    
}







//----------------------------------------------------------------------------
// TODO: not all vertices from the retesselated triangles will be in the output (?)... we do not need to store them!
//Each vector represents the vertices of a layer

//Each vector represents a set of objects in the same layer
//The objects are represented by a set of triangles (defining their boundaries)

double classifyTrianglesAndGenerateOutput(const Nested3DGridWrapper *uniformGrid, 
                                        MeshIntersectionGeometry &geometry, 
                                        const unordered_set<const InputTriangle *> trianglesThatIntersect[2],
                                        vector< pair<const InputTriangle *,vector<BoundaryPolygon>> > polygonsFromRetesselationOfEachTriangle[2],                                                                             
                                        ostream &outputStream) {
	timespec t0,t0ThisFunction,t1;
	clock_gettime(CLOCK_REALTIME, &t0);
	t0ThisFunction = t0;

	//int ctIntersectingTrianglesTotal =0;
  vector<InputTriangle> *inputTriangles = geometry.inputTriangles;

	vector<InputTriangle> outputTriangles[2]; //output triangles generated from non retesselated input triangles
	vector<RetesselationTriangle> outputTrianglesFromRetesselation[2];
	for(int meshId=0;meshId<2;meshId++){
    timespec t0,t1;

    vector<ObjectId> locationOfEachNonIntersectingTrianglesInOtherMesh;

    clock_gettime(CLOCK_REALTIME, &t0);

    locateTrianglesAndPolygonsInOtherMesh(uniformGrid,geometry,trianglesThatIntersect,polygonsFromRetesselationOfEachTriangle,meshId,locationOfEachNonIntersectingTrianglesInOtherMesh);
     //TODO: remove this....
    

    clock_gettime(CLOCK_REALTIME, &t1);

    cerr << "# Time to locate vertices in other mesh: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
    clock_gettime(CLOCK_REALTIME, &t0);

		const int numTrianglesThisMesh = inputTriangles[meshId].size();
		int ctTrianglesProcessed = 0;

		for(int i=0;i<numTrianglesThisMesh;i++) {
			const InputTriangle &t=inputTriangles[meshId][i];
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

   
    //for each polygon   
    bool hadTriangleDontKnow = false;
    int numTrianglesRetesselated = polygonsFromRetesselationOfEachTriangle[meshId].size();    
    for(int i=0;i<numTrianglesRetesselated;i++) {
      //for each boundary polygon from this triangle...
      for(const BoundaryPolygon&p:polygonsFromRetesselationOfEachTriangle[meshId][i].second) {
        ObjectId objWhereTriangleIs = p.getPolyhedronWherePolygonIs(); //locationOfTrianglesFromRetesselationInTheOtherMesh[i];//locationOfEachVertexInOtherMesh[ctTrianglesProcessed++];          
        
        if (objWhereTriangleIs!=OUTSIDE_OBJECT) {   
          for(array<const Vertex *,3> tris:p.triangulatedPolygon) {
            assert(tris[0]->getMeshId()<3);
            outputTrianglesFromRetesselation[meshId].push_back( RetesselationTriangle(*tris[0], *tris[1], *tris[2],p.above, p.below));
          }
        } 
        if(objWhereTriangleIs==DONT_KNOW_FLAG || objWhereTriangleIs==DONT_KNOW_ID) {
          hadTriangleDontKnow = true;
        }   
      }           
    }

    assert(!hadTriangleDontKnow); //this shouldn't happen...


    clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "# Time to classify triangles vertices in other mesh: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
    clock_gettime(CLOCK_REALTIME, &t0);
	}
	clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to locate vertices and classify triangles: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  
  clock_gettime(CLOCK_REALTIME, &t0); 


  #ifdef DEBUGGING_MODE
  for(int meshId=0;meshId<2;meshId++){
    const int szOutputTrianglesFromIntersection = outputTrianglesFromRetesselation[meshId].size();
    for(int i=0;i<szOutputTrianglesFromIntersection;i++) {
      const RetesselationTriangle &t = outputTrianglesFromRetesselation[meshId][i]; 
     // cerr << t.getVertex(0)->getMeshId() << " " << t.getVertex(1)->getMeshId() << endl;
      assert(t.getVertex(0)->getMeshId()<3);
    }  
  }
  #endif

  //Currently, let's write all vertices to the output
	int numVerticesInEachMesh[3] = {geometry.getNumVertices(0),geometry.getNumVertices(1),geometry.getNumVertices(2)};

  int baseIdOfVertices[3];
  baseIdOfVertices[0] = 0; //vertices originally from the first mesh will be counted starting on 1
  baseIdOfVertices[1] = baseIdOfVertices[0] + numVerticesInEachMesh[0]; 
  baseIdOfVertices[2] = baseIdOfVertices[1] + numVerticesInEachMesh[1]; 

	//compute the new ids of the shared vertices...



  int totalNumberOutputVertices = numVerticesInEachMesh[0] + numVerticesInEachMesh[1] + numVerticesInEachMesh[2];
  int totalNumberOutputVerticesFromNonIntersectingTriangles = numVerticesInEachMesh[0] + numVerticesInEachMesh[1] ;
  int totalNumberOutputTriangles = outputTriangles[0].size() + outputTriangles[1].size() + outputTrianglesFromRetesselation[0].size() + outputTrianglesFromRetesselation[1].size();

  
  unordered_map<pair<int,int>,int> edgesIds; //maybe use unordered_map for performance (if necessary...)
  vector<pair<int,int> > outputEdges;
 
	//vector<pair<int,int> > outputEdgesWithRepetition[2];

	for(int meshId=0;meshId<2;meshId++) {
		int numNewEdgesToAdd =0;
		#pragma omp parallel
		{			
			
			vector<pair<int,int> > myEdgesFound;	
			
      
      const int sz =  outputTriangles[meshId].size();
			#pragma omp for
		  for(int i=0;i<sz;i++) {
		  	const Triangle &t = outputTriangles[meshId][i];
        int a = t.getVertex(0)->getId() + baseIdOfVertices[meshId];
        int b = t.getVertex(1)->getId() + baseIdOfVertices[meshId];
        int c = t.getVertex(2)->getId() + baseIdOfVertices[meshId];

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
		  	const RetesselationTriangle &t = outputTrianglesFromRetesselation[meshId][i];
				

        int a = t.getVertex(0)->getId() + baseIdOfVertices[t.getVertex(0)->getMeshId()];
        int b = t.getVertex(1)->getId() + baseIdOfVertices[t.getVertex(1)->getMeshId()];
        int c = t.getVertex(2)->getId() + baseIdOfVertices[t.getVertex(2)->getMeshId()];


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
	__gnu_parallel::sort(outputEdges.begin(),outputEdges.end());
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

  
  double timeWithoutIO =  convertTimeMsecs(diff(t0ThisFunction,t1))/1000;
  cerr << "Total time before writing output: " << timeWithoutIO << endl;
 

  //now, let's write everything in the output!
  outputStream << "OFF\n";
  //outputStream << totalNumberOutputVertices << " " << totalNumberOutputEdges << " " << totalNumberOutputTriangles << '\n';
  outputStream << totalNumberOutputVertices << " " <<  totalNumberOutputTriangles << " 0\n";

  //print the coordinates of the vertices...
  
  geometry.storeAllVertices(outputStream);

	//print edges...
	//for(const pair<int,int> &p:outputEdges) {
	//	outputStream << p.first+1 << " " << p.second+1 << "\n"; //in a GTS file we start counting from 1...
	//}

	//print triangles...
	for(int meshId=0;meshId<2;meshId++) 
		for(Triangle &t : outputTriangles[meshId]) {
      int a = t.getVertex(0)->getId() + baseIdOfVertices[meshId];
      int b = t.getVertex(1)->getId() + baseIdOfVertices[meshId];
      int c = t.getVertex(2)->getId() + baseIdOfVertices[meshId];


			if(t.above != OUTSIDE_OBJECT) {
				//according to the right hand rule, the ordering of the vertices should be (c,b,a)
				swap(a,c);
			} 
      outputStream << "3 " << a << " " << b << " " << c << "\n";
      /*
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
			outputStream << mapEdgesIds2[e.first][e.second]+1 << "\n";*/
		}

    
		for(int meshId=0;meshId<2;meshId++) 
			for(RetesselationTriangle &t : outputTrianglesFromRetesselation[meshId]) {
				int a = t.getVertex(0)->getId() + baseIdOfVertices[t.getVertex(0)->getMeshId()];
				int b = t.getVertex(1)->getId() + baseIdOfVertices[t.getVertex(1)->getMeshId()];
				int c = t.getVertex(2)->getId() + baseIdOfVertices[t.getVertex(2)->getMeshId()];

				

				if(t.above != OUTSIDE_OBJECT) {
					//according to the right hand rule, the ordering of the vertices should be (c,b,a)
					swap(a,c);
				} 

        outputStream << "3 " << a << " " << b << " " << c << "\n";

        /*
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
				outputStream << mapEdgesIds2[e.first][e.second]+1 << "\n";*/
			}
    
	
		cerr << "Output vertices         : " << totalNumberOutputVertices << endl;
		cerr << "Output edges            : " << totalNumberOutputEdges << endl;
		cerr << "Output triangles non int: " << totalNumberOutputTriangles << endl;
		//cerr << "Intersecting triangles  : " << ctIntersectingTrianglesTotal << endl;


    return timeWithoutIO;

}