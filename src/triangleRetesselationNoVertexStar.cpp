#include "triangleRetesselation.h"

long long ctEdgeIntersect  = 0; //how many of the intersection tests performed are true?
long long ctEdgeDoNotIntersect = 0;
int numVerticesInEdges = 0;
//stores some counters (statistics)
struct StatisticsAboutRetesseation {
  StatisticsAboutRetesseation() {
    ctTrianglesRetesselate = 0;
    ctEdgesActuallyInsertedInRetesselation =0;
    ctEdgesInsertedBecauseOfConstraint=0;
    ctVerticesInsertedTrianglesToRetesselate=0;
    ctEdgesTestedToInsertInRetesselation=0;
    simple = difficult = 0;
    ctConvexPolygonsInTriangleRetesselations  = ctConcavePolygonsInTriangleRetesselations = 0;
    numConnectedPolygonSubdivisionOfTriangles = numDisconnectedPolygonSubdivisionOfTriangles = 0;
  }

  void printStatsSoFar() {
    cerr << numVerticesInEdges << endl;
    cerr << "Num triangles to retesselate                                 : "  << ctTrianglesRetesselate << "\n";
    cerr << "Total number of edges we tried to insert during retesselation: "  << ctEdgesTestedToInsertInRetesselation << "\n";  
    cerr << "Total edges tried insert and actually inserted               : " << ctEdgesActuallyInsertedInRetesselation  <<"\n";
    cerr << "Total edges inserted because are constraint edges            : " << ctEdgesInsertedBecauseOfConstraint <<"\n";
    cerr << "Total number of edges in retesselated triangles (add prev. 2): " << ctEdgesInsertedBecauseOfConstraint+ctEdgesActuallyInsertedInRetesselation <<"\n";
    cerr << "Total number of vertices in retesselated triangles           : "  << ctVerticesInsertedTrianglesToRetesselate << "\n";
    cerr << "Simple case triangles                                        : " << simple <<"\n";
    cerr << "Difficult case triangles                                     : " << difficult << "\n";
    cerr << "Number of convex polygons in triangle retesselations         : " << ctConvexPolygonsInTriangleRetesselations << "\n";
    cerr << "Number of concave polygons in triangle retesselations        : " << ctConcavePolygonsInTriangleRetesselations << "\n";
    cerr << "Number of connected polygon subdivision in triangles         : " << numConnectedPolygonSubdivisionOfTriangles << "\n";
    cerr << "Number of disconnected polygon subdivision in triangles      : " << numDisconnectedPolygonSubdivisionOfTriangles << "\n";
  }
  long long simple, difficult;
  long long ctTrianglesRetesselate;
  long long ctEdgesActuallyInsertedInRetesselation;
  long long ctEdgesInsertedBecauseOfConstraint;
  long long ctVerticesInsertedTrianglesToRetesselate;
  long long ctEdgesTestedToInsertInRetesselation;
  long long ctConvexPolygonsInTriangleRetesselations;
  long long ctConcavePolygonsInTriangleRetesselations;
  long long numDisconnectedPolygonSubdivisionOfTriangles;
  long long numConnectedPolygonSubdivisionOfTriangles;
};

//temporary variables... we avoid reallocating them every time this function is calles
struct TempVarsRetesselateTriangleFunction {
  vector<const Vertex *> polygons;
  vector<array<const Vertex *,3> > wedgesTemp;
  vector<pair<const Vertex *,const Vertex *> > raggedArraySortedEdges;
  vector<bool> usedWedgesTemp;
  VertCoord tempVertCoords[8];
  MeshIntersectionGeometry::TempVarsIsCloser tempVarsIsCloser;
  MeshIntersectionGeometry::TempVarsIsAngleWith0Greater tempVarsIsAngleWith0Greater;
  MeshIntersectionGeometry::TempVarsGetPlaneTriangleIsNotPerpendicular tempVarsGetPlaneTriangleIsNotPerpendicular;
  MeshIntersectionGeometry::TempVarsIsBoundaryClockwise tempVarsIsBoundaryClockwise;
  MeshIntersectionGeometry::TempVarsIsTriangleClockwisedOriented tempVarsIsTriangleClockwisedOriented;
};


#define ERROR_CODE -2123123123
int binarySearchFirstOccurrency(const vector<pair<const Vertex *,const Vertex *> > &v,const Vertex * key) {
  int mid;
  int lo = 0;
  int hi = v.size()-1;
  while (lo < hi) {
      mid = lo + (hi-lo)/2;   // note: division truncates
      if (v[mid].first>=key)
         hi = mid;
      else
         lo = mid+1;
  }
          
  if (v[lo].first!=key)
      return ERROR_CODE;                // p(x) is true for all x in S!
      
  return lo;         // lo is the least x for which p(x) is true
}

bool checkIfPolygonsAreConnected(const vector<pair<const Vertex *,const Vertex *> > &raggedArraySortedEdges) {
  queue<const Vertex *> verticesToProcess;
  set<const Vertex *> processedVertices;
  verticesToProcess.push(raggedArraySortedEdges[0].first);
  processedVertices.insert(raggedArraySortedEdges[0].first);

  
  

  //visit all vertices reachable from v...
  while(!verticesToProcess.empty()) {
    const Vertex * v = verticesToProcess.front();
    verticesToProcess.pop();

    int posStartVRaggedArray =  binarySearchFirstOccurrency(raggedArraySortedEdges,v);
    //cerr << "v and pos v " << v << " " << posStartVRaggedArray << endl;
    for(int i=posStartVRaggedArray;i<raggedArraySortedEdges.size()&&raggedArraySortedEdges[i].first==v;i++) {
      const Vertex * adj = raggedArraySortedEdges[i].second; //adj is one of the vertices adjacent to v
      if(processedVertices.count(adj)==0) {
        processedVertices.insert(adj);
        verticesToProcess.push(adj);
      }
    }
  }

  //cerr << "Checking connectivity... " << endl;
  for(const pair<const Vertex *,const Vertex *> &p:raggedArraySortedEdges) {
    if(processedVertices.count(p.first)==0) return false; //at least one vertex was not visited...
    //cerr << p.first << " " << p.second << endl;

  } 
 /* cerr << "Visited: " << endl;
  for(int i:processedVertices) cerr << i << endl;  
  */

  return true;
}



//returns the index i of the vector v, where v[i][0] = key.first and v[i][1] = key.second 
int binarySearch(vector<array<const Vertex *,3> > &v,const pair<const Vertex *,const Vertex *> &key) {
  int lo = 0, hi = v.size()-1, mid;
  while (lo<=hi) {
    mid = lo + (hi-lo)/2;
    if(v[mid][0]==key.first && v[mid][1]==key.second) 
      return mid;
    else if ( (v[mid][0] < key.first) || ((v[mid][0] == key.first) && (v[mid][1] < key.second)) ) //is the element v[mid] smaller than the key?
      lo = mid+1;
    else 
      hi = mid-1;
  }
  return -1;
}






//given a triangle t, this function will split t at the intersection edges (intersection with other triangles)
//and retesselate t, creating more triangles.
//each edge is represented by a pair of ids of the vertices connected by that edge
//t is retesselated by:
//- using an algorithm that sorts the edges creating wedges to extract the polygons (whose union is t) that need to be retesselated
//- it is faster to retesselate these polygons individually than to retesselate the whole t (this is true mainly when there are several edges from intersection in t)

//Output: polygons
//-- Polygons will be a list of vertices ids representing the polygons
//-- Example: v1,v2,v3,v1,v9,v12,v10,v9 represents two polygons: v1,v2,v3,v1 and v9,v12,v10.v9

//raggedArraySortedEdges,usedWedgesTemp and wedgesTemp are temporary vectors (we pass them to the function to avoid too much memory reallocation)
//The first numEdgesFromIntersection edges in edges are from the intersections
//The last ones are edges from the boundary of the original rectangle being retesselated

//raggedArraySortedEdges will be filled with the edges in the triangle (if we have an edge (a,b) , the array will have both (a,b) and (b,a)). Also,
//at the end of the function the edges in raggedArraySortedEdges will be sorted by the first vertex and, if there is a tie, by the slope of the edge (w.r.t. the projection to whatPlaneProjectTo)
//at the end of this function raggedArraySortedEdges will contain the edges (two directed edges for each edge) sorted
//by the first vertex and by the angle of the edge 
void sortEdgesAndExtractPolygonsFromEdgeListUsingWedges(const MeshIntersectionGeometry &meshIntersectionGeometry, const vector<pair<const Vertex *,const Vertex *> > &edges,const int numEdgesFromIntersection, const int planeProjectTriangleTo,  TempVarsRetesselateTriangleFunction &tempVars ) {
  //Algorithm: http://ac.els-cdn.com/016786559390104L/1-s2.0-016786559390104L-main.pdf?_tid=7586b72e-a059-11e6-afdd-00000aacb35d&acdnat=1478021856_7213cc56dd2148587a06891102f5d518
   
  
  vector<array<const Vertex *,3> > &wedgesTemp = tempVars.wedgesTemp;
  vector<pair<const Vertex *,const Vertex *> > &raggedArraySortedEdges = tempVars.raggedArraySortedEdges;
  vector<bool> &usedWedgesTemp = tempVars.usedWedgesTemp;
  vector<const Vertex *> &polygons = tempVars.polygons;

  //First, we need to find all the wedges:

  //Duplicate each undirected edge such that we will have two directed edges (u,v) and (v,u)
  //sort the edges basing first on the first vertex and second on the angle the edge makes with the horizontal line (supposing the triangle is projected to planeProjectTriangleTo)
  //Process the groups in the sorted list:
  //-- A Group is a set of edges with the same first vertex
  //-- Within each group, combine each pair of consecutive edges (including the last with the first) to form a wedge: (a,b) (a,c) --> (c,a,b)
  const int numDirectedEdges = (edges.size())*2 ;
  raggedArraySortedEdges.resize(numDirectedEdges);
  const int numOrigEdges = edges.size();

  for(int i=0;i<numOrigEdges;i++) { //duplicate the edges
    raggedArraySortedEdges[2*i] = edges[i]; //add (u,v)
    raggedArraySortedEdges[2*i+1].first = edges[i].second; //add (v,u)
    raggedArraySortedEdges[2*i+1].second = edges[i].first;
  }
  
  int x = 0;
  //now we need to sort the edges basing on the first vertex and, then, on the angle with horizon (considering whatPlaneProjectTo)
  
  sort(raggedArraySortedEdges.begin(),raggedArraySortedEdges.end(), [&](const pair<const Vertex *,const Vertex *> &e1, const pair<const Vertex *,const Vertex *> &e2) {
            if(e1.first!=e2.first) return e1.first<e2.first; //sort first by the first vertex..
            return meshIntersectionGeometry.isAngleWith0Greater(*e1.first, *e1.second, *e2.second, planeProjectTriangleTo, tempVars.tempVarsIsAngleWith0Greater);
        } );
  

  wedgesTemp.resize(0);
  int numWedgesFound = 0;
  //Now let's process the groups of edges (set of edges with the same first vertex..) to extract the wedges  
  for(int firstElementGroup =0;firstElementGroup<numDirectedEdges;) {
    for(int currElemnt = firstElementGroup;currElemnt<numDirectedEdges;currElemnt++) {
      //combine the currentElement with the next one to form a wedge...
      int nextElement = currElemnt+1;
      //(a,b)         (a,c) --> (c,a,b)
      //currElement    nextElement

      wedgesTemp.resize(numWedgesFound+1);
      numWedgesFound++;

      array<const Vertex *,3> &wedgeToAdd = wedgesTemp.back();
      //did we reach the end of a group?
      if(nextElement >= numDirectedEdges || raggedArraySortedEdges[nextElement].first != raggedArraySortedEdges[currElemnt].first) {
        //if yes, we will create a wedge from currElement to the first one! (to complete the wedges from this group..)
        //the next element is actually the first element in this group...
        wedgeToAdd[0] = raggedArraySortedEdges[firstElementGroup].second; //c
        wedgeToAdd[1] = raggedArraySortedEdges[firstElementGroup].first; //a
        wedgeToAdd[2] = raggedArraySortedEdges[currElemnt].second;  //b

        firstElementGroup = currElemnt+1; //let's process the next group!
        break;
      } else {
        //else, we will create a wedge from the current to the next...
        wedgeToAdd[0] = raggedArraySortedEdges[nextElement].second; //c
        wedgeToAdd[1] = raggedArraySortedEdges[nextElement].first; //a
        wedgeToAdd[2] = raggedArraySortedEdges[currElemnt].second;  //b
      }
    }
  }

  //Now, wedgesTemp contain all the wedges we need to extract the polygons!

  //Once we have the wedges, we need to group the wedges to form the polygons:
  //Sort the list of wedges by the first and, then, second parameters
  //For each unused wedge w=(a,b,c):
  //- Mark w as used
  //- Start a new region with w
  //- Binary search the list to find w2 = (x1,x2,x3)
  //- Assert the binary search does not fail (it shouldn't!)
  //- If x2==a and x3==b (i.e., we returned to w): stop... we finished constructing a polygon...
  //- Go to the first step...

  //sort wedges:
  sort(wedgesTemp.begin(),wedgesTemp.end());
  const int numWedges = wedgesTemp.size();
  usedWedgesTemp.resize(numWedges);
  for(int i=0;i<numWedges;i++) usedWedgesTemp[i] = false;

  //Let's use a ragged array for performance...
  polygons.resize(0);  



  for(int wedgeStart=0;wedgeStart<numWedges;wedgeStart++) {
    if(usedWedgesTemp[wedgeStart]) continue;
    usedWedgesTemp[wedgeStart] = true; //mark w as used...
    polygons.push_back(wedgesTemp[wedgeStart][0]);//first vertex of the polygon...

    pair<const Vertex *,const Vertex *> keyForBinarySearch(wedgesTemp[wedgeStart][1],wedgesTemp[wedgeStart][2]);
    while(true) {
      int next = binarySearch(wedgesTemp,keyForBinarySearch);// binary search the list to find w2 = (x1,x2,x3) where (x1,x2) is the keyForBinarySearch
      //assert...
      assert(next!=-1); //this shouldn't happen..

      usedWedgesTemp[next] = true;
      polygons.push_back(wedgesTemp[next][0]);
      keyForBinarySearch.first = wedgesTemp[next][1];
      keyForBinarySearch.second = wedgesTemp[next][2];

      if(wedgesTemp[next][1] == wedgesTemp[wedgeStart][0]) { //have I completed a polygon?
        polygons.push_back(wedgesTemp[wedgeStart][0]);
        break; //let's find a (possible) next vertex to start a new polygon...
      }
    }
  }

  
}



int doesPolygonContainAllEdges(const vector<const Vertex *> &polygons,const int firstElementPolygon,const int lastElementPolygon,
                                const set<pair<const Vertex *,const Vertex *> > & edgesFromTriangleBoundary) {

  if(edgesFromTriangleBoundary.count(make_pair(polygons[firstElementPolygon],polygons[firstElementPolygon+1]))) {
    //check if all edges in same order are there...
    for(int i=firstElementPolygon+1;i<lastElementPolygon;i++) {
      if(edgesFromTriangleBoundary.count(make_pair(polygons[i],polygons[i+1]))==0) return 0;
    } 
    return 1;
  }
  else 
    if(edgesFromTriangleBoundary.count(make_pair(polygons[firstElementPolygon+1],polygons[firstElementPolygon]))) {
      //check if all edges in same order are there...
      for(int i=firstElementPolygon+1;i<lastElementPolygon;i++) {
        if(edgesFromTriangleBoundary.count(make_pair(polygons[i+1],polygons[i]))==0) return 0;
      } 
      return -1;
    }
    else return 0;
}





int ctIt = 0;



//tempVars should have at least 8 slots
//the new polygons generated from the retesselation will be added to the back of "newPolygonsGeneratedFromRetesselation"
void retesselateTriangleUsingWedgeSorting(MeshIntersectionGeometry & meshIntersectionGeometry, 
                                          const vector< pair<VertexFromIntersection, VertexFromIntersection> >  &edgesFromIntersection,
                                          const vector<int> &edgesFromIntersectionThisTriangle,
                                          const InputTriangle &t, 
                                          vector<BoundaryPolygon> &newPolygonsGeneratedFromRetesselation,
                                          TempVarsRetesselateTriangleFunction &tempVars, 
                                          StatisticsAboutRetesseation &statistics) {
  /*
  //First we need to determine, for each vertex v, the edges incident in v
  //each vertex will be represented twice (as two directed edges): (u,v) and (v,u)
  //The edges incident in each vertex will be sorted based on their slope (supposing they are projected to the plane "whatPlaneProjectTriangleTo")
  */

  //we need to "retesselate" t and orient all the new triangles properly
  //this set will store the edges we used in the retesselated triangle
  //we need to choose what edges to create and, then, use these new edges to reconstruct the triangulation..
  //This set will have the edges that will form the retriangulation of ts..
  vector<pair<const Vertex *,const Vertex *> > edgesUsedInThisTriangle; 

  int ct =0;
  //The edges from intersection will, necessarelly, be in the triangulation...
  for(int i:edgesFromIntersectionThisTriangle) {
    const auto &edgeFromIntersection =  edgesFromIntersection[i];
    //cerr << edges[edgeId].first << " " << edges[edgeId].second << endl;
    //assert(edgeFromIntersection.first != edgeFromIntersection.second);
    edgesUsedInThisTriangle.push_back(make_pair(&edgeFromIntersection.first, &edgeFromIntersection.second));   

    #ifdef COLLECT_STATISTICS
      #pragma omp atomic
      statistics.ctEdgesInsertedBecauseOfConstraint++;   
    #endif
  }


  //Now, we need to add more edges to fill the parts of ts that still do not form triangle
  vector<const Vertex *> verticesToTesselate;
  //what vertices will we have in the new triangulation of ts? (the original vertices + the vertices of the edges from intersections)
  for(const auto &p:edgesUsedInThisTriangle) {
    verticesToTesselate.push_back(p.first);
    verticesToTesselate.push_back(p.second);
  }

  sort(verticesToTesselate.begin(),verticesToTesselate.end(),[](const Vertex *a, const Vertex *b){
                                return (*a)<(*b);
                              });
  auto newEnd = unique(verticesToTesselate.begin(),verticesToTesselate.end(),[](const Vertex *a, const Vertex *b){
                                return (*a)==(*b);
                              });
  verticesToTesselate.resize(newEnd-verticesToTesselate.begin());

  //Initially we will insert only the vertices from the intersection (the original vertices of the triangle will be inserted later)

    
  const int numVerticesToTesselateNow = verticesToTesselate.size();
  const InputVertex* tv0 = t.getInputVertex(0);
  const InputVertex* tv1 = t.getInputVertex(1);
  const InputVertex* tv2 = t.getInputVertex(2);
  const int meshContainingT = t.getMeshId();

  //Computes the vertices incident to each edge of the original triangle...
  //verticesIncidentEachEdgeOriginalTriangle[0] is for vertex t.p[0]-t.p[1]
  //verticesIncidentEachEdgeOriginalTriangle[1] is for vertex t.p[1]-t.p[2]
  //verticesIncidentEachEdgeOriginalTriangle[2] is for vertex t.p[2]-t.p[0]
  vector<const Vertex *> verticesIncidentEachEdgeOriginalTriangle[3];
  for(int i=0;i<numVerticesToTesselateNow;i++) {
    const VertexFromIntersection &v = *((VertexFromIntersection*) verticesToTesselate[i]); //since we have not inserted the original vertices yet, all vertices are InputVertex...
    //for each vertex v that is not one of the original vertices of the triangle, let's see if v is in one of the edges of t...

    //to be incident to an edge, v should be the intersection of an edge of t with another triangle!
    //the other intersectionVertices are from the intersection of an edge of other triangles with t!
    //VertexFromIntersection are composed of an edge and a triangle. If the edge defining v is from this mesh 
    //, then we know that the edge is necessarely an edge of this triangle (t)
    if(v.getMeshOfEdgeDefiningVertex()!=meshContainingT) continue;

    //now we know v is defined as the intersection of an edge of t with a triangle of the other mesh

    if( (v.edge[0] == (*tv0) && v.edge[1]  == (*tv1)) || (v.edge[0]  == (*tv1) && v.edge[1]  == (*tv0)) ) {
      verticesIncidentEachEdgeOriginalTriangle[0].push_back(&v);
      continue;
    }   
    if( (v.edge[0]  == (*tv1) && v.edge[1]  == (*tv2)) || (v.edge[0]  == (*tv2) && v.edge[1]  == (*tv1)) ) {
      verticesIncidentEachEdgeOriginalTriangle[1].push_back(&v);
      continue;
    }
    if( (v.edge[0]  == (*tv2) && v.edge[1]  == (*tv0)) || (v.edge[0]  == (*tv0) && v.edge[1]  == (*tv2)) ) {
      verticesIncidentEachEdgeOriginalTriangle[2].push_back(&v);
      continue;
    }
    assert(false); //this should never be true...
  }

  verticesToTesselate.push_back(t.getVertex(0));
  verticesToTesselate.push_back(t.getVertex(1));
  verticesToTesselate.push_back(t.getVertex(2));
  int numVTriangle = verticesToTesselate.size();


  #ifdef COLLECT_STATISTICS
    #pragma omp atomic
    statistics.ctVerticesInsertedTrianglesToRetesselate += numVTriangle;
  #endif



  
  set<pair<const Vertex *,const Vertex *> > edgesFromTriangleBoundary;

  bool allSimple = true;
  //no vertex intersect the edge t.p[0]-t.p[1] --> it will be in the triangulation!!!
  for(int edge =0;edge<3;edge++) {
    const int v1Pos = edge;
    const int v2Pos = (edge+1)%3;

    const int verticesIncidentEdge = verticesIncidentEachEdgeOriginalTriangle[edge].size();
    if(verticesIncidentEdge==0) {
      const Vertex * v1 = t.getVertex(v1Pos);
      const Vertex * v2 = t.getVertex(v2Pos);
              
      pair<const Vertex *,const Vertex *> candidateEdge(v1,v2);  

      #ifdef COLLECT_STATISTICS
        #pragma omp atomic
        statistics.ctEdgesActuallyInsertedInRetesselation++;
      #endif

      assert(candidateEdge.first != candidateEdge.second);
      edgesUsedInThisTriangle.push_back(candidateEdge);
      edgesFromTriangleBoundary.insert(candidateEdge);
    } else {
      if(verticesIncidentEdge==1) {
        const Vertex * v1 = t.getVertex(v1Pos);
        const Vertex * v2 = t.getVertex(v2Pos);

        const Vertex * v =  *(verticesIncidentEachEdgeOriginalTriangle[edge].begin());          
        
        pair<const Vertex *,const Vertex *>  candidateEdge;    

        candidateEdge.first = v1; 
        candidateEdge.second = v;


        #ifdef COLLECT_STATISTICS
          #pragma omp atomic
          statistics.ctEdgesActuallyInsertedInRetesselation++;
        #endif

        assert(candidateEdge.first != candidateEdge.second);
        edgesUsedInThisTriangle.push_back(candidateEdge);
        edgesFromTriangleBoundary.insert(candidateEdge);

        candidateEdge.first = v; 
        candidateEdge.second = v2;        

        #ifdef COLLECT_STATISTICS
          #pragma omp atomic
          statistics.ctEdgesActuallyInsertedInRetesselation++;
        #endif
        assert(candidateEdge.first != candidateEdge.second);
        edgesUsedInThisTriangle.push_back(candidateEdge);
        edgesFromTriangleBoundary.insert(candidateEdge);
      } else {
        allSimple= false;
      }
    }
  }




  #ifdef COLLECT_STATISTICS        
      if(allSimple) {
        #pragma omp atomic
        statistics.simple++;
      }
      else {
        #pragma omp atomic
        statistics.difficult++;
      }
  #endif
  //Now we have a set of internal edges and possibly some boundary edges
  //let's split the boundary edges (with more than 2 vertices) at the intersection points
  for(int i=0;i<3;i++) {
    const InputVertex * edgeOrig = t.getInputVertex(i); //the two vertices of the edge i...
    const InputVertex * edgeDest = t.getInputVertex((i+1)%3); 

    vector<const Vertex *> &verticesToCreateEdges = verticesIncidentEachEdgeOriginalTriangle[i];
    //cerr << "Num of vertices in edge: " << i << " : " << verticesToCreateEdges.size() << endl;
    if(verticesToCreateEdges.size()<=1) continue; //edges were already inserted in the previous step...
    //here we need to deal only with the situation where 2 or more vertices are incident to an edge.

    //cerr << "Ids of candidate vertices: " << endl;
    //for(int v: verticesToCreateEdges) cerr << v << endl;

    //sort the vertices basing on the distance from v1

    sort(verticesToCreateEdges.begin(),verticesToCreateEdges.end(), [&](const Vertex *v1, const Vertex *v2) {
                              return meshIntersectionGeometry.isCloser(*edgeOrig,*v1,*v2,tempVars.tempVarsIsCloser);
                            });
    
    verticesToCreateEdges.push_back(edgeDest); //v2 is the vertex that is furthest from v1 along edge v1,v2 
    //we will add edges:
    //v1 to verticesToCreateEdges[0]
    //verticesToCreateEdges[0] to verticesToCreateEdges[1]
    //....... verticesToCreateEdges[n-2] to verticesToCreateEdges[n-1] (=v2, that was pushed_back above..)

    pair<const Vertex *,const Vertex *> candidateEdge;   
    candidateEdge.first = edgeOrig; 
    candidateEdge.second = verticesToCreateEdges[0];
    #ifdef COLLECT_STATISTICS
      #pragma omp atomic
      statistics.ctEdgesActuallyInsertedInRetesselation++;
    #endif
    //cerr << "Candidate edge: " << candidateEdge.first << " " << candidateEdge.second << endl;
    //cerr <<"V1 v2 : " << v1 << " " << v2 << endl;
    assert(candidateEdge.first != candidateEdge.second);
    //assert(edgesUsedInThisTriangle.count(candidateEdge)==0);
    //cerr << "Inserting edge: " << candidateEdge.first << " " << candidateEdge.second << endl;
    //printVertexForDebugging(*getPointFromVertexId(candidateEdge.first,meshWhereTriangleIs));
    //printVertexForDebugging(*getPointFromVertexId(candidateEdge.second,meshWhereTriangleIs));

    edgesUsedInThisTriangle.push_back(candidateEdge);
    edgesFromTriangleBoundary.insert(candidateEdge);

    const int numVerticesToCreateEdges = verticesToCreateEdges.size();
    for(int i=0;i<numVerticesToCreateEdges-1;i++) {
      candidateEdge.first = verticesToCreateEdges[i]; 
      candidateEdge.second = verticesToCreateEdges[i+1];      
      #ifdef COLLECT_STATISTICS
        #pragma omp atomic
        statistics.ctEdgesActuallyInsertedInRetesselation++;
      #endif
      assert(candidateEdge.first != candidateEdge.second);
      //assert(edgesUsedInThisTriangle.count(candidateEdge)==0);

      //cerr << "Inserting edge: " << candidateEdge.first << " " << candidateEdge.second << endl;
      //printVertexForDebugging(*getPointFromVertexId(candidateEdge.first,meshWhereTriangleIs));
      //printVertexForDebugging(*getPointFromVertexId(candidateEdge.second,meshWhereTriangleIs));
      edgesUsedInThisTriangle.push_back(candidateEdge);
      edgesFromTriangleBoundary.insert(candidateEdge);
    }
  }


  
  // (at least) three of the last edges added to the triangle are, necessarely, in the boundary of the triangle
  //they also have the same orientation of the triangle.
  //so, we can use any of them as the "seed edge"
  //we need a seed output edge that will be oriented in the same way t is oriented..
  //we will use this seed to reorient all the triangles resulting from the retesselation of t
  pair<const Vertex *,const Vertex *> seedEdge = edgesUsedInThisTriangle.back(); 

  //Now, edgesUsedInThisTriangle will contain the original edges of the triangle (split because of intersections) and new edges
  //created because of intersections with other triangles. Thus, edgesUsedInThisTriangle will define a planar partitioning..

 // cerr << "End..." << allSimple << endl;
  
  vector<pair<const Vertex *,const Vertex *> > edgesUsedInThisTriangleV(edgesUsedInThisTriangle.begin(),edgesUsedInThisTriangle.end());

  const int whatPlaneProjectTriangleTo = meshIntersectionGeometry.getPlaneTriangleIsNotPerpendicular(t,tempVars.tempVarsGetPlaneTriangleIsNotPerpendicular);


  sortEdgesAndExtractPolygonsFromEdgeListUsingWedges(meshIntersectionGeometry, edgesUsedInThisTriangleV,edgesFromIntersection.size(), whatPlaneProjectTriangleTo,  tempVars );

  vector<const Vertex *> &polygons = tempVars.polygons; 

  //at the end of the previous function raggedArrayEdges will contain the directed edges sorted by their first vertex...
  bool arePolygonsConnected = checkIfPolygonsAreConnected(tempVars.raggedArraySortedEdges);



  //at the end of the previous function, raggedArrayEdges will be filled with the edges sorted by the first vertex of each edge
  //(if there is a tie the edges are sorted by the slope w.r.t. the plane we projected the triangle)
  int numConvexPolygonsFound=0;
  int numConcavePolygonsFound = 0;

  #ifdef COLLECT_STATISTICS
  if(arePolygonsConnected) {
    #pragma omp atomic
    statistics.numConnectedPolygonSubdivisionOfTriangles++;
  } else {
    #pragma omp atomic
    statistics.numDisconnectedPolygonSubdivisionOfTriangles++;
  }

  #endif

  const int sizePolygonsVector = polygons.size();

  int numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes = newPolygonsGeneratedFromRetesselation.size();
 
  if(arePolygonsConnected) {
    //since the graph is connected, all polygons in "polygons" are oriented similarly
    //we need to check if this orientation is equal to the orientation of the original triangle
    //if it is not, all the polygons will have an orientation that is contrary to the triangle's orientation

    //Let's copy the polygons (except the exterior polygon) to the newPolygonsGeneratedFromRetesselation vector 
    //once we find the exterior polygon, we check if the orientation is correct

    //cerr << "Number of vertices in polygons in triangle (including exterior): " << sizePolygonsVector << endl;
    //cerr << "Triangle: " << t.p[0] << " " << t.p[1] << " " << t.p[2] << endl;
    //for(int v:polygons) {
    //  Point &a = *getPointFromVertexId(v,meshWhereTriangleIs);
    //  cerr << v << " ( " << a[0].get_d() << " , " << a[1].get_d() << " ) \n";
    //}
    //cout << endl;
    

    bool isOrientationOfRetesselatedEqualToTriangle = true;
    //newPolygonsGeneratedFromRetesselation.reserve(numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes+sizePolygonsVector-1);

    bool exteriorPolygonAlreadyFound = false; //have we already found the exterior polygon?

    //for each polygon p inside the triangle, we will triangulate p
    for(int firstElementPolygon=0;firstElementPolygon<sizePolygonsVector;) { 
      int lastElementPolygon = firstElementPolygon+1;
      while(polygons[lastElementPolygon]!= polygons[firstElementPolygon]) lastElementPolygon++;

      //let's process the polygon defined by the vertices in polygons[firstElementPolygon..lastElementPolygon]

      //First, let's check if this polygon is the exterior polygon (since all polygons are connected,
      //it will be iff the three vertices of the original triangle are in the polygon)
      if(!exteriorPolygonAlreadyFound) {
        int numVerticesOriginalTriangleInThisPolygon = 0;
        for(int i=firstElementPolygon;i<lastElementPolygon;i++) {
          if(polygons[i]==t.getInputVertex(0) || polygons[i]==t.getInputVertex(1) || polygons[i]==t.getInputVertex(2)) numVerticesOriginalTriangleInThisPolygon++;
        }
        if(numVerticesOriginalTriangleInThisPolygon>=3) {          
          //this will be 0 if false, 1 if has all edges and in same orientation, -1 if has all edges but in reverse orientation..
          int polygonContainAllEdges = doesPolygonContainAllEdges(polygons,firstElementPolygon,lastElementPolygon,edgesFromTriangleBoundary);
          if(polygonContainAllEdges!=0) { //we found the boundary!
            exteriorPolygonAlreadyFound = true;
            isOrientationOfRetesselatedEqualToTriangle = polygonContainAllEdges==-1; //if the exterior triangle is oriented in the countrary way --> the polygons will be oriented similarly to the triangle
            firstElementPolygon = lastElementPolygon+1;
            continue; //we will not add the exterior polygon to the list of polygons in retesselation...
          }     
        }
      }

      //let's copy the polygon to the list of polygons from retesselation...
      newPolygonsGeneratedFromRetesselation.push_back(BoundaryPolygon(whatPlaneProjectTriangleTo));
      BoundaryPolygon &polygon = newPolygonsGeneratedFromRetesselation.back();
      polygon.vertexSequence.assign(polygons.begin()+firstElementPolygon,polygons.begin()+lastElementPolygon+1);
      polygon.whatPlaneProjectTriangleTo = whatPlaneProjectTriangleTo;

      firstElementPolygon = lastElementPolygon+1;
    }

    int numberPolygonsFromRetesselationNow = newPolygonsGeneratedFromRetesselation.size();
    //let's assign the polygon orientation...
    assert(exteriorPolygonAlreadyFound);
    if(isOrientationOfRetesselatedEqualToTriangle)
      for(int i=numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes;i<numberPolygonsFromRetesselationNow;i++) { //the polygons we are adding are only in this interval (the vector may have other polygons before we added new ones...)
        newPolygonsGeneratedFromRetesselation[i].above = t.above;
        newPolygonsGeneratedFromRetesselation[i].below = t.below;
      }
    else 
      for(int i=numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes;i<numberPolygonsFromRetesselationNow;i++) {
        newPolygonsGeneratedFromRetesselation[i].above = t.below;
        newPolygonsGeneratedFromRetesselation[i].below = t.above;
      }

    //if the cross product of the vectors: (vertex1-vertex0) and (vertex2-vertex0) is negative --> the triangle is oriented
    //in clockwise orientation --> the seed is also in clockwise orientation...
    bool isBoundaryOriginalTriangleClockwisedOriented = meshIntersectionGeometry.isTriangleClockwisedOriented(t,whatPlaneProjectTriangleTo,tempVars.tempVarsIsTriangleClockwisedOriented);

    bool areInteriorPolygonsClockwiseOriented = true;
    if(isBoundaryOriginalTriangleClockwisedOriented && !isOrientationOfRetesselatedEqualToTriangle) areInteriorPolygonsClockwiseOriented = false;
    if(!isBoundaryOriginalTriangleClockwisedOriented && isOrientationOfRetesselatedEqualToTriangle) areInteriorPolygonsClockwiseOriented = false;

   // cerr << "Polygons before and after: " << numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes << " " << newPolygonsGeneratedFromRetesselation.size() << endl;
    //cerr << "Orienting polygon.." << endl;
    //Now, we will reorient the polygons such that all polygons will be oriented clockwise
    //orientPolygonsClockwise(newPolygonsGeneratedFromRetesselation.begin()+numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes, 
   //                        newPolygonsGeneratedFromRetesselation.end(),seedEdge,isSeedEdgeClockwisedOriented,
    //                        isOrientationOfRetesselatedEqualToTriangle,meshWhereTriangleIs);

    if(!areInteriorPolygonsClockwiseOriented) {
      //we need to reverse the orientation of all polygons...
      for(int i=numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes;i<numberPolygonsFromRetesselationNow;i++) 
        newPolygonsGeneratedFromRetesselation[i].reverseVerticesOrder();
    }

    #ifdef SANITY_CHECKS
    
    for(int i=numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes;i<numberPolygonsFromRetesselationNow;i++) {
      bool isOrientedCorrectly = newPolygonsGeneratedFromRetesselation[i].isInClockwiseDirection(vertices, meshWhereTriangleIs);
      if(!isOrientedCorrectly) {
        #pragma omp critical
          {
          cerr << "Error... polygon " << i << " generated from triangle "  << " is not oriented correctly..\n";
          cerr << numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes << " " << numberPolygonsFromRetesselationNow << endl;
          cerr << "isBoundaryOrigintalTriOrientedClockwise " << isBoundaryOriginalTriangleClockwisedOriented << endl;
          cerr << "isOrientationOfRetesselatedEqualToTriangle " << isOrientationOfRetesselatedEqualToTriangle << endl;
          cerr << "Original triangle: " << endl;
          cerr << t.p[0] << " " << t.p[1] << " " << t.p[2] << endl;
          printVertexForDebugging(getPointFromVertexId(t.p[0],meshWhereTriangleIs)->data());
          printVertexForDebugging(getPointFromVertexId(t.p[1],meshWhereTriangleIs)->data());
          printVertexForDebugging(getPointFromVertexId(t.p[2],meshWhereTriangleIs)->data());
          cerr << "Polygon: " << endl;
          for(int v:newPolygonsGeneratedFromRetesselation[i].vertexSequence) cerr << v << " "; cerr << endl;
          for(int v:newPolygonsGeneratedFromRetesselation[i].vertexSequence) printVertexForDebugging(getPointFromVertexId(v,meshWhereTriangleIs)->data());
          cerr << "What plane: " << whatPlaneProjectTriangleTo << endl;
          cerr << "Edges from planar graph " << endl;
          for(auto a:edgesUsedInThisTriangleV) {
            cerr << a.first << " " << a.second << endl;
            printVertexForDebugging(getPointFromVertexId(a.first,meshWhereTriangleIs)->data());
            printVertexForDebugging(getPointFromVertexId(a.second,meshWhereTriangleIs)->data());
            cerr << endl;
          }
          cerr << "Polygons extracted " << endl;
          for(int j=numberPolygonsFromRetesselationInVectorBeforeWeAddedNewOnes;j<numberPolygonsFromRetesselationNow;j++) {
            for(int v:newPolygonsGeneratedFromRetesselation[j].vertexSequence) cerr << v << " "; cerr << endl;
            for(int v:newPolygonsGeneratedFromRetesselation[j].vertexSequence) printVertexForDebugging(getPointFromVertexId(v,meshWhereTriangleIs)->data());
          }
        }
      }
      assert(isOrientedCorrectly);
    }
    #endif
   // cerr << "Oriented \n";
  } /*  else {
    //TODO
    //TODO: remember to set whatPlaneProjectTriangleTo of the polygon...
    //let's print the triangles with disconnected polygons in its interior for debugging purposes...
    //edgesUsedInThisTriangle
    cerr << "Disconnected triangle:::: this situation is not treated yet..." << endl << endl << endl;
    int numDisconnectedTrianglesNow = 0;
    numDisconnectedTrianglesNow = statistics.numDisconnectedPolygonSubdivisionOfTriangles;
    stringstream trianglesDisconnectedPath;
    trianglesDisconnectedPath << "triangleDisconnectedPath_" << numDisconnectedTrianglesNow << ".gts";
    string path = trianglesDisconnectedPath.str();
    saveEdgesAsGTS(edgesUsedInThisTriangle, meshWhereTriangleIs,  path);
  }
*/


  #ifdef COLLECT_STATISTICS
    #pragma omp atomic
    statistics.ctConvexPolygonsInTriangleRetesselations += numConvexPolygonsFound;
    #pragma omp atomic
    statistics.ctConcavePolygonsInTriangleRetesselations += numConcavePolygonsFound;
  #endif

  
}













//edges represent the edges generated from the intersection of the intersecting triangles.
//intersectingTrianglesThatGeneratedEdges[i] contains the pair of triangles whose intersection generated edges[i]
void retesselateIntersectingTriangles(MeshIntersectionGeometry & meshIntersectionGeometry, 
                                      const vector< pair<VertexFromIntersection, VertexFromIntersection> >  &edgesFromIntersection, 
                                      const vector< pair<InputTriangle *,InputTriangle *> > &intersectingTrianglesThatGeneratedEdges,
                                      vector<BoundaryPolygon> polygonsFromRetesselation[2]) {

  assert(edgesFromIntersection.size()==intersectingTrianglesThatGeneratedEdges.size());
  timespec t0,t1;
  cerr << "Retesselating triangles..." << "\n";
  clock_gettime(CLOCK_REALTIME, &t0);


  StatisticsAboutRetesseation statisticsAboutRetesseation;
  

  unordered_map<const InputTriangle *, vector<int> > intersectingEdgesInEachTriangle[2];
  //trianglesIntersectingEachTriangle[0] --> for each triangle t from map 0 (that intersects other triangles), 
  //trianglesIntersectingEachTriangle[0][t] is a vector of triangles from mesh 1 intersecting t...

  const int numEdges = edgesFromIntersection.size();
  for(int i=0;i<numEdges;i++) { 
    //the intersection of tA (from mesh 0) with tB generated the i-th edge...
    const InputTriangle *tA = intersectingTrianglesThatGeneratedEdges[i].first;
    const InputTriangle *tB = intersectingTrianglesThatGeneratedEdges[i].second;

    intersectingEdgesInEachTriangle[0][tA].push_back(i);
    intersectingEdgesInEachTriangle[1][tB].push_back(i);
  }
 
  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to extract the edges intersecting each triangle: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 


  //for each triangle t from mesh meshId, intersectingEdgesInEachTriangle[meshId][t] contains a vector of ids of edges formed by the intersection
  //of t with other triangles...


  //TODO: special cases, collinear, etc..

  vector<pair<int,int> > newTriEdgesFromEachMap[2];

  for(int meshIdToProcess=0;meshIdToProcess<2;meshIdToProcess++) {
    vector< const pair<const InputTriangle * const, vector<int> >* > intersectingEdgesInEachTriangleToProcess;
    for(const auto &elem:intersectingEdgesInEachTriangle[meshIdToProcess]) intersectingEdgesInEachTriangleToProcess.push_back(&elem);

    //for each triangle in this mesh, we have a vector of the edges in this triangle that need to be processed during the retesselation...
    //the vector of the edges actually have the ids (position in the edgesFromIntersection vector)

    int numTrianglesToProcess = intersectingEdgesInEachTriangleToProcess.size();

    cerr << "Number of triangles to retesselate in this mesh: " << numTrianglesToProcess << "\n";

    statisticsAboutRetesseation.ctTrianglesRetesselate += numTrianglesToProcess;

    #pragma omp parallel 
    {
      //vector<pair<int,int> > myNewTriEdgesFromEachMap;

      vector<BoundaryPolygon> myNewPolygonsFromRetesselation;
      //myNewPolygonsFromRetesselation.reserve(numTrianglesToProcess);

      TempVarsRetesselateTriangleFunction tempVars;

      const int percentShowLog = 1;
      int onePercentNumTriangles = (percentShowLog*numTrianglesToProcess)/100;
      if(onePercentNumTriangles==0) onePercentNumTriangles = 1;

      #pragma omp for
      for(int i=0;i<numTrianglesToProcess;i++) {      
        if((i%onePercentNumTriangles)==0) {
          #pragma omp critical
          clog << "Retesselating " << i << " of " << numTrianglesToProcess << " Percent= " << (i*100)/numTrianglesToProcess << "\n";      
        }
        //for each triangle ts, let's process the triangles intersecting t and the edges
        //formed by the intersection of ts and these triangles

        const auto &ts  = *intersectingEdgesInEachTriangleToProcess[i];
        //ts.first = a triangle from map 0
        //ts.second = list of edges formed by the intersection of ts.first with other triangles...
        const auto &t = *ts.first;
        const auto &edgesFromIntersectionThisTriangle =  ts.second;

        
        #ifdef COLLECT_STATISTICS_PRINT_TRIANGLES_INTERSECTIONS
          #pragma omp critical
          //cout << "Triangle from mesh " << meshIdToProcess << " intersects: " << ts.second.size() << "\n";
          if(ts.second.size()==6) {
            #pragma omp critical
            printTriangleWithManyIntersection(t,meshIdToProcess,ts.second,intersectingTrianglesThatGeneratedEdges);
          }
          if(ts.second.size()==6) {
            #pragma omp critical
            saveEdgesAsGTS(edges,ts.second,"edgesFromTriangleIntersetingManyTriangles6.gts");
          }
         if(meshIdToProcess==0) {
            #pragma omp critical
            saveEdgesAsGTS(edges,ts.second,"edgesFromTriangleIntersetingMesh0.gts");
          }
        #endif

        //retesselateTriangle(edges, edgesFromIntersection,t, meshIdToProcess, myNewTrianglesFromRetesselation, tempVars, statisticsAboutRetesseation);
        retesselateTriangleUsingWedgeSorting(meshIntersectionGeometry,edgesFromIntersection,edgesFromIntersectionThisTriangle,t, myNewPolygonsFromRetesselation, tempVars, statisticsAboutRetesseation);

      }

      #pragma omp critical 
      {
        //cerr << "Number of polygons in vector of retesselated: " << myNewPolygonsFromRetesselation.size() << endl;
        //cerr << "Capacity: " << myNewPolygonsFromRetesselation.capacity() << endl;
        //newTriEdgesFromEachMap[meshIdToProcess].insert( newTriEdgesFromEachMap[meshIdToProcess].end(),myNewTriEdgesFromEachMap.begin(),myNewTriEdgesFromEachMap.end());
        //trianglesFromRetesselation[meshIdToProcess].insert(trianglesFromRetesselation[meshIdToProcess].end(),myNewTrianglesFromRetesselation.begin(),myNewTrianglesFromRetesselation.end());
        polygonsFromRetesselation[meshIdToProcess].insert(polygonsFromRetesselation[meshIdToProcess].end(),myNewPolygonsFromRetesselation.begin(),myNewPolygonsFromRetesselation.end() );
      }

      //break;
    }

    
    for(BoundaryPolygon p:polygonsFromRetesselation[meshIdToProcess]) {
      //assert(p.isInClockwiseDirection(vertices,meshIdToProcess));
    }
  }


  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to retesselate creating new polygons: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n";  

  clock_gettime(CLOCK_REALTIME, &t0);
  cerr << "Triangulating polygons from retesselation...\n";
  for(int meshIdToProcess=0;meshIdToProcess<2;meshIdToProcess++) {
    int numPolygons = polygonsFromRetesselation[meshIdToProcess].size();

    const int percentShowLog = 10;
    int onePercentNumPolygons = (numPolygons*percentShowLog)/100;
    if(onePercentNumPolygons==0) onePercentNumPolygons = 1;

    #pragma omp parallel
    {
      VertCoord tempCoords[2];
      #pragma omp for
      for(int i=0;i<numPolygons;i++) {
        //cerr << "Triangulating " << i << " of " << numPolygons << " Percent= " << i*100/numPolygons << endl;
        if((i%onePercentNumPolygons)==0) {
          clog << "Triangulating " << i << " of " << numPolygons << " Percent= " << (i*100)/numPolygons << "\n";
          clog << "TODO: finish implementation of triangulate polygons...\n";
        }
        //polygonsFromRetesselation[meshIdToProcess][i].triangulatePolygon(meshIntersectionGeometry,meshIdToProcess,tempCoords);
      }
    }
  }

  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to triangulate polygons: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 

  
  cerr << "Counts computed during retesselation: " << "\n";
  statisticsAboutRetesseation.printStatsSoFar();
  cerr << "Number of inters. tests (for insertin tris.)that are true    : " << ctEdgeIntersect << "\n";
  cerr << "Number of inters. tests (for insertin tris.)that are false    : " << ctEdgeDoNotIntersect << "\n";



 
}