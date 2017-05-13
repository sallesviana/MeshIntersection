#include "triangleRetesselation.h"


bool lessThan(const Vertex* const a, const Vertex* const b) {
  //cerr << "Testing pointers for less than..." << endl;
  //cerr << "Answer: " << endl;
  //cerr << a->compare(*b) << endl;
  return a->compare(*b)<0;
}

bool equalTo(const Vertex* const a, const Vertex* const b) {
  return a->compare(*b)==0;
}

bool equalTo(const InputTriangle* const a, const InputTriangle* const b) {
  return a->compare(*b)==0;
}

struct VertexPtrLessComparator {
    bool operator() (const Vertex* const a, const Vertex* const b) const{
        return lessThan(a,b);
    }
};

struct VertexPairPtrLessComparator {
    bool operator() (const pair<const Vertex*,const Vertex*> &a, const pair<const Vertex*,const Vertex*> &b) const{
      if(!equalTo(a.first,b.first)) return lessThan(a.first,b.first);
      return lessThan(a.second,b.second);
    }
};

struct VertexPtrEqualComparator {
    bool operator() (const Vertex* const a, const Vertex* const b) const{
        return equalTo(a,b);
    }
};



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
  //MeshIntersectionGeometry::TempVarsIsBoundaryClockwise tempVarsIsBoundaryClockwise;
  MeshIntersectionGeometry::TempVarsIsTriangleClockwisedOriented tempVarsIsTriangleClockwisedOriented;
  MeshIntersectionGeometry::TempVarsSortEdgesByAngle tempVarsSortEdgesByAngle;
  MeshIntersectionGeometry::TempVarsDoIntersect tempVarsDoIntersect;
};


#define ERROR_CODE -2123123123
int binarySearchFirstOccurrency(const vector<pair<const Vertex *,const Vertex *> > &v,const Vertex * key) {
  int mid;
  int lo = 0;
  int hi = v.size()-1;
  while (lo < hi) {
      mid = lo + (hi-lo)/2;   // note: division truncates
      if (!lessThan(v[mid].first,key))
         hi = mid;
      else
         lo = mid+1;
  }
          
  if (!equalTo(v[lo].first,key))
      return ERROR_CODE;                // p(x) is true for all x in S!
      
  return lo;         // lo is the least x for which p(x) is true
}

bool checkIfPolygonsAreConnected(const vector<pair<const Vertex *,const Vertex *> > &raggedArraySortedEdges) {
  queue<const Vertex *> verticesToProcess;
  set<const Vertex *,VertexPtrLessComparator> processedVertices;
  verticesToProcess.push(raggedArraySortedEdges[0].first);
  processedVertices.insert(raggedArraySortedEdges[0].first);

  
  //visit all vertices reachable from v...
  while(!verticesToProcess.empty()) {
    const Vertex * v = verticesToProcess.front();
    verticesToProcess.pop();

    int posStartVRaggedArray =  binarySearchFirstOccurrency(raggedArraySortedEdges,v);
    //cerr << "v and pos v " << v << " " << posStartVRaggedArray << endl;
    for(int i=posStartVRaggedArray;i<raggedArraySortedEdges.size()&& equalTo(raggedArraySortedEdges[i].first,v);i++) {
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
 

  return true;
}



//returns the index i of the vector v, where v[i][0] = key.first and v[i][1] = key.second 
int binarySearch(vector<array<const Vertex *,3> > &v,const pair<const Vertex *,const Vertex *> &key) {
  int lo = 0, hi = v.size()-1, mid;
  while (lo<=hi) {
    mid = lo + (hi-lo)/2;
    if( equalTo(v[mid][0],key.first) && equalTo(v[mid][1],key.second)) 
      return mid;
    else if ( lessThan(v[mid][0] , key.first) || (equalTo(v[mid][0] , key.first) && lessThan(v[mid][1] , key.second)) ) //is the element v[mid] smaller than the key?
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
void sortEdgesAndExtractPolygonsFromEdgeListUsingWedges(const MeshIntersectionGeometry &meshIntersectionGeometry, const vector<pair<const Vertex *,const Vertex *> > &edges, const int planeProjectTriangleTo,  TempVarsRetesselateTriangleFunction &tempVars ) {
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
  
  
  //edges sharing the first vertex will be together after this sorting...
  sort(raggedArraySortedEdges.begin(),raggedArraySortedEdges.end(), [&](const pair<const Vertex *,const Vertex *> &e1, const pair<const Vertex *,const Vertex *> &e2) {
            if(!equalTo(e1.first,e2.first)) return lessThan(e1.first,e2.first); //sort first by the first vertex..
            else return false;
            //return meshIntersectionGeometry.isAngleWith0Greater(*e1.first, *e1.second, *e2.second, planeProjectTriangleTo, tempVars.tempVarsIsAngleWith0Greater);
        } );
  

  //now, we have to sort edges sharing the first vertex by angle.
  vector<pair<const Vertex *,const Vertex *> >::iterator startSortingByAngle = raggedArraySortedEdges.begin();
  while(startSortingByAngle!=raggedArraySortedEdges.end()) {
    //cerr << "getting end.." << endl;
    vector<pair<const Vertex *,const Vertex *> >::iterator endSortingByAngle = startSortingByAngle;
    while(endSortingByAngle != raggedArraySortedEdges.end() && equalTo(startSortingByAngle->first,endSortingByAngle->first))
      endSortingByAngle++;

    /*cerr << "Edges by before sorting...: " << endl;
    for(auto i=startSortingByAngle;i<endSortingByAngle;i++) {
      auto a = *i;
      a.first->print(); cerr << endl;
      a.second->print(); cerr << endl;
      cerr << endl; 
    }*/

    meshIntersectionGeometry.sortEdgesSharingStartingVertexByAngle(startSortingByAngle,endSortingByAngle,planeProjectTriangleTo,tempVars.tempVarsSortEdgesByAngle);
    
    /*cerr << "Edges by after sorting...: " << endl;
    for(auto i=startSortingByAngle;i<endSortingByAngle;i++) {
      auto a = *i;
      a.first->print();  cerr << endl;
      a.second->print();  cerr << endl;
      cerr << endl; 
    }*/


    startSortingByAngle = endSortingByAngle;
  }

  /*cerr << "Edges by angle...: " << endl;
  for(auto a:raggedArraySortedEdges) {
    a.first->print(); cerr << endl;
    a.second->print(); cerr << endl;
    cerr << endl; 
  }
  
  cerr << "Extracting wedges..." << endl;*/
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

      //cerr << "Equal to? " << nextElement << " " << currElemnt << " " << equalTo(raggedArraySortedEdges[nextElement].first , raggedArraySortedEdges[currElemnt].first) << endl;
      //did we reach the end of a group?
      if(nextElement >= numDirectedEdges || !equalTo(raggedArraySortedEdges[nextElement].first , raggedArraySortedEdges[currElemnt].first)) {
        //if yes, we will create a wedge from currElement to the first one! (to complete the wedges from this group..)
        //the next element is actually the first element in this group...
        wedgeToAdd[0] = raggedArraySortedEdges[firstElementGroup].second; //c
        wedgeToAdd[1] = raggedArraySortedEdges[firstElementGroup].first; //a
        wedgeToAdd[2] = raggedArraySortedEdges[currElemnt].second;  //b

        /*cerr << "END: first group, curr, next: " << firstElementGroup << " " << currElemnt << " " << nextElement << endl;
        wedgeToAdd[0]->print(); cerr << endl;
        wedgeToAdd[1]->print(); cerr << endl;
        wedgeToAdd[2]->print(); cerr << endl;
        cerr << endl;*/

        firstElementGroup = currElemnt+1; //let's process the next group!
        break;
      } else {
        //else, we will create a wedge from the current to the next...
        wedgeToAdd[0] = raggedArraySortedEdges[nextElement].second; //c
        wedgeToAdd[1] = raggedArraySortedEdges[nextElement].first; //a
        wedgeToAdd[2] = raggedArraySortedEdges[currElemnt].second;  //b

        /*cerr << "first group, curr, next: " << firstElementGroup << " " << currElemnt << " " << nextElement << endl;
        wedgeToAdd[0]->print(); cerr << endl;
        wedgeToAdd[1]->print(); cerr << endl;
        wedgeToAdd[2]->print(); cerr << endl;
        cerr << endl;*/
      }
    }
  }

  /*cerr << "Wedges before sorting...: " << endl;
  for(auto a:wedgesTemp) {
    a[0]->print(); cerr << endl;
    a[1]->print(); cerr << endl;
    a[2]->print(); cerr << endl;
    cerr << endl; 
  }*/

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
  sort(wedgesTemp.begin(),wedgesTemp.end(), [](const array<const Vertex *,3> &a, const array<const Vertex *,3> &b) {
                                              if(!equalTo(a[0],b[0])) return lessThan(a[0],b[0]);
                                              if(!equalTo(a[1],b[1])) return lessThan(a[1],b[1]);
                                              return lessThan(a[2],b[2]);
                                            });
  const int numWedges = wedgesTemp.size();
  usedWedgesTemp.resize(numWedges);
  for(int i=0;i<numWedges;i++) usedWedgesTemp[i] = false;

  //Let's use a ragged array for performance...
  polygons.resize(0);  

  /*cerr << "Wedges: " << endl;
  for(auto a:wedgesTemp) {
    a[0]->print(); cerr << endl;
    a[1]->print(); cerr << endl;
    a[2]->print(); cerr << endl;
    cerr << endl; 
  }*/


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

      if(equalTo(wedgesTemp[next][1] , wedgesTemp[wedgeStart][0])) { //have I completed a polygon?
        polygons.push_back(wedgesTemp[wedgeStart][0]);
        break; //let's find a (possible) next vertex to start a new polygon...
      }
    }
  }

  
}








//returns true iff the candidate edge (represented using the ids of the vertices in meshIdToProcess -- negative vertices are in the "common layer")
//intersects an edge from the set edgesToTest (these edges are in the same layer).
//we do not consider intersections in endpoints..
//TODO: consider special cases: vertical triangle, parallel edges, etc...
bool intersects(const MeshIntersectionGeometry &meshIntersectionGeometry,const pair<const Vertex *,const Vertex *> &candidateEdge,
                  const set<pair<const Vertex *,const Vertex *>,VertexPairPtrLessComparator > &edgesToTest,
                  int meshIdToProcess, int whatPlaneProjectTriangleTo, TempVarsRetesselateTriangleFunction &tempVars) {
  //cerr << "Testing intersection: " << candidateEdge.first << " " << candidateEdge.second << endl;

  for(const pair<const Vertex *,const Vertex *> &edgeToTest:edgesToTest) {
    bool inter = meshIntersectionGeometry.doIntersect(candidateEdge,edgeToTest,whatPlaneProjectTriangleTo,tempVars.tempVarsDoIntersect);
    
    if(inter) 
      return true;
  }
  return false;
}


void addEdgesToTriangleAndMakeItConnected(const MeshIntersectionGeometry &meshIntersectionGeometry, 
                                          const vector<const Vertex *> vertices, 
                                          vector<pair<const Vertex *,const Vertex *> > &edges, 
                                          const int whatPlaneProjectTriangleTo, const int meshWhereTriangleIs, 
                                          TempVarsRetesselateTriangleFunction &tempVars ) {
  //TODO: count statistics ...



  //this set will store the edges we used in the retesselated triangle
  //we need to choose what edges to create and, then, use these new edges to reconstruct the triangulation..
  //This set will have the edges that will form the retriangulation of ts..
  set<pair<const Vertex *,const Vertex *>,VertexPairPtrLessComparator > setEdgesUsedInTriangle(edges.begin(),edges.end());


  //cerr << "Edges before adding extra edges using brute force: " << setEdgesUsedInTriangle.size() << endl;

  const int numVTriangle =  vertices.size();
  for(int i=0;i<numVTriangle;i++)
    for(int j=i+1;j<numVTriangle;j++) {
      const Vertex * v1 = vertices[i];
      const Vertex * v2 = vertices[j];
            
      //let's try to insert edge (v1,v2)...
      pair<const Vertex *,const Vertex *> candidateEdge(v1,v2);
      if(setEdgesUsedInTriangle.count(candidateEdge)!=0) continue; //the edge was already used...
      if(setEdgesUsedInTriangle.count( pair<const Vertex *,const Vertex *>(v2,v1))!=0) continue; //the edge was already used...


      //if edge e=(v1,v2) does not intersect other edges already inserted in this triangle,
      //we will add e to the new triangulation...
      //any triangulation is fine as soon as the edges from the intersection of the meshes are there...
      if(!intersects(meshIntersectionGeometry,candidateEdge,setEdgesUsedInTriangle,meshWhereTriangleIs,whatPlaneProjectTriangleTo,tempVars)) {
        setEdgesUsedInTriangle.insert(candidateEdge); //if it does not intersect, we can safely add this edge to the triangulation...
      } 
    }       

  //cerr << "Edges after adding more: " << setEdgesUsedInTriangle.size() << endl << endl;
    
  edges = vector<pair<const Vertex *,const Vertex *> >(setEdgesUsedInTriangle.begin(),setEdgesUsedInTriangle.end());
}





int doesPolygonContainAllEdges(const vector<const Vertex *> &polygons,const int firstElementPolygon,const int lastElementPolygon,
                                const set<pair<const Vertex *,const Vertex *>, VertexPairPtrLessComparator > & edgesFromTriangleBoundary) {

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

  //cerr << "Adding edges from intersection" << endl;
  /*for(auto i:edgesFromIntersectionThisTriangle) {
    cerr << edgesFromIntersection[i].first.getMeshId() << " " << edgesFromIntersection[i].first.getId() << endl;
    cerr << edgesFromIntersection[i].second.getMeshId() << " " << edgesFromIntersection[i].second.getId() << endl << endl;
  }*/
  //cerr << "End print..." << endl;


  //Now, we need to add more edges to fill the parts of ts that still do not form triangle
  vector<const Vertex *> verticesToTesselate;
  //what vertices will we have in the new triangulation of ts? (the original vertices + the vertices of the edges from intersections)
  for(const auto &p:edgesUsedInThisTriangle) {
    verticesToTesselate.push_back(p.first);
    verticesToTesselate.push_back(p.second);
  }

  //cerr << "Sorting..." << endl;
  sort(verticesToTesselate.begin(),verticesToTesselate.end(),[](const Vertex *a, const Vertex *b){
                                return lessThan(a,b); // (*a)<(*b);
                              });
  //cerr << "Unique..." << endl;
  auto newEnd = unique(verticesToTesselate.begin(),verticesToTesselate.end(),[](const Vertex *a, const Vertex *b){
                                return equalTo(a,b); // (*a)==(*b)
                              });
  verticesToTesselate.resize(newEnd-verticesToTesselate.begin());

  //Initially we will insert only the vertices from the intersection (the original vertices of the triangle will be inserted later)

    
  const int numVerticesToTesselateNow = verticesToTesselate.size();
  const InputVertex* tv0 = t.getInputVertex(0);
  const InputVertex* tv1 = t.getInputVertex(1);
  const InputVertex* tv2 = t.getInputVertex(2);
  const int meshContainingT = t.getMeshId();

  //Computes the vertices incident to each edge of the original triangle...
  //verticesIncidentEachEdgeOriginalTriangle[0] is for edge t.p[0]-t.p[1]
  //verticesIncidentEachEdgeOriginalTriangle[1] is for edge t.p[1]-t.p[2]
  //verticesIncidentEachEdgeOriginalTriangle[2] is for edge t.p[2]-t.p[0]
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

    //if( (v.edge[0] == (*tv0) && v.edge[1]  == (*tv1)) || (v.edge[0]  == (*tv1) && v.edge[1]  == (*tv0)) ) {
    if(  (equalTo(&(v.edge[0]),tv0) && equalTo(&(v.edge[1]),tv1)) || (equalTo(&(v.edge[0]),tv1) && equalTo(&(v.edge[1]),tv0)) ) {
      verticesIncidentEachEdgeOriginalTriangle[0].push_back(&v);
      continue;
    }   
    //if( (v.edge[0]  == (*tv1) && v.edge[1]  == (*tv2)) || (v.edge[0]  == (*tv2) && v.edge[1]  == (*tv1)) ) {
    if(  (equalTo(&(v.edge[0]),tv1) && equalTo(&(v.edge[1]),tv2)) || (equalTo(&(v.edge[0]),tv2) && equalTo(&(v.edge[1]),tv1)) ) {
      verticesIncidentEachEdgeOriginalTriangle[1].push_back(&v);
      continue;
    }
    //if( (v.edge[0]  == (*tv2) && v.edge[1]  == (*tv0)) || (v.edge[0]  == (*tv0) && v.edge[1]  == (*tv2)) ) {
    if(  (equalTo(&(v.edge[0]),tv2) && equalTo(&(v.edge[1]),tv0)) || (equalTo(&(v.edge[0]),tv0) && equalTo(&(v.edge[1]),tv2)) ) {
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


  //cerr << "Adding 3 edges" << endl;
  /*for(auto e:edgesUsedInThisTriangle) {
    cerr << e.first->getMeshId() << " " << e.first->getId() << endl;
    cerr << e.second->getMeshId() << " " << e.second->getId() << endl << endl;
  }*/
  //cerr << "End print..." << endl;
  
  set<pair<const Vertex *,const Vertex *>, VertexPairPtrLessComparator > edgesFromTriangleBoundary;


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

      assert(!equalTo(candidateEdge.first , candidateEdge.second));
      edgesUsedInThisTriangle.push_back(candidateEdge);
      edgesFromTriangleBoundary.insert(candidateEdge);
    } else {
      if(verticesIncidentEdge==1) {
        const Vertex * v1 = t.getVertex(v1Pos);
        const Vertex * v2 = t.getVertex(v2Pos);

        const Vertex * v =  *(verticesIncidentEachEdgeOriginalTriangle[edge].begin());  

        /*cerr << "We have 1 vertex in edge... " << endl;
        cerr << v->getMeshId() << " " << v->getId() << endl;
        cerr << v1->getMeshId() << " " << v1->getId() << endl;
        cerr << v2->getMeshId() << " " << v2->getId() << endl << endl;*/

                
        
        pair<const Vertex *,const Vertex *>  candidateEdge;    

        candidateEdge.first = v1; 
        candidateEdge.second = v;


        #ifdef COLLECT_STATISTICS
          #pragma omp atomic
          statistics.ctEdgesActuallyInsertedInRetesselation++;
        #endif

        assert(!equalTo(candidateEdge.first , candidateEdge.second));
        edgesUsedInThisTriangle.push_back(candidateEdge);
        edgesFromTriangleBoundary.insert(candidateEdge);

        candidateEdge.first = v; 
        candidateEdge.second = v2;        

        #ifdef COLLECT_STATISTICS
          #pragma omp atomic
          statistics.ctEdgesActuallyInsertedInRetesselation++;
        #endif
        assert(!equalTo(candidateEdge.first , candidateEdge.second));
        edgesUsedInThisTriangle.push_back(candidateEdge);
        edgesFromTriangleBoundary.insert(candidateEdge);

        /*for(auto e:edgesUsedInThisTriangle) {
          cerr << e.first->getMeshId() << " " << e.first->getId() << endl;
          cerr << e.second->getMeshId() << " " << e.second->getId() << endl << endl;
        }*/
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

  /*for(auto e:edgesUsedInThisTriangle) {
    cerr << e.first->getMeshId() << " " << e.first->getId() << endl;
    cerr << e.second->getMeshId() << " " << e.second->getId() << endl << endl;
  }*/

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
                              //return meshIntersectionGeometry.isCloser(*edgeOrig,*v1,*v2,tempVars.tempVarsIsCloser);
                              //until now all the vertices in the vector are from intersection...
                              //thus, we can static cast...
                              return meshIntersectionGeometry.isCloser(*edgeOrig,*static_cast<const VertexFromIntersection*>(v1),*static_cast<const VertexFromIntersection*>(v2),tempVars.tempVarsIsCloser);
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
    assert(!equalTo(candidateEdge.first , candidateEdge.second));
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
      assert(!equalTo(candidateEdge.first , candidateEdge.second));
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

  /*for(auto e:edgesUsedInThisTriangle) {
    cerr << e.first->getMeshId() << " " << e.first->getId() << endl;
    cerr << e.second->getMeshId() << " " << e.second->getId() << endl << endl;
  }*/

  //-----------------------------------------------------------------------------------------------------------------------
  //Now, it is finally time to start creating the polygons, checking if they are connected, etc...

  //TODO: why do we need a copy of edgesUsedInThisTriangle?
  vector<pair<const Vertex *,const Vertex *> > edgesUsedInThisTriangleV(edgesUsedInThisTriangle.begin(),edgesUsedInThisTriangle.end());

  const int whatPlaneProjectTriangleTo = meshIntersectionGeometry.getPlaneTriangleIsNotPerpendicular(t,tempVars.tempVarsGetPlaneTriangleIsNotPerpendicular);

  sortEdgesAndExtractPolygonsFromEdgeListUsingWedges(meshIntersectionGeometry, edgesUsedInThisTriangleV, whatPlaneProjectTriangleTo,  tempVars );


  //at the end of the previous function raggedArrayEdges will contain the directed edges sorted by their first vertex...
  bool arePolygonsConnected = checkIfPolygonsAreConnected(tempVars.raggedArraySortedEdges);

  if(!arePolygonsConnected) {
    static int numDisconnectedTrianglesNow = 0;
    
    /*
    {
      numDisconnectedTrianglesNow++;
      stringstream trianglesDisconnectedPath;
      trianglesDisconnectedPath << "triangleDisconnectedPath_before_" << numDisconnectedTrianglesNow << ".gts";
      string path = trianglesDisconnectedPath.str();
      meshIntersectionGeometry.saveEdgesAsGTS(edgesUsedInThisTriangleV,  path);
    }*/


    //the graph is not connected... let's add edges using "brute force" and make it connected...
    //since this should happen rarely, we can use a simple algorithm to connect the graph...


    addEdgesToTriangleAndMakeItConnected(meshIntersectionGeometry,verticesToTesselate, edgesUsedInThisTriangleV,whatPlaneProjectTriangleTo,meshContainingT,tempVars);

    //now the graph is connected... let's run the polygon extraction algorithm again 
    sortEdgesAndExtractPolygonsFromEdgeListUsingWedges(meshIntersectionGeometry, edgesUsedInThisTriangleV,whatPlaneProjectTriangleTo,  tempVars );

 
    //TODO: remove...
    //now the graph have to be connected!
    arePolygonsConnected = checkIfPolygonsAreConnected(tempVars.raggedArraySortedEdges);
    assert(arePolygonsConnected);
    //TODO: update are polygons connected....


    /*
    {
      stringstream trianglesDisconnectedPath;
      trianglesDisconnectedPath << "triangleDisconnectedPath_after_" << numDisconnectedTrianglesNow << ".gts";
      string path = trianglesDisconnectedPath.str();
      meshIntersectionGeometry.saveEdgesAsGTS(edgesUsedInThisTriangleV,  path);
    }*/
  }

  //these are the polygons extracted by the polygon extraction algorithm...
  vector<const Vertex *> &polygons = tempVars.polygons; 


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
 

  //cerr << "Checking orientation.." << endl;
  //TOOD: remove the "if"
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

    //for(auto i:polygons) cerr << i->getMeshId() << " " << i->getId()  << " , "; cerr << endl;

    //for each polygon p inside the triangle, we will triangulate p
    for(int firstElementPolygon=0;firstElementPolygon<sizePolygonsVector;) { 
      int lastElementPolygon = firstElementPolygon+1;
      //cerr << "Bef. while" << endl;
      while(!equalTo(polygons[lastElementPolygon], polygons[firstElementPolygon])) lastElementPolygon++;
     // cerr << "After while\n\n";
      //let's process the polygon defined by the vertices in polygons[firstElementPolygon..lastElementPolygon]

      //First, let's check if this polygon is the exterior polygon (since all polygons are connected,
      //it will be iff the three vertices of the original triangle are in the polygon)
      if(!exteriorPolygonAlreadyFound) {
        int numVerticesOriginalTriangleInThisPolygon = 0;
        for(int i=firstElementPolygon;i<lastElementPolygon;i++) {
          if( equalTo(polygons[i],t.getInputVertex(0)) || equalTo(polygons[i],t.getInputVertex(1)) || equalTo(polygons[i],t.getInputVertex(2))) numVerticesOriginalTriangleInThisPolygon++;
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

 // cerr << "End" << endl;
  #ifdef COLLECT_STATISTICS
    #pragma omp atomic
    statistics.ctConvexPolygonsInTriangleRetesselations += numConvexPolygonsFound;
    #pragma omp atomic
    statistics.ctConcavePolygonsInTriangleRetesselations += numConcavePolygonsFound;
  #endif

  
}





//each boundary polygon stores the information about the polygons adjacent to each edge and the
//polyhedra on each side of its edges from intersection. Let's update this information now!
void updateBoundaryPolygonWithEdgeAdjacencyInformation( vector<BoundaryPolygon> &boundaryPolygonsFromThisTriangle , // the boundary polygons we will update
                                                        const InputTriangle &triangleBeingRetesselated, // the triangle being retesselated...
                                                        const vector< pair<VertexFromIntersection, VertexFromIntersection> >  &edgesFromIntersection, //the edges from the intersection (all... we are interested only in the ones in edgesFromIntersectionThisTriangle)
                                                        const vector< pair<InputTriangle *,InputTriangle *> > &intersectingTrianglesThatGeneratedEdges, //the pairs of triangles that generated each edge from intersection.
                                                        const vector<int> &edgesFromIntersectionThisTriangle) // the positions (in the previous vectors) of the edges from intersection that are in this triangle...
{
  //we need to things in order to process the data:
  //1) Given a pair of vertices pointers, we should be able to determine which boundary polygon contain it...
  //2) Given an edge from the intersection, we should be able to determine which polyhedra from the other mesh bounds the (other) triangle that generated this edge...

  map<pair<const Vertex *,const Vertex * >, BoundaryPolygon *, VertexPairPtrLessComparator> orientedEdgeToBoundaryPolygon;
  for(BoundaryPolygon&bp:boundaryPolygonsFromThisTriangle) {
    const auto &vertexSequence = bp.vertexSequence;
    const int numEdges = vertexSequence.size()-1;
    for(int i=0;i<numEdges;i++) { //edge i connects vertexSequence[i] to vertexSequence[i+1] (in the other polygon we will have [i+1] to [i])
      orientedEdgeToBoundaryPolygon[make_pair(vertexSequence[i],vertexSequence[i+1])] = &bp;
    }
  }

  int meshIdRetesselatedTriangle = triangleBeingRetesselated.getMeshId();

  //this map will not store oriented edges...
  //for each edge from intersection, it will store the objects of the other mesh bounded by that edge...
  map<pair<const Vertex *,const Vertex * >, pair<ObjectId,ObjectId>, VertexPairPtrLessComparator > objectsOtherMeshBoundedByEachEdge;

  for(int edgeId:edgesFromIntersectionThisTriangle) {
    const pair<VertexFromIntersection, VertexFromIntersection> &pairVerticesOfEdge = edgesFromIntersection[edgeId];
    const pair<InputTriangle *,InputTriangle *> &pairTrianglesGeneratedEdge = intersectingTrianglesThatGeneratedEdges[edgeId];

    const InputTriangle* triangleFromOtherMeshGeneratedEdge = (meshIdRetesselatedTriangle==0)?pairTrianglesGeneratedEdge.second:pairTrianglesGeneratedEdge.first;
    const InputTriangle* triangleThisMeshGeneratedEdge = (meshIdRetesselatedTriangle==1)?pairTrianglesGeneratedEdge.second:pairTrianglesGeneratedEdge.first;
    assert( equalTo(triangleThisMeshGeneratedEdge,&triangleBeingRetesselated)); // (*triangleThisMeshGeneratedEdge) == triangleBeingRetesselated );
    

    const pair<const Vertex *,const Vertex *> edge(&pairVerticesOfEdge.first,&pairVerticesOfEdge.second);
    objectsOtherMeshBoundedByEachEdge[edge] = make_pair(triangleFromOtherMeshGeneratedEdge->above,triangleFromOtherMeshGeneratedEdge->below);
  }


  //int nunNonNullNeighbors = 0, numEdgesFromIntersection = 0, numTotEdges = 0;
  for(BoundaryPolygon &polygon:boundaryPolygonsFromThisTriangle) {
    //  vector<BoundaryPolygon *> boundaryPolygonOtherSideEdge;
    //vector<pair<ObjectId,ObjectId> > objectsOtherMeshBoundedByThisEdge;
    const auto &vertexSequence = polygon.vertexSequence;
    const int numEdges = vertexSequence.size()-1;
    //numTotEdges+= numEdges;

    polygon.boundaryPolygonOtherSideEdge.resize(numEdges, NULL);
    polygon.objectsOtherMeshBoundedByThisEdge.resize(numEdges,make_pair(DONT_KNOW_ID,DONT_KNOW_ID));

    for(int i=0;i<numEdges;i++) {
       const pair<const Vertex *,const Vertex *> edgeVU(vertexSequence[i+1],vertexSequence[i]);
      //if this boundary polygon has an edge (u,v) --> the neighbor one has an edge (v,u)
      BoundaryPolygon *p = orientedEdgeToBoundaryPolygon[edgeVU];
      polygon.boundaryPolygonOtherSideEdge[i] = p;
      //nunNonNullNeighbors += (p!=NULL);

      map<pair<const Vertex *,const Vertex * >, pair<ObjectId,ObjectId> >::iterator it;
      it = objectsOtherMeshBoundedByEachEdge.find(edgeVU);
      if(it==objectsOtherMeshBoundedByEachEdge.end()) {
        it = objectsOtherMeshBoundedByEachEdge.find(make_pair(vertexSequence[i],vertexSequence[i+1])); //try to find other edge...
      }
      if(it!=objectsOtherMeshBoundedByEachEdge.end()) { //is uv (or vu) an edge from the intersection?
         polygon.objectsOtherMeshBoundedByThisEdge[i] = it->second;
         //numEdgesFromIntersection++;
      }
    }
  }


  /*#pragma omp critical
  {
    cerr << "Boundary polygons: " << boundaryPolygonsFromThisTriangle.size() << endl;
    cerr << "Tot edges: " << numTotEdges << endl;
    cerr << "Double edges from intersection: " << numEdgesFromIntersection << " " << 2*edgesFromIntersectionThisTriangle.size() << endl;
    cerr << "nunNonNullNeighbors: " << nunNonNullNeighbors << endl << endl;
  }*/
}







//edges represent the edges generated from the intersection of the intersecting triangles.
//intersectingTrianglesThatGeneratedEdges[i] contains the pair of triangles whose intersection generated edges[i]
void retesselateIntersectingTriangles(MeshIntersectionGeometry & meshIntersectionGeometry, 
                                      const vector< pair<VertexFromIntersection, VertexFromIntersection> >  &edgesFromIntersection, 
                                      const vector< pair<InputTriangle *,InputTriangle *> > &intersectingTrianglesThatGeneratedEdges,
                                      vector< pair<const InputTriangle *,vector<BoundaryPolygon>> > polygonsFromRetesselationOfEachTriangle[2]) {

  assert(edgesFromIntersection.size()==intersectingTrianglesThatGeneratedEdges.size());
  timespec t0,t1;

  #ifdef VERBOSE
  cerr << "Retesselating triangles..." << "\n";
  #endif

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

  #ifdef VERBOSE
    cerr << "Time to extract the edges intersecting each triangle: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 
  #endif

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
    polygonsFromRetesselationOfEachTriangle[meshIdToProcess].resize(numTrianglesToProcess);

    #ifdef VERBOSE
      cerr << "Number of triangles to retesselate in this mesh: " << numTrianglesToProcess << "\n";
    #endif

    statisticsAboutRetesseation.ctTrianglesRetesselate += numTrianglesToProcess;

    #pragma omp parallel 
    {
      //vector<pair<int,int> > myNewTriEdgesFromEachMap;

      

      TempVarsRetesselateTriangleFunction tempVars;

      const int percentShowLog = 1;
      int onePercentNumTriangles = (percentShowLog*numTrianglesToProcess)/100;
      if(onePercentNumTriangles==0) onePercentNumTriangles = 1;

      #pragma omp for
      for(int i=0;i<numTrianglesToProcess;i++) {      

        #ifdef VERBOSE
          if((i%onePercentNumTriangles)==0) {
            #pragma omp critical
            clog << "Retesselating " << i << " of " << numTrianglesToProcess << " Percent= " << (i*100)/numTrianglesToProcess << "\n";      
          }
        #endif  
        //for each triangle ts, let's process the triangles intersecting t and the edges
        //formed by the intersection of ts and these triangles

        const auto &ts  = *intersectingEdgesInEachTriangleToProcess[i];
        //ts.first = a triangle from map 0
        //ts.second = list of edges formed by the intersection of ts.first with other triangles...
        const auto &edgesFromIntersectionThisTriangle =  ts.second;

        polygonsFromRetesselationOfEachTriangle[meshIdToProcess][i].first = ts.first;
        vector<BoundaryPolygon> &boundaryPolygonsFromThisTriangle  = polygonsFromRetesselationOfEachTriangle[meshIdToProcess][i].second;

        
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

        /*cerr << "Edges from intersection: " << endl;
        for(auto e:edgesFromIntersection) {
          cerr << e.first.getMeshId() << e.first.getId() << endl;
          cerr << e.second.getMeshId() << e.second.getId() << endl << endl;
        }

        cerr << "Triangles intersecting this one: " << endl;
        for(int i=0;i<intersectingTrianglesThatGeneratedEdges.size();i++) {
          if(equalTo(intersectingTrianglesThatGeneratedEdges[i].first,ts.first))
            intersectingTrianglesThatGeneratedEdges[i].second->print(); cerr << endl;
        }

        cerr << "Wedge sorting... " << endl;*/
        //retesselateTriangle(edges, edgesFromIntersection,t, meshIdToProcess, myNewTrianglesFromRetesselation, tempVars, statisticsAboutRetesseation);
        retesselateTriangleUsingWedgeSorting(meshIntersectionGeometry,edgesFromIntersection,edgesFromIntersectionThisTriangle,*polygonsFromRetesselationOfEachTriangle[meshIdToProcess][i].first, boundaryPolygonsFromThisTriangle, tempVars, statisticsAboutRetesseation);
        
        //cerr << "Update polygon adjacency..." << endl;
        //each boundary polygon stores the information about the polygons adjacent to each edge and the
        //polyhedra on each side of its edges from intersection. Let's update this information now!
        updateBoundaryPolygonWithEdgeAdjacencyInformation(boundaryPolygonsFromThisTriangle,*polygonsFromRetesselationOfEachTriangle[meshIdToProcess][i].first, edgesFromIntersection,intersectingTrianglesThatGeneratedEdges,edgesFromIntersectionThisTriangle);
        //cerr << "End\n\n";
      }

      

      //break;
    }

    
    /*for(BoundaryPolygon p:polygonsFromRetesselation[meshIdToProcess]) {
      //assert(p.isInClockwiseDirection(vertices,meshIdToProcess));
    }*/
  }




  #ifdef VERBOSE
    clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "Time to retesselate creating new polygons: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n";  

    clock_gettime(CLOCK_REALTIME, &t0);
    cerr << "Triangulating polygons from retesselation...\n";
  #endif
    
  for(int meshIdToProcess=0;meshIdToProcess<2;meshIdToProcess++) {
    int numTriangles = polygonsFromRetesselationOfEachTriangle[meshIdToProcess].size();

    const int percentShowLog = 10;
    int onePercentNumPolygons = (numTriangles*percentShowLog)/100;
    if(onePercentNumPolygons==0) onePercentNumPolygons = 1;

    #pragma omp parallel
    {
      BoundaryPolygon::TempVarsTriangulatePolygon tempVarsTriangulatePolygon;

      #pragma omp for
      for(int i=0;i<numTriangles;i++) {
        //Let's triangulate the polygons on this triangle....
        #ifdef VERBOSE
          if((i%onePercentNumPolygons)==0) {
            clog << "Triangulating " << i << " of " << numTriangles << " Percent= " << (i*100)/numTriangles << "\n";
          }
        #endif
        for(BoundaryPolygon &polygon:polygonsFromRetesselationOfEachTriangle[meshIdToProcess][i].second)        
          polygon.triangulatePolygon(meshIntersectionGeometry,tempVarsTriangulatePolygon);
      }
    }
  }

  #ifdef VERBOSE
  clock_gettime(CLOCK_REALTIME, &t1);
    cerr << "Time to triangulate polygons: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 

    
    cerr << "Counts computed during retesselation: " << "\n";
    statisticsAboutRetesseation.printStatsSoFar();
    cerr << "Number of inters. tests (for insertin tris.)that are true    : " << ctEdgeIntersect << "\n";
    cerr << "Number of inters. tests (for insertin tris.)that are false    : " << ctEdgeDoNotIntersect << "\n";
  #endif


 
}