#include "3dGeometry.h"





#include "nested3DGrid.h"
#include <sstream>


#ifdef COLLECT_GEOMETRY_STATISTICS
struct GeometryStatistics {
  int ctDegeneraciesIntersectTwoTriangles =0,
      ctDegeneraciesIntersectTwoTrianglesCoplanar= 0,
            ctDegeneraciesIsCloser=0,
            ctDegeneraciesIsAngleWith0Greater=0,
            ctDegeneraciesIsVertexInTriangleProjection=0,
            ctDegeneraciesIsVertexConvex=0,
            ctDegeneraciesIsVertexInInputTriangleProjection=0,
            ctDegeneraciesIsTriangleNormalPointingPositiveZ=0,
            ctDegeneraciesIsTriangleAbovePointSoS=0,
            ctDegeneraciesGetBestTrianglePointInObjectSoS=0;

  //orientation2D_original_original_original
  //orientation2D_original_original_intersect
  //....
  int orientation1DOO=0,orientation1DOI=0,orientation1DII=0;
  int orientation2DOOO = 0, orientation2DOOI = 0, orientation2DOII = 0, orientation2DIII = 0;
  int orientation3DOO = 0, orientation3DOI = 0 ;
  int signVector = 0;

  void printStatsSoSCalls() {
    cerr << "SoS calls: \n" ;
    cerr << "Orientation 1D:" << "\n";
    cerr << "OO :"  << orientation1DOO << "\n";
    cerr << "OI :"  << orientation1DOI << "\n";
    cerr << "II :"  << orientation1DII << "\n";
    cerr << "Orientation 2D:"  <<"\n";
    cerr << "OOO: " << orientation2DOOO <<"\n";
    cerr << "OOI: " << orientation2DOOI <<"\n";
    cerr << "OII: " << orientation2DOII <<"\n";
    cerr << "III: " << orientation2DIII <<"\n";
    cerr << "Orientation 3D:"  <<"\n";
    cerr << "OO: " << orientation3DOO <<"\n";
    cerr << "OI: " << orientation3DOI <<"\n";
    cerr << "Sign of vector coord: " << signVector << "\n";
  }

   void printStats() {
     cerr << "ctDegeneraciesIntersectTwoTriangles: " << ctDegeneraciesIntersectTwoTriangles<< "\n";
     cerr << "ctDegeneraciesIntersectTwoTrianglesCoplanar: " << ctDegeneraciesIntersectTwoTrianglesCoplanar << "\n";
     cerr << "ctDegeneraciesIsCloser: " << ctDegeneraciesIsCloser<< "\n";
     cerr << "ctDegeneraciesIsAngleWith0Greater: " << ctDegeneraciesIsAngleWith0Greater<< "\n";
     cerr << "ctDegeneraciesIsVertexInTriangleProjection: " << ctDegeneraciesIsVertexInTriangleProjection<< "\n";
     cerr << "ctDegeneraciesIsVertexConvex: " << ctDegeneraciesIsVertexConvex << "\n";
     cerr << "ctDegeneraciesIsVertexInInputTriangleProjection: " << ctDegeneraciesIsVertexInInputTriangleProjection<< "\n";
     cerr << "ctDegeneraciesIsTriangleNormalPointingPositiveZ: " << ctDegeneraciesIsTriangleNormalPointingPositiveZ<< "\n";
     cerr << "ctDegeneraciesIsTriangleAbovePointSoS: " << ctDegeneraciesIsTriangleAbovePointSoS<< "\n";     
     cerr << "ctDegeneraciesGetBestTrianglePointInObjectSoS: " << ctDegeneraciesGetBestTrianglePointInObjectSoS<< "\n";  
   }
};
GeometryStatistics geometryStatisticsDegenerateCases, geometryStatisticsNonDegenerateCases;
#endif



#include "geometry/3dGeometryGeometricalPredicatesSoSImpl.cpp"
#include "geometry/3dGeometryGeometricalPredicatesMainImpl.cpp"
#include "geometry/3dGeometryGeometricalPredicatesMainImplOrig.cpp"
#include "geometry/3dGeometryGeometricalPredicates.cpp"



#include "boundaryPolygon.cpp"





void MeshIntersectionGeometry::sortEdgesSharingStartingVertexByAngle(vector<pair<const Vertex *,const Vertex *> >::iterator begin,
                                            vector<pair<const Vertex *,const Vertex *> >::iterator end,
                                            const int planeProjectTriangleTo,  TempVarsSortEdgesByAngle &tempVars)  {

  

  //sort(begin,end, [&](const pair<const Vertex *,const Vertex *> &e1, const pair<const Vertex *,const Vertex *> &e2) {
  //          return meshIntersectionGeometry.isAngleWith0Greater(*e1.first, *e1.second, *e2.second, planeProjectTriangleTo, tempVars.tempVarsIsAngleWith0Greater);
  //      } );
  //cerr << "Sorting..." << endl;
  bool hasDegenerate = false;
  vector<pair<const Vertex *,const Vertex *> >::iterator firstNonZero = begin;

  vector<pair<const Vertex *,const Vertex *> >::iterator it = begin;
  while(it!=end) {
    int onZeroPlusAxis = isOnZeroPlusAxisNoSoS(*it->first,*it->second,planeProjectTriangleTo,tempVars.tempVarsIsOnZeroPlusAxisNoSoS);
    if(onZeroPlusAxis==1) {
      //vertex pointed by it is on zeroplus axis...
      swap(*it,*firstNonZero);
      firstNonZero++;
    } else if(onZeroPlusAxis ==0) {
      hasDegenerate = true;
      break;
    }

    it++;
  }
  //cerr << "End sorting" << endl;

  if(hasDegenerate) {
    //we will need an special function to sort everything...
    //a degenerate edge may be anywhere...

    //if this happens, let's just call the SoS implementation
    //cerr << "Degenerate" << endl;
    //cerr << "Plane: " << planeProjectTriangleTo << endl;
    sort(begin,end, [&](const pair<const Vertex *,const Vertex *> &e1, const pair<const Vertex *,const Vertex *> &e2) {
      return 1==isAngleWith0GreaterSoSImpl(*e1.first, *e1.second, *e2.second, planeProjectTriangleTo, tempVars.tempVarsIsAngleWith0Greater);
    } );
    //cerr << "end Degenerate" << endl;
    
  } else {
    //we do not have any degenerate edge

    //cerr << "non degenerate..." << endl;
    //the vertices [begin,firstNonZero) are on the zero+ axis
    //we can simply sort them by orientation using SoS (since they are coincident)
    //since they are all coincident, we do not have to care with the axis (just with the orientation)
    //so, everything is sorted now (but we are not considering SoS...)
    //thus, we have to sort everything that is coincident using SoS

    //first, let's sort the edges with angle 0...
    //we have to sort them using SoS...
    sort(begin,firstNonZero, [&](const pair<const Vertex *,const Vertex *> &e1, const pair<const Vertex *,const Vertex *> &e2) {
        return isOrientationPositiveSoSImpl(*e1.first, *e1.second, *e2.second, planeProjectTriangleTo, tempVars.tempVarsIsAngleWith0Greater);
      } );

    #ifdef COLLECT_GEOMETRY_STATISTICS
      int numDegenerates = 0;
      if(firstNonZero-begin > 1) numDegenerates = firstNonZero-begin-1;
      #pragma omp atomic
        geometryStatisticsDegenerateCases.ctDegeneraciesIsAngleWith0Greater+= numDegenerates;      
    #endif



    //the vertices [firstNonZero,end) are not on the zero+ axis
    //we can sort them using the non-SoS version of the function...
    sort(firstNonZero,end, [&](const pair<const Vertex *,const Vertex *> &e1, const pair<const Vertex *,const Vertex *> &e2) {
        return isAngleWith0GreaterNoSoSNonZeroAngle(*e1.first, *e1.second, *e2.second, planeProjectTriangleTo, tempVars.tempVarsIsAngleWith0Greater);
      } );

  
    //now, we have to sort the coincident vertices (that don't have angle 0)
    vector<pair<const Vertex *,const Vertex *> >::iterator itB = firstNonZero;
    vector<pair<const Vertex *,const Vertex *> >::iterator itE;
    while(itB!=end) {
      itE = itB+1;    

      //check if there is a coincidence...
      while(itE!=end && 0==isAngleWith0GreaterNoSoSNonZeroAngle(*(itB->first), *(itB->second),*(itE->second), planeProjectTriangleTo, tempVars.tempVarsIsAngleWith0Greater)){
        itE++;
      }

      //edges in [itB,itE) have equal angles... we have to sort them by orientation...
      sort(itB,itE, [&](const pair<const Vertex *,const Vertex *> &e1, const pair<const Vertex *,const Vertex *> &e2) {
        return isOrientationPositiveSoSImpl(*e1.first, *e1.second, *e2.second, planeProjectTriangleTo, tempVars.tempVarsIsAngleWith0Greater);
      } );

      #ifdef COLLECT_GEOMETRY_STATISTICS
        int numDegenerates = 0;
        if(itE-itB > 1) numDegenerates = itE-itB-1;
        #pragma omp atomic
          geometryStatisticsDegenerateCases.ctDegeneraciesIsAngleWith0Greater+= numDegenerates;      
      #endif

      itB = itE;
    }

    //cerr << "end non degenerate..." << endl;
  }


  #ifdef DEBUGGING_MODE
    int nv = end-begin;

    int ctDiff = 0;
    for(int i=0;i<nv;i++) {
        const pair<const Vertex *,const Vertex *> &e1 = *(begin+i);
        const pair<const Vertex *,const Vertex *> &e2 = *(begin+(i+1)%nv);
        int isAngle = isAngleWith0GreaterSoSImpl(*e1.first, *e1.second, *e2.second, planeProjectTriangleTo, tempVars.tempVarsIsAngleWith0Greater);
        ctDiff += isAngle==0;
      } 


    if(ctDiff!=1) {
      #pragma omp critical
      {
        for(int i=0;i<10;i++)
          cerr << "======================================================" << endl;
        cerr << "Comparing: ";      
        for(int i=0;i<nv;i++) {
          const pair<const Vertex *,const Vertex *> &e1 = *(begin+i);
          const pair<const Vertex *,const Vertex *> &e2 = *(begin+(i+1)%nv);
          int isAngle = isAngleWith0GreaterSoSImpl(*e1.first, *e1.second, *e2.second, planeProjectTriangleTo, tempVars.tempVarsIsAngleWith0Greater);
          cerr << isAngle << " ";
        } cerr << endl;
      }
    }
  #endif
}














MeshIntersectionGeometry::~MeshIntersectionGeometry() {
  #ifdef COLLECT_GEOMETRY_STATISTICS
    cerr << "\nCounts degenerate cases:\n";
    geometryStatisticsDegenerateCases.printStats();
    cerr << "\nCounts non-degenerate cases:\n";
    geometryStatisticsNonDegenerateCases.printStats();
    cerr << "\nCounts SoS calls: \n";
    geometryStatisticsDegenerateCases.printStatsSoSCalls();

    cerr << "\nFreeing memory..." << endl;
  #endif
}

const MeshIntersectionGeometry::PlaneEquation &MeshIntersectionGeometry::getPlaneEquationInputTriangle(int meshId, int triId,TempVarsComputePlaneEquation &tempVars) {
	assert(meshId>=0 && meshId<=1);
  assert(triId>=0 && triId<planeEquationsInputTriangles[meshId].size());
  assert(isPlaneEquationInputTrianglesInitialized[meshId][triId]);
  //initPlaneEquationInputTriangle(meshId,triId,tempVars);
	return planeEquationsInputTriangles[meshId][triId];
}

void MeshIntersectionGeometry::initPlaneEquationInputTriangle(int meshId, int triId,TempVarsComputePlaneEquation &tempVars) {
  if(!isPlaneEquationInputTrianglesInitialized[meshId][triId]) {
    const InputTriangle &t = inputTriangles[meshId][triId];
    computePlaneEquation(planeEquationsInputTriangles[meshId][triId], getCoordinates(*t.getInputVertex(0)), getCoordinates(*t.getInputVertex(1)),getCoordinates(*t.getInputVertex(2)), tempVars);
    isPlaneEquationInputTrianglesInitialized[meshId][triId] = true;
  }
}



void MeshIntersectionGeometry::storeIntersectionVerticesCoordinatesAndUpdateVerticesIds(vector< pair<VertexFromIntersection, VertexFromIntersection> >  &edgesFromIntersection,const vector< pair<Point, Point> > &coordsVerticesOfEdges, const vector< pair<InputTriangle *,InputTriangle *> >  &intersectingTrianglesThatGeneratedEdges) {
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
  cerr << "After assigning vertices ids (for rational coords): " <<  convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 
  clock_gettime(CLOCK_REALTIME, &t0);


  //TODO: maybe it would be better to use an unordered map or to sort and binary search...
  struct VertInterComparator {
    bool operator()(const VertexFromIntersection &a,const VertexFromIntersection &b)  const {
        return a.compare(b)<0;
    }
  };
  map<VertexFromIntersection,int,VertInterComparator> idsOfVerticesFromIntersection;
  for(int i=0;i<numEdges;i++) {
    VertexFromIntersection &v1 = edgesFromIntersection[i].first;
    VertexFromIntersection &v2 = edgesFromIntersection[i].second;
    int idV1,idV2;

    if(idsOfVerticesFromIntersection.count(v1)==0) {
      idV1 = idsOfVerticesFromIntersection.size();
      idsOfVerticesFromIntersection[v1] = idV1;
    } else {
      idV1 = idsOfVerticesFromIntersection[v1];
    } 
    if(idsOfVerticesFromIntersection.count(v2)==0) {
      idV2 = idsOfVerticesFromIntersection.size();
      idsOfVerticesFromIntersection[v2] = idV2;
    } else {
      idV2 = idsOfVerticesFromIntersection[v2];
    }

    v1.idForEps = idV1;
    v2.idForEps = idV2;
  }
  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "After assigning vertices ids (for eps coords): " <<  convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 
  cerr << "Num unique vertices from intersection: " << idsOfVerticesFromIntersection.size() << "\n";


  clock_gettime(CLOCK_REALTIME, &t0);


  int numPairsInterTri = intersectingTrianglesThatGeneratedEdges.size();
  int numUniqueInterTriSeenSoFar = 0;

  struct InputTriComparator {
    bool operator()(const InputTriangle &a,const InputTriangle &b)  const {
        return a.compare(b)<0;
    }
  };
  map<InputTriangle,int,InputTriComparator> idForEpsEachTriangle;
  for(int i=0;i<numPairsInterTri;i++) { //only triangles intersecting other triangles need to have their locks initialized...
    int triId = intersectingTrianglesThatGeneratedEdges[i].first->id;
    if(triId==-1) { //not initialized yet...
      intersectingTrianglesThatGeneratedEdges[i].first->setIdForEps(numUniqueInterTriSeenSoFar);
      idForEpsEachTriangle[*intersectingTrianglesThatGeneratedEdges[i].first] = numUniqueInterTriSeenSoFar;
      numUniqueInterTriSeenSoFar++;
    }

    triId = intersectingTrianglesThatGeneratedEdges[i].second->id;
    if(triId==-1) { //not initialized yet...
      intersectingTrianglesThatGeneratedEdges[i].second->setIdForEps(numUniqueInterTriSeenSoFar);
      idForEpsEachTriangle[*intersectingTrianglesThatGeneratedEdges[i].second] = numUniqueInterTriSeenSoFar;
      numUniqueInterTriSeenSoFar++;
    }  
  }
  
  //unfortunately in our current implementation the vertices from intersection do not store pointers to input triangles, 
  //it would be be better to use pointers and keep the instances unique...
  for(int i=0;i<numEdges;i++) {
    VertexFromIntersection &v1 = edgesFromIntersection[i].first;
    VertexFromIntersection &v2 = edgesFromIntersection[i].second;
    v1.triangle.id = idForEpsEachTriangle[v1.triangle];
    v2.triangle.id = idForEpsEachTriangle[v2.triangle];
  }


  cerr << "Number of unique triangles: " << numUniqueInterTriSeenSoFar << endl;

  tEpsDeterminanCommonTerms.resize(numUniqueInterTriSeenSoFar);
  epsCoefficientsVertexIntersection.resize(idsOfVerticesFromIntersection.size() );



  /*int n = tEpsDeterminanCommonTerms.size();
  #pragma omp parallel for
  for(int i=0;i<n;i++) {
    if(tEpsDeterminanCommonTerms[i].willBeUsed) //we only need to init the lock if the triangle intersects other triangles...
      omp_init_lock(&(tEpsDeterminanCommonTerms[i].lock));
  }*/

  int n = tEpsDeterminanCommonTerms.size();
  #pragma omp parallel for
  for(int i=0;i<n;i++) {
    tEpsDeterminanCommonTerms[i].init = false;
    omp_init_lock(&(tEpsDeterminanCommonTerms[i].lock));  
  }

  n = epsCoefficientsVertexIntersection.size();
  #pragma omp parallel for
  for(int i=0;i<n;i++) {
    epsCoefficientsVertexIntersection[i].init = false;
  }

  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to init vectors for numerator/denominator for t: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
}

void MeshIntersectionGeometry::computeIntersections(const vector<pair<InputTriangle *,InputTriangle *> > &inputTrianglesToConsider, 
                          vector< pair<InputTriangle *,InputTriangle *> >  &intersectingTrianglesThatGeneratedEdges, 
                          vector< pair<VertexFromIntersection, VertexFromIntersection> >  &edgesFromIntersection, 
                          unsigned long long &numIntersectionTests){ 

  const long long numPairsTrianglesToProcess = inputTrianglesToConsider.size();
  numIntersectionTests = numPairsTrianglesToProcess;

  timespec t0,t1;
  clock_gettime(CLOCK_REALTIME, &t0);

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

    vector<int> trianglesToConsider(inputTriangles[meshId].size(),0);
    if(meshId==0) {
      #pragma omp parallel for
      for(int i=0;i<numPairsTrianglesToProcess;i++) {
        int id = inputTrianglesToConsider[i].first-&inputTriangles[0][0];
        trianglesToConsider[id] = 1;
      }
    } else {
      #pragma omp parallel for
      for(int i=0;i<numPairsTrianglesToProcess;i++) {
        int id = inputTrianglesToConsider[i].second-&inputTriangles[1][0];
        trianglesToConsider[id] = 1;
      }
    }


    int numTri = inputTriangles[meshId].size();
    cerr << "Triangles to compute equations determined\n";

  
    #pragma omp parallel
    {
      TempVarsComputePlaneEquation tempVars;

      #pragma omp for
      for(int i=0;i<numTri;i++) {
        if(trianglesToConsider[i])
          initPlaneEquationInputTriangle(meshId, i,tempVars);
      }
    }

  }
  //for(int i:isPlaneEquationInputTrianglesInitialized[1]) cerr << i << " "; cerr << endl;
  cerr << "Computing the intersections...\n"; 

  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to only compute plane equations: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 

  clock_gettime(CLOCK_REALTIME, &t0);

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

  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to only compute intersections: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 

  cerr << "Registering new vertices" << endl;
  //After the intersections are computed, we need to "register" the coordinates of the new vertices and store their ids...
  storeIntersectionVerticesCoordinatesAndUpdateVerticesIds(edgesFromIntersection,coordsVerticesOfEdges,intersectingTrianglesThatGeneratedEdges);


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
  assert(inputFileExtension=="gts" || inputFileExtension=="ium" || inputFileExtension=="off");
  bool isGtsFile = inputFileExtension=="gts";
  bool isOffFile = inputFileExtension=="off";

  int numTrianglesPreviouslyRead = (meshId==0)?0:inputTriangles[0].size();
  if(isGtsFile)
    readGTSFile(path,verticesCoordinates[meshId],inputTriangles[meshId],meshBoundingBoxes[meshId],meshId,numTrianglesPreviouslyRead);
  else if (isOffFile) 
    readOFFFile(path,verticesCoordinates[meshId],inputTriangles[meshId],meshBoundingBoxes[meshId],meshId,numTrianglesPreviouslyRead);
  else 
    readLiumFile(path,verticesCoordinates[meshId],inputTriangles[meshId],meshBoundingBoxes[meshId],meshId,numTrianglesPreviouslyRead);

  cerr << "Initializing bounding-boxes of triangles\n"; 
  timespec t0,t1;
  clock_gettime(CLOCK_REALTIME, &t0);

  for(int meshId=0;meshId<2;meshId++) {
  	const int numInputTrianglesThisMesh = inputTriangles[meshId].size();

  	inputTrianglesBoundingBox[meshId].resize(numInputTrianglesThisMesh);
  	#pragma omp parallel for
  	for(int tid=0;tid<numInputTrianglesThisMesh;tid++) {
  		initTriangleBoundingBox(meshId,tid);
  	}

  	planeEquationsInputTriangles[meshId].resize(numInputTrianglesThisMesh);
  	isPlaneEquationInputTrianglesInitialized[meshId] = vector<int>(numInputTrianglesThisMesh,0);
  }
  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to init bounding boxes: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
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
void readGTSFile(string fileName, vector<Point> &vertices,vector<InputTriangle> &triangles, Point boundingBox[2],const int meshId,const int numTrianglesPreviouslyRead) {
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


/*
//Parsing rationals exactly
//source: https://mortoray.com/2013/03/14/parsing-an-exact-decimal-value-using-gmp/
mpq_class parse_decimal( std::string const & str )
{
    //meta = pdf_is_exact;
 
    //this is intentionally far more lenient than the one in node_parser to retain flexibility
    static boost::regex re_num( "([+-]?)([0-9]*)(\\.([0-9]*))?([eE]([+-]?)([0-9]+))?" );
    boost::match_results<std::string::const_iterator> what;
    if( !boost::regex_match( str, what, re_num ) )
        PRE_FAIL( std::string( "invalid number: " + str ) );
 
    mpz_class part;
    mpq_class combine;
    bool negative = false;
 
    if( what[1].matched && *what[1].first == '-' )
        negative = true;
 
    if( what[2].matched )
        combine = parse_bitz( what[2].first, what[2].second );
 
    if( what[3].matched ) {
        //meta |= pdf_had_decimal;
        unsigned len = unsigned( what[4].second - what[4].first );
        if( len ) {
            part = parse_bitz( what[4].first, what[4].second );
            mpz_class div(10);
            mpz_pow_ui( div.get_mpz_t(), div.get_mpz_t(), len );
            combine += mpq_class( part, div );
        }
    }
 
    if( what[5].matched ) {
       // meta |= pdf_had_exponent;
        mpz_class e = parse_bitz( what[7].first, what[7].second );
        mpz_class base(10), base_exp;
 
        //TODO: check that "e" is constrained to unsigned long
        mpz_pow_ui( base_exp.get_mpz_t(), base.get_mpz_t(), e.get_ui() );
 
        if( *what[6].first == '-' )
            combine /= base_exp;
        else
            combine *= base_exp;
    }
 
    if( negative )
        combine *= -1;
    combine.canonicalize();
    return number( combine );
}
*/

void readOFFFile(string fileName, vector<Point> &vertices,vector<InputTriangle> &triangles, Point boundingBox[2],const int meshId,const int numTrianglesPreviouslyRead) {
  ifstream fin(fileName.c_str());


  if (!fin) {
    cerr << "ERROR: failed to open file " << fileName << endl;
    exit(1);
  }
  cerr << "Reading file " << fileName << endl;
  int numVertices, numTriangles,numObjects;
  string off;
  getline(fin,off);
  if(off!="OFF") {
    cerr << "Warning: Is this file really an OFF file? " << off << endl;
  }

  string temp;
  //fin.ignore('\n');
  while(getline(fin,temp)) {
    cerr <<  "tmep: " << temp << endl;
    if(temp.find("#")==std::string::npos) break;
  }
  stringstream ss;
  ss<< temp;
  cerr << temp << endl;
  assert(ss >> numVertices >> numTriangles >> numObjects);



  cerr << "Vertices, triangles, objects: " << numVertices << " " << numTriangles << " " << numObjects << endl;
  vertices.resize(numVertices);
  triangles.resize(numTriangles);

  //rational xr,yr,zr;
  for(int i=0;i<numVertices;i++) {
    getline(fin,temp);
    stringstream ss;
    ss<< temp;
    /*if(!(ss >> vertices[i][0] >> vertices[i][1] >> vertices[i][2])) {
      cerr << "Error reading vertex coord:" << temp << endl;
      exit(1);
    }*/
    double x,y,z;
    ss >> x >> y >> z;

    //ss << temp;
    //ss >> vertices[i][0] >> vertices[i][1] >> vertices[i][2];
    //cerr << x << " " << vertices[i][0].get_d() << " " << vertices[i][0] << endl;


    //double x,y,z;
    //assert(fscanf(inputFile,"%lf %lf %lf",&x,&y,&z)==3);

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
    int nv;
    
    getline(fin,temp);
    stringstream ss;
    ss<< temp;
    if(!(ss >> nv >> a >> b >> c)) {
      cerr << "Error reading: " << temp << endl;
      exit(1);
    }
    assert(nv==3);

    //assert(fscanf(inputFile,"%d %d %d %d %d",&a,&b,&c,&objNormal,&objContrNormal)==5);
    //a--;b--;c--; //we count from 0, not from 1...
    

    //cerr << objNormal << " " << objContrNormal << endl;
    triangles[i] = InputTriangle(InputVertex(meshId,a),InputVertex(meshId,b),InputVertex(meshId,c),OUTSIDE_OBJECT,1);
  }
}


//Reads a Lium file, fills the boundingBox with the boundingBox of the triangles read
void readLiumFile(string fileName, vector<Point> &vertices,vector<InputTriangle> &triangles, Point boundingBox[2],const int meshId,const int numTrianglesPreviouslyRead) {
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
      out << std::setprecision (std::numeric_limits<double>::digits10 + 2) << p[0].get_d() << " " << p[1].get_d() << " " << p[2].get_d() << "\n";    
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


void MeshIntersectionGeometry::saveEdgesAsGTS(const vector<pair<const Vertex *,const Vertex *>>  &edges,const string &path)  {
  vector<pair<array<double,3>,array<double,3>> > edgesToStore;
  for(const pair<const Vertex *,const Vertex *> &edge:edges) {

    array<double,3> v0 = getCoordinatesForDebugging(*edge.first);
    array<double,3> v1 = getCoordinatesForDebugging(*edge.second);
    edgesToStore.push_back({v0,v1});
  }

  storeEdgesAsGts(path,edgesToStore );
}

void MeshIntersectionGeometry::storeEdgesAsGts(const string &path,const vector<pair<array<double,3>,array<double,3>> > &edgesToStore)  {
  /*map<array<double,3>, int> vertexToId;

  vector<array<double,3> > vertices;
  for(auto &e:edgesToStore) {
    if(vertexToId.count(e.first)==0) {
      int id = vertexToId.size();
      vertexToId[e.first] = id;
      vertices.push_back(e.first);
    }
    if(vertexToId.count(e.second)==0) {
      int id = vertexToId.size();
      vertexToId[e.second] = id;
      vertices.push_back(e.second);
    }
  }*/

  ofstream fout(path.c_str());
  int numVert = 3*edgesToStore.size();
  int numEdges = 3*edgesToStore.size();
  int numFaces = edgesToStore.size();

  fout << numVert << " " << numEdges << " " << numFaces << "\n";
  for(pair<array<double,3>,array<double,3>> e:edgesToStore) {
    fout << e.first[0] << " " << e.first[1] << " " << e.first[2] << "\n";
    fout << e.second[0] << " " << e.second[1] << " " << e.second[2] << "\n";
    fout << e.second[0] << " " << e.second[1] << " " << e.second[2] << "\n";
  }
  int start = 1;
  for(int i=0;i<numFaces;i++) {
    fout << start << " " << start+1 << "\n";
    fout << start+1 << " " << start+2 << "\n";
    fout << start+2 << " " << start << "\n";
    start+= 3;
  }  
  start = 1;
  for(int i=0;i<numFaces;i++) {
    fout << start << " " << start+1 << " " << start+2 <<  "\n";
    start+= 3;
  }
}