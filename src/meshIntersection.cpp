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
#include <map>
#include <cstdlib>
#include <unordered_set>
#include <unordered_map>
#include <time.h>
#include "rationals.h"
#include <omp.h>



using namespace std;


//===============================================================
// General templates (from WRF source codes)
//===============================================================

template <class T>
T max(const T a, const T b, const T c, const T d) {
  return max(max(a,b), max(c,d));
}

template <class T>
T min(const T a, const T b, const T c, const T d) {
  return min(min(a,b), min(c,d));
}

template <class T, class U>
void accum_max(T &maxsofar, const U x) {
  if (x>maxsofar) maxsofar = x;
}

template <class T>
void accum_min(T &minsofar, const T x) {
  if (x<minsofar) minsofar = x;
}

template <class T>
T square(const T x) { return x*x; }

template <class T>    // Squared length from dx and dy
T sqlen(const T dx, const T dy) { return square(dx)+square(dy); }


//From Boost, for using pairs in unordered_sets:
template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
  template<typename S, typename T> struct hash<pair<S, T>>
  {
    inline size_t operator()(const pair<S, T> & v) const
    {
      size_t seed = 0;
      ::hash_combine(seed, v.first);
      ::hash_combine(seed, v.second);
      return seed;
    }
  };
}


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

double timeReadData,timeCreateGrid,timeDetectIntersections,timeRetesselate,timeClassifyTriangles;

//===============================================================
// Constants...
//===============================================================



//===============================================================
#include "3d_objects.cpp"

#include "nested3DGrid.cpp"

#include "tritri_isectline.c"


//Each vector represents the vertices of a layer
//The first layear is mesh 0, the second is mesh1 and the third contains vertices generated from the intersection of the two meshes...
vector<Point> vertices[3]; 

//Each vector represents a set of objects in the same layer
//The objects are represented by a set of triangles (defining their boundaries)
vector<Triangle> triangles[2]; 
//vector<Triangle> trianglesFromRetesselation[2]; //triangles formed by retesselating triangles from each mesh...

Point bBox[2][2] ;//bounding box of each mesh (each boundingbox has two vertices)
Point boundingBoxTwoMeshesTogetter[2]; //bounding box considering both meshes togetter

Nested3DGridWrapper uniformGrid;

timespec t0BeginProgram, t0AfterDatasetRead;

//=======================================================================================================================


#include "PinMesh.cpp"





void extractPairsTrianglesInGridCell(const Nested3DGrid *grid,int i,int j, int k, int gridSize,vector<pair<Triangle *,Triangle *> > &pairsTrianglesToProcess) {
	int numTrianglesMesh0 = grid->numTrianglesInGridCell(0, gridSize,i,j,k);//cell.triangles[0].size();
  int numTrianglesMesh1 = grid->numTrianglesInGridCell(1, gridSize,i,j,k);



  Triangle **ptrTriMesh0 = grid->getPointerStartListTriangles(0,gridSize,i,j,k);
  Triangle **ptrTriMesh1Temp = grid->getPointerStartListTriangles(1,gridSize,i,j,k);

  for(int tA = 0;tA < numTrianglesMesh0; tA++) {
  	Triangle **ptrTriMesh1 = ptrTriMesh1Temp;
    for(int tB=0;tB<numTrianglesMesh1;tB++) {
      pairsTrianglesToProcess.push_back(pair<Triangle *,Triangle *>(*ptrTriMesh0, *ptrTriMesh1));
      ptrTriMesh1++;
    }
    ptrTriMesh0++;
  }
}

//the vector will be filled with pointers to the triangles in the same uniform grid cells
//pairs will appear maximum once in the vector
//the first element in the pair is a triangle from mesh 0 and the second one is from mesh 1
void getPairsTrianglesInSameUnifGridCells(const Nested3DGridWrapper *uniformGrid,vector<pair<Triangle *,Triangle *> > &pairsTrianglesToProcess) {
	//timespec t0,t1;
	//clock_gettime(CLOCK_REALTIME, &t0);

	pairsTrianglesToProcess.reserve(min(uniformGrid->trianglesInGrid[0]->size(),uniformGrid->trianglesInGrid[1]->size()));

  int gridSizeLevel1 =  uniformGrid->gridSizeLevel1;
	int gridSizeLevel2 =  uniformGrid->gridSizeLevel2;
  for(int i=0;i<gridSizeLevel1;i++) 
    for(int j=0;j<gridSizeLevel1;j++) 
      for(int k=0;k<gridSizeLevel1;k++) {
        if (uniformGrid->grid.hasSecondLevel(i,j,k) ) {
        	Nested3DGrid *secondLevelGrid = uniformGrid->grid.getChildGrid(i,j,k); 
			    for(int iLevel2=0;iLevel2<gridSizeLevel2;iLevel2++) 
			    	for(int jLevel2=0;jLevel2<gridSizeLevel2;jLevel2++) 
			      	for(int kLevel2=0;kLevel2<gridSizeLevel2;kLevel2++) {
			      		  extractPairsTrianglesInGridCell(secondLevelGrid,iLevel2,jLevel2,kLevel2,gridSizeLevel2,pairsTrianglesToProcess);   
			      	}    
        } else {
        	extractPairsTrianglesInGridCell(&(uniformGrid->grid),i,j,k,gridSizeLevel1,pairsTrianglesToProcess);         
        }
      }

//  clock_gettime(CLOCK_REALTIME, &t1);
  //cerr << "T before sort: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 

  sort(pairsTrianglesToProcess.begin(),pairsTrianglesToProcess.end());
  vector<pair<Triangle *,Triangle *> >::iterator it = std::unique (pairsTrianglesToProcess.begin(), pairsTrianglesToProcess.end());
  pairsTrianglesToProcess.resize( std::distance(pairsTrianglesToProcess.begin(),it) ); // 10 20 30 20 10 


 // clock_gettime(CLOCK_REALTIME, &t1);
  //cerr << "T after sort: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n";    
}


array<double,3> toDoublePoint(const Point &p) {
	array<double,3> pnew;
  for(int i=0;i<3;i++) {
    pnew[i] = p[i].get_d();
  }
  return pnew;
}

/*
void saveEdges(const string &path) {
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
*/
void storeEdgesAsGts(const string &path,const vector<pair<array<double,3>,array<double,3>> > &edgesToStore) {
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



void storeTriangleIntersections(const vector<pair<Triangle *,Triangle *> >  &pairsIntersectingTriangles) {
	cerr << "Storing triangles that intersect (and intersection edges...)" << endl;
	map<Triangle *, vector<Triangle *> > trianglesFromOtherMeshIntersectingThisTriangle[2];
	
	for(auto &p:pairsIntersectingTriangles) {
		trianglesFromOtherMeshIntersectingThisTriangle[0][p.first].push_back(p.second);
		trianglesFromOtherMeshIntersectingThisTriangle[1][p.second].push_back(p.first);
	}

	VertCoord p0InterLine[3],p1InterLine[3];
  VertCoord tempRationals[100];
  int coplanar;

  
	for(int meshId=0;meshId<2;meshId++) {
		vector<pair<array<double,3>,array<double,3>> > edgesToStore;
		int ctTriangles = 0;
		for(auto &p:trianglesFromOtherMeshIntersectingThisTriangle[meshId]) {
			Triangle *t0 = p.first;
			vector<Triangle *> &trianglesIntersectingT0 = p.second;

			
			edgesToStore.push_back(make_pair(toDoublePoint(vertices[meshId][t0->p[0]]),toDoublePoint(vertices[meshId][t0->p[1]])));
			edgesToStore.push_back(make_pair(toDoublePoint(vertices[meshId][t0->p[1]]),toDoublePoint(vertices[meshId][t0->p[2]])));
			edgesToStore.push_back(make_pair(toDoublePoint(vertices[meshId][t0->p[2]]),toDoublePoint(vertices[meshId][t0->p[0]])));
			//stringstream outputTrianglePath,outputTriangleCutPath;
			//outputTrianglePath << "out/triangle_" << meshId << "_" << ctTriangles++  << ".gts";
			//outputTriangleCutPath << "out/cut_triangle_" << meshId << "_" << ctTriangles << ".gts";

			//storeEdgesAsGts(outputTrianglePath.str(),edgesToStore);

			Triangle &a = *t0;
			for(auto t1:trianglesIntersectingT0) {				
      	Triangle &b = *t1;
      
      	int ans = tri_tri_intersect_with_isectline(vertices[meshId][a.p[0]].data(),vertices[meshId][a.p[1]].data(),vertices[meshId][a.p[2]].data()     ,     
                                                  vertices[1-meshId][b.p[0]].data(),vertices[1-meshId][b.p[1]].data(),vertices[1-meshId][b.p[2]].data(),
                                                  &coplanar, p0InterLine ,p1InterLine,tempRationals);

      	assert(ans);
      	if(!coplanar) {
      		array<double,3> p0,p1;
      		for(int i=0;i<3;i++) {
      			p0[i] = p0InterLine[i].get_d();
      			p1[i] = p1InterLine[i].get_d();
      		}
      		edgesToStore.push_back(make_pair(p0,p1));
      	} else {
      		cerr << "Coplanar found!" << endl << endl << endl;
      	}
			}
			
		}
		if(meshId==0)
			storeEdgesAsGts("out/triangles0Intersect.gts",edgesToStore);
		else 
			storeEdgesAsGts("out/triangles1Intersect.gts",edgesToStore);
	}
	
}










//returns true iff the candidate edge (represented using the ids of the vertices in meshIdToProcess -- negative vertices are in the "common layer")
//intersects an edge from the set edgesToTest (these edges are in the same layer).
//we do not consider intersections in endpoints..
bool intersects(const pair<int,int> &candidateEdge,const set<pair<int,int> > &edgesToTest,int meshIdToProcess) {
	return false;
}




//edges represent the edges generated from the intersection of the intersecting triangles.
//intersectingTrianglesThatGeneratedEdges[i] contains the pair of triangles whose intersection generated edges[i]

void retesselateIntersectingTriangles(const vector< pair< array<VertCoord,3>,array<VertCoord,3> > > &edges, const vector< pair<Triangle *,Triangle *> > &intersectingTrianglesThatGeneratedEdges) {
	assert(edges.size()==intersectingTrianglesThatGeneratedEdges.size());
	timespec t0,t1;
  cerr << "Retesselating triangles..." << "\n";
  clock_gettime(CLOCK_REALTIME, &t0);


	unordered_map<const Triangle *, vector<int> > intersectingEdgesInEachTriangle[2];
  //trianglesIntersectingEachTriangle[0] --> for each triangle t from map 0 (that intersects other triangles), 
  //trianglesIntersectingEachTriangle[0][t] is a vector of triangles from mesh 1 intersecting t...

  map< Point, int > vertexIdPlusOne; //lets add 1 to the ids (this will simplify the creation of new entries...)
  
  int numEdges = edges.size();
  vector<pair<int,int> > edgesUsingVertexId(numEdges); //if we have an edge (a,b)  --> a<b
  for(int i=0;i<numEdges;i++) {
  	//cerr << "i " << i << " " << numEdges << endl;
  	const auto &e = edges[i];
  	int idA = vertexIdPlusOne[e.first];
  	if(idA==0) {
  		idA = vertexIdPlusOne[e.first] = vertices[2].size()+1; //when we create a new vertex we start counting from 1 (only here!!!)
  		vertices[2].push_back(e.first);
  	}
  	int idB = vertexIdPlusOne[e.second];
  	if(idB==0) {
  		idB = vertexIdPlusOne[e.second] = vertices[2].size()+1;  //when we create a new vertex we start counting from 1 (only here!!!)
  		vertices[2].push_back(e.second);
  	}
  	idA = -idA; //for simplicity, let's use a negative id to referrence a "shared" vertex (from layer 2)
  	idB = -idB;
  	if(idA > idB) swap(idA,idB); 
  	edgesUsingVertexId[i] = make_pair(idA,idB); 

  	//the intersection of tA (from mesh 0) with tB generated the i-th edge...
  	const Triangle *tA = intersectingTrianglesThatGeneratedEdges[i].first;
  	const Triangle *tB = intersectingTrianglesThatGeneratedEdges[i].second;

  	intersectingEdgesInEachTriangle[0][tA].push_back(i);
  	intersectingEdgesInEachTriangle[1][tB].push_back(i);
  }

  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to extract the edges intersecting each triangle and create the new vertices: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 
 

  //edgesUsingVertexId : each edge contains the ids of its two vertices (idA,idB)
  //idA,idB start with 0 and they represent the position of the vertex in vertices[2].
  //i.e., if an edge is = (idA,idB) --> the two vertices that form this edge are vertices[2][idA],vertices[2][idB]

  //for each triangle t from mesh meshId, intersectingEdgesInEachTriangle[meshId][t] contains a vector of ids of edges formed by the intersection
  //of t with other triangles...


  //TODO: special cases, collinear, etc..

  vector<pair<int,int> > newTriEdgesFromEachMap[2];

  for(int meshIdToProcess=0;meshIdToProcess<2;meshIdToProcess++)
	  for(const auto &ts:intersectingEdgesInEachTriangle[meshIdToProcess]) {
	  	//ts.first = a triangle from map 0
	  	//ts.second = list of edges formed by the intersection of ts.first with other triangles...
	  	const auto &t = *ts.first;
	  	const auto &edgesFromIntersection =  ts.second;

	  	//we need to "retesselate" t and orient all the new triangles properly

	  	//this set will store the edges we used in the retesselated triangle
	  	//we need to choose what edges to create and, then, use these new edges to reconstruct the triangulation..
	  	set<pair<int,int> > edgesUsedInThisTriangle; //TODO: maybe use unordered_set (see overhead difference...)

	  	//The edges from intersection will, necessarelly, be in the triangulation...
	  	//edgesUsedInThisTriangle.insert(edgesFromIntersection.begin(),edgesFromIntersection.end());
	  	for(int edgeId:edgesFromIntersection) {
	  		//cerr << edgesUsingVertexId[edgeId].first << " " << edgesUsingVertexId[edgeId].second << endl;
	  		edgesUsedInThisTriangle.insert(edgesUsingVertexId[edgeId]);
	  	}

	  	vector<int> verticesToTesselate;
	  	for(const auto &p:edgesUsedInThisTriangle) {
	  		verticesToTesselate.push_back(p.first);
	  		verticesToTesselate.push_back(p.second);
	  	}
	  	verticesToTesselate.push_back(t.p[0]);
	  	verticesToTesselate.push_back(t.p[1]);
	  	verticesToTesselate.push_back(t.p[2]);
	  	sort(verticesToTesselate.begin(),verticesToTesselate.end());
	  	auto newEnd = unique(verticesToTesselate.begin(),verticesToTesselate.end());
	  	verticesToTesselate.resize(newEnd-verticesToTesselate.begin());
	  	int numVTriangle = verticesToTesselate.size();
	  	//cerr << "Num vertices to tesselate: " << numVTriangle << "\n";
	  	for(int i=0;i<numVTriangle;i++)
	  		for(int j=i+1;j<numVTriangle;j++) {
	  			int v1 = verticesToTesselate[i];
	  			int v2 = verticesToTesselate[j];
	  						
	  			assert(v1<v2); //the vector was previously sorted! also, all elements should be unique!
	  			//let's try to insert edge (v1,v2)...
	  			pair<int,int> candidateEdge(v1,v2);
	  			if(edgesUsedInThisTriangle.count(candidateEdge)!=0) continue; //the edge was already used...
	  			if(!intersects(candidateEdge,edgesUsedInThisTriangle,meshIdToProcess)) {
	  				edgesUsedInThisTriangle.insert(candidateEdge); //if it does not intersect, we can safely add this edge to the triangulation...
	  			}	  			
	  		}
	  	//edgesUsedInThisTriangle.insert(pair<int,int>(t.p[0],t.p[1]));		
	  	//edgesUsedInThisTriangle.insert(pair<int,int>(t.p[1],t.p[2]));	
	  	//edgesUsedInThisTriangle.insert(pair<int,int>(t.p[2],t.p[0]));	
	  	for(auto &elem:edgesUsedInThisTriangle)
	  		newTriEdgesFromEachMap[meshIdToProcess].push_back(elem);		 
	  } 


	clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to retesselate creating new tri edges: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n";  

  for(int meshIdToProcess=0;meshIdToProcess<2;meshIdToProcess++) {
  	vector<pair<array<double,3>,array<double,3>> > edgesToStore;

  	for(const pair<int,int> &edge:newTriEdgesFromEachMap[meshIdToProcess]) {
  		const Point *v1 = (edge.first>=0)?&vertices[meshIdToProcess][edge.first]:&vertices[2][-edge.first -1];
  		const Point *v2 = (edge.second>=0)?&vertices[meshIdToProcess][edge.second]:&vertices[2][-edge.second -1];
  		edgesToStore.push_back(make_pair(toDoublePoint(*v1),toDoublePoint(*v2)));
  	}

  	
  	if(meshIdToProcess==0)
  		storeEdgesAsGts("out/retesselatedMesh0.gts",edgesToStore);
  	else 
  		storeEdgesAsGts("out/retesselatedMesh1.gts",edgesToStore);
  }
  
}



//returns the number of intersections found
//the sets are filled with the triangles (from the corresponding mesh) that intersect
unsigned long long  computeIntersections(const Nested3DGridWrapper *uniformGrid, unordered_set<const Triangle *> trianglesThatIntersect[2]) {
  timespec t0,t1;
  clock_gettime(CLOCK_REALTIME, &t0);
   

  unsigned long long totalIntersections = 0;

  vector<pair<Triangle *,Triangle *> > vtPairsTrianglesToProcess;
  
  
  getPairsTrianglesInSameUnifGridCells(uniformGrid,vtPairsTrianglesToProcess);
  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time creating list of pairs of triangles to process (intersection): " << convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 
  //pairsTrianglesToProcess.reserve(1804900);
  

  //vector<pair<Triangle *,Triangle *> > vtPairsTrianglesToProcess(pairsTrianglesToProcess.begin(),pairsTrianglesToProcess.end());

  int numPairsToTest = vtPairsTrianglesToProcess.size();
  cerr << "Num pairs to test: " << numPairsToTest << endl;
  vector<bool> pairsIntersect(numPairsToTest,false); //true if the corresponding pair intersects..

  unsigned long long totalTests = 0;

  //TODO: try to reduce amount of computation here...

  vector< pair< array<VertCoord,3>,array<VertCoord,3> > > edges;
  vector< pair<Triangle *,Triangle *> >  intersectingTrianglesThatGeneratedEdges;
  #pragma omp parallel
  {
    VertCoord p0InterLine[3],p1InterLine[3];
    VertCoord tempRationals[100];

    vector< pair< array<VertCoord,3>,array<VertCoord,3> > > edgesTemp;
    vector< pair<Triangle *,Triangle *> >  intersectingTrianglesTemp;
    //edgesTemp.reserve(numPairsToTest/10);
    
    unsigned long long totalIntersectionsTemp = 0;
    unsigned long long totalTestsTemp = 0;

    #pragma omp for
    for(int i=0;i<numPairsToTest;i++) {   
      pair<Triangle *,Triangle *> &pairTriangles = vtPairsTrianglesToProcess[i];
      int coplanar;
      Triangle &a = *pairTriangles.first;
      Triangle &b = *pairTriangles.second;

      
      int ans = tri_tri_intersect_with_isectline(vertices[0][a.p[0]].data(),vertices[0][a.p[1]].data(),vertices[0][a.p[2]].data()     ,     
                                                  vertices[1][b.p[0]].data(),vertices[1][b.p[1]].data(),vertices[1][b.p[2]].data(),
                                                  &coplanar, p0InterLine ,p1InterLine,tempRationals);

      totalTestsTemp++;
      if (ans && !coplanar) {
        edgesTemp.push_back( pair<array<VertCoord,3>,array<VertCoord,3>>( {p0InterLine[0],p0InterLine[1],p0InterLine[2]} , {p1InterLine[0],p1InterLine[1],p1InterLine[2]}   ) );
      	intersectingTrianglesTemp.push_back(pairTriangles);
      }

      if(ans) {
        totalIntersectionsTemp++;        
        pairsIntersect[i] = true;
      } 
    }

    #pragma omp critical
    {
      totalIntersections += totalIntersectionsTemp;
      totalTests += totalTestsTemp;
      edges.insert(edges.end(), edgesTemp.begin(), edgesTemp.end());
      intersectingTrianglesThatGeneratedEdges.insert(intersectingTrianglesThatGeneratedEdges.end(), intersectingTrianglesTemp.begin(), intersectingTrianglesTemp.end());
    }
  }
  
  clock_gettime(CLOCK_REALTIME, &t1);
  timeDetectIntersections = convertTimeMsecs(diff(t0,t1))/1000; 




  clock_gettime(CLOCK_REALTIME, &t0);
  retesselateIntersectingTriangles(edges,intersectingTrianglesThatGeneratedEdges);
  clock_gettime(CLOCK_REALTIME, &t1);
  timeRetesselate = convertTimeMsecs(diff(t0,t1))/1000; 



  //Some statistics...

  for(int i=0;i<numPairsToTest;i++)  
  	if(pairsIntersect[i]) {
  		trianglesThatIntersect[0].insert(vtPairsTrianglesToProcess[i].first);
  		trianglesThatIntersect[1].insert(vtPairsTrianglesToProcess[i].second);
  	}



  //TODO: remove this from timing data (this is just a statistic)
  map<const Triangle *,int> ctIntersectionsEachTriangleFromMap[2];
  for(int i=0;i<numPairsToTest;i++)  
  	if(pairsIntersect[i]) {
  		ctIntersectionsEachTriangleFromMap[0][vtPairsTrianglesToProcess[i].first]++;
  		ctIntersectionsEachTriangleFromMap[1][vtPairsTrianglesToProcess[i].second]++;
  	}  

  int ctTrianglesIntersect[2]  = {0,0};	
  int sumIntersections[2] = {0,0};	
  int maxIntersections[2] = {0,0};
  for(int meshId=0;meshId<2;meshId++) {
  	for(auto &a:ctIntersectionsEachTriangleFromMap[meshId]) {
  		if(a.second>maxIntersections[meshId])
  			maxIntersections[meshId] = a.second;
  		ctTrianglesIntersect[meshId]++;
  		sumIntersections[meshId]+= a.second;
  	}
  	cerr << "Mesh: " << meshId << endl;
  	cerr << "Max intersections                 : " << maxIntersections[meshId] << endl;
  	cerr << "Total intersections               : " << sumIntersections[meshId] << endl;
  	cerr << "Mum triangles intersecting        : " << ctTrianglesIntersect[meshId] << endl;
  	cerr << "Average intersections per triangle: " << sumIntersections[meshId]*1.0/ctTrianglesIntersect[meshId] << endl;
  }


  cerr << "Total tests: " << totalTests << endl;
  cerr << "Total intersections: " << totalIntersections << endl;


  //for debugging/visualization, let's store in a file for each triangle the original edges of that triangle and 
  //the edges formed by the intersection of other triangles...

  vector<pair<Triangle *,Triangle *> > pairsIntersectingTriangles;
  for(int i=0;i<numPairsToTest;i++)  
  	if(pairsIntersect[i]) {
  		pairsIntersectingTriangles.push_back(vtPairsTrianglesToProcess[i]);
  	}
  storeTriangleIntersections(pairsIntersectingTriangles); // TODO: use the edges vector (we already have the edges from intersections!!!)


           
  return totalIntersections;
}





//----------------------------------------------------------------------------

/*//Each vector represents the vertices of a layer
vector<Point> vertices[2];

//Each vector represents a set of objects in the same layer
//The objects are represented by a set of triangles (defining their boundaries)
vector<Triangle> triangles[2]; //
*/
void classifyTrianglesAndGenerateOutput(const Nested3DGridWrapper *uniformGrid, const unordered_set<const Triangle *> trianglesThatIntersect[2],ostream &outputStream) {
	timespec t0,t0ThisFunction,t1;
	clock_gettime(CLOCK_REALTIME, &t0);
	t0ThisFunction = t0;

	int ctIntersectingTrianglesTotal =0;
	vector<Triangle> outputTriangles[2]; //output triangles generated from triangles of each mesh
	for(int meshId=0;meshId<2;meshId++){
		vector<Point *> verticesToLocateInOtherMesh;

		vector<int> global_x_coord_vertex_to_locate,global_y_coord_vertex_to_locate,global_z_coord_vertex_to_locate;
		for(const Triangle&t:triangles[meshId]) {
			if(trianglesThatIntersect[meshId].count(&t)==0) { //this triangle does not intersect the other mesh...
				verticesToLocateInOtherMesh.push_back(&(vertices[meshId][t.p[0]]));			
				global_x_coord_vertex_to_locate.push_back(uniformGrid->get_global_x_coord_mesh_vertex(meshId,t.p[0]));	
				global_y_coord_vertex_to_locate.push_back(uniformGrid->get_global_y_coord_mesh_vertex(meshId,t.p[0]));	
				global_z_coord_vertex_to_locate.push_back(uniformGrid->get_global_z_coord_mesh_vertex(meshId,t.p[0]));	
			} else {
				ctIntersectingTrianglesTotal++;
			}
		}
		//TODO: locate unique vertices???
		vector<ObjectId> locationOfEachVertexInOtherMesh(verticesToLocateInOtherMesh.size());
		//vertices of mesh "meshId" will be located in mesh "1-meshId"

		timespec t0,t1;
		clock_gettime(CLOCK_REALTIME, &t0);
		locateVerticesInObject(uniformGrid,  verticesToLocateInOtherMesh,global_x_coord_vertex_to_locate,global_y_coord_vertex_to_locate,global_z_coord_vertex_to_locate,locationOfEachVertexInOtherMesh,1-meshId);
		clock_gettime(CLOCK_REALTIME, &t1);
  	cerr << "Total time to locate: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
		//now, we know in what object of the other mesh each triangle that does not intersect other triangles is...


		int numTrianglesThisMesh = triangles[meshId].size();
		int ctNonIntersectingTrianglesProcessed = 0;
		for(int i=0;i<numTrianglesThisMesh;i++) {
			const Triangle &t=triangles[meshId][i];
			if(trianglesThatIntersect[meshId].count(&t)==0) { //this triangle does not intersect the other mesh...
				//this will (probably) be an output triangle...
				ObjectId objWhereTriangleIs = locationOfEachVertexInOtherMesh[ctNonIntersectingTrianglesProcessed++];		
				//cerr << "obj: " << objWhereTriangleIs << endl;		
				if (objWhereTriangleIs!=OUTSIDE_OBJECT) {
					//if the triangle is not outside the other mesh, it will be in the output (we still need to update the left/right objects correctly...)
					outputTriangles[meshId].push_back(t);
				}				
			}
			
		}		
	}
	clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to locate vertices and classify triangles: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  Print_Current_Process_Memory_Used();	
  clock_gettime(CLOCK_REALTIME, &t0); 
  
	vector<bool> verticesToWriteOutputMesh0(vertices[0].size(),false);
	vector<int> newIdVerticesMesh0(vertices[0].size(),-1);
	for(const Triangle &t : outputTriangles[0]) {
		verticesToWriteOutputMesh0[t.p[0]] = verticesToWriteOutputMesh0[t.p[1]] = verticesToWriteOutputMesh0[t.p[2]] = true;
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
	int numVerticesMesh1InOutput = 0;
	sz = vertices[1].size();
	for(int i=0;i<sz;i++) 
		if(verticesToWriteOutputMesh1[i]) {			
			newIdVerticesMesh1[i] = numVerticesMesh0InOutput + numVerticesMesh1InOutput;
			numVerticesMesh1InOutput++;
		}

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

	clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to update ids of new vertices: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  Print_Current_Process_Memory_Used();	
  clock_gettime(CLOCK_REALTIME, &t0); 

  int totalNumberOutputVertices = numVerticesMesh0InOutput+numVerticesMesh1InOutput;
  int totalNumberOutputTriangles = outputTriangles[0].size() + outputTriangles[1].size();

  
  unordered_map<pair<int,int>,int> edgesIds; //maybe use unordered_map for performance (if necessary...)
  vector<pair<int,int> > outputEdges;
  /*
  for(int meshId=0;meshId<2;meshId++) {
	  for(const Triangle &t : outputTriangles[meshId]) {
			int a = t.p[0];
			int b = t.p[1];
			int c = t.p[2];

			pair<int,int> e;
			if (a<b) {e.first = a; e.second = b;}
			else     {e.first = b; e.second = a;}
			if(edgesIds.count(e)==0) {
				outputEdges.push_back(e);
				int sz = edgesIds.size();
				edgesIds[e] = sz;
			}

			if (b<c) {e.first = b; e.second = c;}
			else     {e.first = c; e.second = b;}
			if(edgesIds.count(e)==0) {
				outputEdges.push_back(e);
				int sz = edgesIds.size();
				edgesIds[e] = sz;
			}

			if (c<a) {e.first = c; e.second = a;}
			else     {e.first = a; e.second = c;}
			if(edgesIds.count(e)==0) {
				outputEdges.push_back(e);
				int sz = edgesIds.size();
				edgesIds[e] = sz;
			}
		}
	}
	*/
	
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
			sort(myEdgesFound.begin(),myEdgesFound.end());
			auto newEndItr = unique(myEdgesFound.begin(),myEdgesFound.end());
			myEdgesFound.resize(newEndItr- myEdgesFound.begin());

			

			#pragma omp critical
			{
				outputEdges.insert(outputEdges.end(),myEdgesFound.begin(),myEdgesFound.end());
				/*for (const pair<int,int> &e:myEdgesFound) {
					//if(edgesIds.count(e)==0) {
						outputEdges.push_back(e);
					//	int sz = edgesIds.size();
					//	edgesIds[e] = sz;
					//}
				}*/
			}
		}
	}
	sort(outputEdges.begin(),outputEdges.end());
	auto newEndItr = unique(outputEdges.begin(),outputEdges.end());
	outputEdges.resize(newEndItr- outputEdges.begin());

	clock_gettime(CLOCK_REALTIME, &t1);
	cerr << "T so far: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;

	
	int totalNumberOutputEdges = outputEdges.size();
	vector<map<int,int> > mapEdgesIds2(totalNumberOutputVertices);
	for (int i=0;i<totalNumberOutputEdges;i++) {
		const pair<int,int> &e = outputEdges[i];
		//edgesIds[e] = i;
		mapEdgesIds2[e.first][e.second] = i;
	}

	clock_gettime(CLOCK_REALTIME, &t1);
	timeClassifyTriangles = convertTimeMsecs(diff(t0ThisFunction,t1))/1000;

  cerr << "Time to create edges: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  cerr << "Total time (excluding I/O) so far: " << convertTimeMsecs(diff(t0AfterDatasetRead,t1))/1000 << endl;
  Print_Current_Process_Memory_Used();	
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
	
		cerr << "Output vertices         : " << totalNumberOutputVertices << endl;
		cerr << "Output edges            : " << totalNumberOutputEdges << endl;
		cerr << "Output triangles non int: " << totalNumberOutputTriangles << endl;
		cerr << "Intersecting triangles  : " << ctIntersectingTrianglesTotal << endl;

}


void readInputMesh(int meshId, string path) {
	assert(path.size()>=4);
  string inputFileExtension = path.substr(path.size()-3,3);
  assert(inputFileExtension=="gts" || inputFileExtension=="ium");
  bool isGtsFile = inputFileExtension=="gts";
  if(isGtsFile)
    TIME(readGTSFile(path,vertices[meshId],triangles[meshId],bBox[meshId]));
  else 
    TIME(readLiumFile(path,vertices[meshId],triangles[meshId],bBox[meshId]));
}




// TODO: if(a*b > 0 ) --> do not multiply!!! use test a>>0 && b>0 || a<0 && b<0

int main(int argc, char **argv) {
  if (argc!=7) {
      cerr << "Error... use ./3dIntersection inputMesh0 inputMesh1 gridSizeLevel1 gridSizeLevel2 triggerSecondLevel outputFile.gts" << endl;
      cerr << "The mesh file may be in the gts format or in the lium format (multimaterial)" << endl;
      exit(1);
  }
  clock_gettime(CLOCK_REALTIME, &t0BeginProgram);

  string mesh0Path = argv[1];
  string mesh1Path = argv[2];
  int gridSizeLevel1 = atoi(argv[3]);
  int gridSizeLevel2 = atoi(argv[4]);
  int triggerSecondLevel = atoi(argv[5]);
  int maxTreeDepth = 2;

  Print_Current_Process_Memory_Used();
  cerr << "Reading meshes..." << endl;

  timespec t0,t1;
  clock_gettime(CLOCK_REALTIME, &t0); 

  readInputMesh(0,mesh0Path);
  readInputMesh(1,mesh1Path);

  clock_gettime(CLOCK_REALTIME, &t1);
  t0AfterDatasetRead = t1;
  cerr << "Time to read: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  timeReadData = convertTimeMsecs(diff(t0,t1))/1000; 
  Print_Current_Process_Memory_Used();


  boundingBoxTwoMeshesTogetter[0] = bBox[0][0];
  boundingBoxTwoMeshesTogetter[1] = bBox[0][1];

  for(int i=0;i<3;i++) {
  	if(bBox[1][0][i] < boundingBoxTwoMeshesTogetter[0][i])
  		boundingBoxTwoMeshesTogetter[0][i] = bBox[1][0][i];
  	if(bBox[1][1][i] > boundingBoxTwoMeshesTogetter[1][i])
  		boundingBoxTwoMeshesTogetter[1][i] = bBox[1][1][i];
  }
  

  

  

  

  for(int meshId=0;meshId<2;meshId++)
  	cerr << "Bounding box mesh " << meshId << ": " << setprecision(18) << std::fixed << bBox[meshId][0][0].get_d() << " " << bBox[meshId][0][1].get_d() <<  " " << bBox[meshId][0][2].get_d() << " -- " << bBox[meshId][1][0].get_d() << " " << bBox[meshId][1][1].get_d() << " " << bBox[meshId][1][2].get_d() <<endl;
  cerr << "Bounding box two meshes togetter " << ": " << setprecision(18) << std::fixed << boundingBoxTwoMeshesTogetter[0][0].get_d() << " " << boundingBoxTwoMeshesTogetter[0][1].get_d() <<  " " << boundingBoxTwoMeshesTogetter[0][2].get_d() << " -- " << boundingBoxTwoMeshesTogetter[1][0].get_d() << " " << boundingBoxTwoMeshesTogetter[1][1].get_d() << " " << boundingBoxTwoMeshesTogetter[1][2].get_d() <<endl;


  
  cerr <<"Creating nested grid..." << endl;

  clock_gettime(CLOCK_REALTIME, &t0); 
  vector<Triangle *> trianglesPointers[2];

  int sz = triangles[0].size();
  trianglesPointers[0].resize(sz);
  for(int i=0;i<sz;i++) trianglesPointers[0][i]  = & (triangles[0][i]);

  sz = triangles[1].size();
  trianglesPointers[1].resize(sz);
  for(int i=0;i<sz;i++) trianglesPointers[1][i]  = & (triangles[1][i]);	


//  void init(const vector<Triangle *> trianglesInsert[2], const vector<Point> vertices[2], const int gridSizeLevel1, const int gridSizeLevel2, const Point &p0, const Point &p1,const long long prodThreshold);

  TIME(uniformGrid.init(trianglesPointers, vertices, gridSizeLevel1,gridSizeLevel2,boundingBoxTwoMeshesTogetter[0],boundingBoxTwoMeshesTogetter[1],triggerSecondLevel));


  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to create and refine grid: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  timeCreateGrid = convertTimeMsecs(diff(t0,t1))/1000; 

  //After the uniform grid is initialized, let's compute the intersection between the triangles...
  cerr << "Detecting intersections..." << endl;
  clock_gettime(CLOCK_REALTIME, &t0); 
  unordered_set<const Triangle *> trianglesThatIntersect[2];
  unsigned long long numIntersectionsDetected = computeIntersections(&uniformGrid,trianglesThatIntersect);

  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to detect intersections (includes time for computing statistics and for saving intersections for debugging): " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  Print_Current_Process_Memory_Used();



  clock_gettime(CLOCK_REALTIME, &t0); 

  ofstream outputStream(argv[6]);
  assert(outputStream);
  classifyTrianglesAndGenerateOutput(&uniformGrid,trianglesThatIntersect,outputStream);

  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Total time to classify triangles and generate output: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  Print_Current_Process_Memory_Used();

  cerr << "----------------------------------------------------" << endl;
  cerr << "Summary of ACTUAL times (excluding the times to compute statistics, to write edges for debugging, etc): " << endl;
  cerr << "Time to read the data         : " << timeReadData << endl;
  cerr << "Time to create and refine grid: " << timeCreateGrid << endl;
  cerr << "Time to detect intersections  : " << timeDetectIntersections << endl;  
  cerr << "Time to classify the triangles: " << timeClassifyTriangles << endl;
  cerr << "----------------------------------------------------" << endl;
/*
  if(isGtsFile)
    for(const ObjectId& id:pointIds) {
      cout << id << "\n";
    }
  else
    for(const ObjectId& id:pointIds) {  //when we deal with LIUM files the -1 represents the outside object and objects start with id 0.... (we convert from -1 to 0 and from 0 to 1 in the input... so, we need to revert this change...s)
        cout << id-1 << "\n";
    } 
*/
}