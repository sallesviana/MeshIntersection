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
#include <set>
#include <queue>
#include <sstream>
#include <time.h>
#include "rationals.h"
#include "utils.h"
#include "nested3DGrid.h"
#include "common2.h"
#include <omp.h>
#include <parallel/algorithm>


using namespace std;


//===============================================================
// Constants, flags...
//===============================================================

//#define DEBUGGING_MODE

//This may slow down the algorithm
#define COLLECT_STATISTICS
//#define COLLECT_STATISTICS_PRINT_TRIANGLES_INTERSECTIONS

//#define SANITY_CHECKS

//===============================================================





double timeReadData,timeCreateGrid,timeDetectIntersections,timeRetesselate,timeClassifyTriangles;





Nested3DGridWrapper uniformGrid;

timespec t0BeginProgram, t0AfterDatasetRead;


//=======================================================================================================================


void extractPairsTrianglesInGridCell(const Nested3DGrid *grid,int i,int j, int k, int gridSize,vector<pair<InputTriangle *,InputTriangle *> > &pairsTrianglesToProcess) {
  int numTrianglesMesh0 = grid->numTrianglesInGridCell(0, gridSize,i,j,k);//cell.triangles[0].size();
  int numTrianglesMesh1 = grid->numTrianglesInGridCell(1, gridSize,i,j,k);



  InputTriangle **ptrTriMesh0 = grid->getPointerStartListTriangles(0,gridSize,i,j,k);
  InputTriangle **ptrTriMesh1Temp = grid->getPointerStartListTriangles(1,gridSize,i,j,k);

  for(int tA = 0;tA < numTrianglesMesh0; tA++) {
    InputTriangle **ptrTriMesh1 = ptrTriMesh1Temp;
    for(int tB=0;tB<numTrianglesMesh1;tB++) {
      pairsTrianglesToProcess.push_back(pair<InputTriangle *,InputTriangle *>(*ptrTriMesh0, *ptrTriMesh1));
      ptrTriMesh1++;
    }
    ptrTriMesh0++;
  }

}


unsigned long long sumPairsTrianglesInGridCell(const Nested3DGrid *grid,int i,int j, int k, int gridSize,vector<pair<InputTriangle *,InputTriangle *> > &pairsTrianglesToProcess) {
  int numTrianglesMesh0 = grid->numTrianglesInGridCell(0, gridSize,i,j,k);//cell.triangles[0].size();
  int numTrianglesMesh1 = grid->numTrianglesInGridCell(1, gridSize,i,j,k);

  return ((unsigned long long)numTrianglesMesh0)*numTrianglesMesh1;

}

unsigned long long computeNumPairsTrianglesToProcessBeforeUnique(const Nested3DGridWrapper *uniformGrid,vector<pair<InputTriangle *,InputTriangle *> > &pairsTrianglesToProcess) {
  int gridSizeLevel1 =  uniformGrid->gridSizeLevel1;
  int gridSizeLevel2 =  uniformGrid->gridSizeLevel2;

  unsigned long long ans = 0;
  for(int i=0;i<gridSizeLevel1;i++) 
    for(int j=0;j<gridSizeLevel1;j++) 
      for(int k=0;k<gridSizeLevel1;k++) {
        if (uniformGrid->grid.hasSecondLevel(i,j,k) ) {
          Nested3DGrid *secondLevelGrid = uniformGrid->grid.getChildGrid(i,j,k); 
          for(int iLevel2=0;iLevel2<gridSizeLevel2;iLevel2++) 
            for(int jLevel2=0;jLevel2<gridSizeLevel2;jLevel2++) 
              for(int kLevel2=0;kLevel2<gridSizeLevel2;kLevel2++) {
                  ans += sumPairsTrianglesInGridCell(secondLevelGrid,iLevel2,jLevel2,kLevel2,gridSizeLevel2,pairsTrianglesToProcess);   
              }    
        } else {
          ans += sumPairsTrianglesInGridCell(&(uniformGrid->grid),i,j,k,gridSizeLevel1,pairsTrianglesToProcess);         
        }
      }
  return ans;
}

//the vector will be filled with pointers to the triangles in the same uniform grid cells
//pairs will appear maximum once in the vector
//the first element in the pair is a triangle from mesh 0 and the second one is from mesh 1
void getPairsTrianglesInSameUnifGridCells(const Nested3DGridWrapper *uniformGrid,vector<pair<InputTriangle *,InputTriangle *> > &pairsTrianglesToProcess) {
  //timespec t0,t1;
  //clock_gettime(CLOCK_REALTIME, &t0);

  //TODO: remove this for performance purposess.... just for debugging
  //cerr << "Num pairs to process according to this configuration: " << computeNumPairsTrianglesToProcessBeforeUnique(uniformGrid,pairsTrianglesToProcess) << "\n";
  //exit(0);

  pairsTrianglesToProcess.reserve(min(uniformGrid->trianglesInGrid[0].size(),uniformGrid->trianglesInGrid[1].size()));

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

  cerr << "Pairs before unique: " << pairsTrianglesToProcess.size() << "\n";
  sort(pairsTrianglesToProcess.begin(),pairsTrianglesToProcess.end());
  vector<pair<InputTriangle *,InputTriangle *> >::iterator it = std::unique (pairsTrianglesToProcess.begin(), pairsTrianglesToProcess.end());
  pairsTrianglesToProcess.resize( std::distance(pairsTrianglesToProcess.begin(),it) ); // 10 20 30 20 10 

  cerr << "Pairs after unique: " << pairsTrianglesToProcess.size() << "\n";

 // clock_gettime(CLOCK_REALTIME, &t1);
  //cerr << "T after sort: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n";    
}

//returns the number of intersections found
//the sets are filled with the triangles (from the corresponding mesh) that intersect
unsigned long long  computeIntersections(const Nested3DGridWrapper *uniformGrid, unordered_set<const InputTriangle *> trianglesThatIntersect[2]) {
  timespec t0,t1;
  clock_gettime(CLOCK_REALTIME, &t0);
   

  unsigned long long totalIntersections = 0;

  vector<pair<InputTriangle *,InputTriangle *> > vtPairsTrianglesToProcess;
  
  cerr << "Getting pairs of triangles from grid cells" << "\n";
  getPairsTrianglesInSameUnifGridCells(uniformGrid,vtPairsTrianglesToProcess);
  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time creating list of pairs of triangles to process (intersection): " << convertTimeMsecs(diff(t0,t1))/1000 << "\n"; 
  //pairsTrianglesToProcess.reserve(1804900);
  
/*
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

            
  return totalIntersections;*/
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

  {
  MeshIntersectionGeometry meshIntersectionGeometry(mesh0Path,mesh1Path);

  clock_gettime(CLOCK_REALTIME, &t1);
  t0AfterDatasetRead = t1;
  cerr << "Time to read: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  timeReadData = convertTimeMsecs(diff(t0,t1))/1000; 
  Print_Current_Process_Memory_Used();


  
  meshIntersectionGeometry.printBoundingBoxes(); 

  
  
  cerr <<"Creating nested grid..." << endl;

  clock_gettime(CLOCK_REALTIME, &t0); 



//  void init(const vector<Triangle *> trianglesInsert[2], const vector<Point> vertices[2], const int gridSizeLevel1, const int gridSizeLevel2, const Point &p0, const Point &p1,const long long prodThreshold);

  TIME(uniformGrid.init(meshIntersectionGeometry, gridSizeLevel1,gridSizeLevel2,triggerSecondLevel));


  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to create and refine grid: " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  timeCreateGrid = convertTimeMsecs(diff(t0,t1))/1000; 

   //After the uniform grid is initialized, let's compute the intersection between the triangles...
  cerr << "Detecting intersections..." << endl;
  clock_gettime(CLOCK_REALTIME, &t0); 
  unordered_set<const InputTriangle *> trianglesThatIntersect[2];
  unsigned long long numIntersectionsDetected = computeIntersections(&uniformGrid,trianglesThatIntersect);

  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to detect intersections (includes time for computing statistics and for saving intersections for debugging): " << convertTimeMsecs(diff(t0,t1))/1000 << endl;
  Print_Current_Process_Memory_Used();






  }
  cerr << "Mesh intersection geometry freed" << endl; 



}