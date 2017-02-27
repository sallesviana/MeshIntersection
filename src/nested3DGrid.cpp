#include <stdio.h>
#include <strings.h>
#include <climits>
#include <vector>
#include <array>

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



#include <unordered_map>
#include "nested3DGrid.h"

#include <list>


using namespace std;



void Nested3DGrid::deleteRaggedArrays() {
	for(int imap=0;imap<2;imap++) {
		delete []startPositionRaggedArray[imap];
		delete []triangles[imap];
	}
}

void Nested3DGrid::deleteMemory(int gridSizeLevel1) {
	deleteRaggedArrays();

	if (childGrids==NULL) return;

	for(int i=0;i<gridSizeLevel1;i++) {
  	for(int j=0;j<gridSizeLevel1;j++) {
  		for(int k=0;k<gridSizeLevel1;k++) {
  			if (childGrids[i][j][k] == NULL) continue;  			
  			childGrids[i][j][k]->deleteRaggedArrays();
  			delete childGrids[i][j][k];
  		}
  		delete []childGrids[i][j];
  	}
  	delete []childGrids[i];
	}
	delete []childGrids;

}






void Nested3DGridWrapper::initGridCells() {
  int xCellmin;
  int yCellmin;
  int zCellmin;

  int xCellmax;
  int yCellmax;
  int zCellmax; 
  
  int grid3 = gridSizeLevel1*gridSizeLevel1*gridSizeLevel1;
  int gridSize = gridSizeLevel1;
  vector<int> countTrianglesInEachGridTemp(grid3);

  MeshIntersectionGeometry &meshGeometry = *meshGeometryPtr;

  cerr << "Initializing grid cells first step..." << endl;

  grid.childGrids = new Nested3DGrid ***[gridSizeLevel1];
  for(int i=0;i<gridSizeLevel1;i++) {
  	grid.childGrids[i] = new Nested3DGrid **[gridSizeLevel1];
  	for(int j=0;j<gridSizeLevel1;j++) {
  		grid.childGrids[i][j] = new Nested3DGrid *[gridSizeLevel1];
  		for(int k=0;k<gridSizeLevel1;k++) {
  			grid.childGrids[i][j][k] = NULL;
  		}
  	}
  }


  timespec t0,t1;

  cerr << "Starting first pass..." << endl;


  for(int imap=0;imap<2;imap++) { // for each map...
    const size_t sz = meshGeometry.inputTriangles[imap].size();

    grid.startPositionRaggedArray[imap] = new InputTriangle**[grid3+1];
		for(int i=0;i<=grid3;i++) grid.startPositionRaggedArray[imap][i] = 0;

    for(size_t i=0;i<grid3;i++) countTrianglesInEachGridTemp[i] = 0;
    //we will performa first pass just to count the number of triangles we will insert in each grid cell...

    clock_gettime(CLOCK_REALTIME, &t0);
    for(size_t i = 0;i<sz;i++) { 
      InputTriangle*t = &meshGeometry.inputTriangles[imap][i];
      

      /*
      xCellmin =  gridCellEachPointLevel1[imap][ t->boundingBox[0][0] ][0]  ;
      yCellmin =  gridCellEachPointLevel1[imap][ t->boundingBox[0][1] ][1]  ;
      zCellmin =  gridCellEachPointLevel1[imap][ t->boundingBox[0][2] ][2]  ;

      xCellmax =  gridCellEachPointLevel1[imap][ t->boundingBox[1][0] ][0]  ;
      yCellmax =  gridCellEachPointLevel1[imap][ t->boundingBox[1][1] ][1]  ;
      zCellmax =  gridCellEachPointLevel1[imap][ t->boundingBox[1][2] ][2]  ;*/

      //getVertexIdInputTriangleWithMinCoord(int meshId, int triangleId, int coord)
      xCellmin =  gridCellEachPointLevel1[imap][ meshGeometry.getVertexIdInputTriangleWithMinCoord(imap, i, 0) ][0]  ;
      yCellmin =  gridCellEachPointLevel1[imap][ meshGeometry.getVertexIdInputTriangleWithMinCoord(imap, i, 1) ][1]  ;
      zCellmin =  gridCellEachPointLevel1[imap][ meshGeometry.getVertexIdInputTriangleWithMinCoord(imap, i, 2) ][2]  ;

      xCellmax =  gridCellEachPointLevel1[imap][ meshGeometry.getVertexIdInputTriangleWithMaxCoord(imap, i, 0) ][0]  ;
      yCellmax =  gridCellEachPointLevel1[imap][ meshGeometry.getVertexIdInputTriangleWithMaxCoord(imap, i, 1) ][1]  ;
      zCellmax =  gridCellEachPointLevel1[imap][ meshGeometry.getVertexIdInputTriangleWithMaxCoord(imap, i, 2) ][2]  ;

      

      assert(xCellmin>=0);
      assert(yCellmin>=0);
      assert(zCellmin>=0);
      assert(xCellmax<gridSizeLevel1);
      assert(yCellmax<gridSizeLevel1);
      assert(zCellmax<gridSizeLevel1);

		  for(int x = xCellmin;x<=xCellmax;x++) {
		  	const int xBase = x*gridSize*gridSize;
		    for(int y = yCellmin;y<=yCellmax;y++) {
		    	const int yBase = y*gridSize;
		      for(int z = zCellmin;z<=zCellmax;z++) {
      			countTrianglesInEachGridTemp[xBase+yBase+z]++;
		      }
		    }
		  }
		}
		clock_gettime(CLOCK_REALTIME, &t1);
  	cerr << "Time to perform first pass: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n";

		int totalNumberTriangles = 0;
		for(int i=0;i<grid3;i++) {
			totalNumberTriangles += countTrianglesInEachGridTemp[i];
		}
		grid.triangles[imap] = new InputTriangle*[totalNumberTriangles];


		
		//make the start pointers point to the correct positions....
		int acccumTotalTriangles =0;
		for(int i=0;i<=grid3;i++) {
			grid.startPositionRaggedArray[imap][i] = (grid.triangles[imap] + acccumTotalTriangles);
			if (i<grid3)		
				acccumTotalTriangles += countTrianglesInEachGridTemp[i];
		}



		clock_gettime(CLOCK_REALTIME, &t0);
		for(size_t i=0;i<grid3;i++) countTrianglesInEachGridTemp[i] = 0; //we will reuse this array to count the amount of triangles we've already inserted...
		for(size_t i = 0;i<sz;i++) { 
      InputTriangle*t = &(meshGeometry.inputTriangles[imap][i]);
      xCellmin =  gridCellEachPointLevel1[imap][ meshGeometry.getVertexIdInputTriangleWithMinCoord(imap, i, 0)  ][0]  ;
      yCellmin =  gridCellEachPointLevel1[imap][ meshGeometry.getVertexIdInputTriangleWithMinCoord(imap, i, 1) ][1]  ;
      zCellmin =  gridCellEachPointLevel1[imap][ meshGeometry.getVertexIdInputTriangleWithMinCoord(imap, i, 2) ][2]  ;

      xCellmax =  gridCellEachPointLevel1[imap][ meshGeometry.getVertexIdInputTriangleWithMaxCoord(imap, i, 0)  ][0]  ;
      yCellmax =  gridCellEachPointLevel1[imap][ meshGeometry.getVertexIdInputTriangleWithMaxCoord(imap, i, 1)  ][1]  ;
      zCellmax =  gridCellEachPointLevel1[imap][ meshGeometry.getVertexIdInputTriangleWithMaxCoord(imap, i, 2)  ][2]  ;

      assert(xCellmin>=0);
      assert(yCellmin>=0);
      assert(zCellmin>=0);
      assert(xCellmax<gridSizeLevel1);
      assert(yCellmax<gridSizeLevel1);
      assert(zCellmax<gridSizeLevel1);

		  for(int x = xCellmin;x<=xCellmax;x++) {
		  	const int xBase = x*gridSize*gridSize;
		    for(int y = yCellmin;y<=yCellmax;y++) {
		    	const int yBase = y*gridSize;
		      for(int z = zCellmin;z<=zCellmax;z++) {
		      	InputTriangle** startPositionCell = grid.startPositionRaggedArray[imap][xBase+yBase+z];
		      	*(startPositionCell+countTrianglesInEachGridTemp[xBase+yBase+z]) = t;
      			countTrianglesInEachGridTemp[xBase+yBase+z]++;
		      }
		    }
		  }
		}
		clock_gettime(CLOCK_REALTIME, &t1);
  	cerr << "Time to perform second pass: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n";
    
  }
}

void Nested3DGridWrapper::initNestedGridCells(InputTriangle** trianglesInsert[2],int numTrianglesInsert[2], Nested3DGrid &nestedGrid,int ig,int jg, int kg) {
  int xCellmin;
  int yCellmin;
  int zCellmin;

  int xCellmax;
  int yCellmax;
  int zCellmax; 

  int grid3 = gridSizeLevel2*gridSizeLevel2*gridSizeLevel2;
  int gridSize = gridSizeLevel2;
  vector<int> countTrianglesInEachGridTemp(grid3);



  //cerr << "Here..." << endl;
  nestedGrid.childGrids =  NULL; //we will use only up to two levels of the grid...

  const MeshIntersectionGeometry &meshGeometry = *meshGeometryPtr;

  //cerr << "Inserting triangles..." << endl;
  for(int imap=0;imap<2;imap++) { // for each map...
    const size_t sz = numTrianglesInsert[imap];

    //cerr << "Map, triangles: " << imap << " " << sz << endl;

    nestedGrid.startPositionRaggedArray[imap] = new InputTriangle**[grid3+1];
		for(int i=0;i<=grid3;i++) nestedGrid.startPositionRaggedArray[imap][i] = 0;

    for(size_t i=0;i<grid3;i++) countTrianglesInEachGridTemp[i] = 0;


   //cerr << "First pass.." << endl;
    //we will performa first pass just to count the number of triangles we will insert in each grid cell...
    for(size_t triangle = 0;triangle<sz;triangle++) { //for(const Chain &c : chains[imap]) {  // Iterate over edges in this map...
      InputTriangle*t = trianglesInsert[imap][triangle];

      int idInputTriangle = t- &meshGeometry.inputTriangles[imap][0];

      int vertexWithXMin = meshGeometry.getVertexIdInputTriangleWithMinCoord(imap,idInputTriangle,0);//t->boundingBox[0][0];
      int vertexWithYMin = meshGeometry.getVertexIdInputTriangleWithMinCoord(imap,idInputTriangle,1);//t->boundingBox[0][1];
      int vertexWithZMin = meshGeometry.getVertexIdInputTriangleWithMinCoord(imap,idInputTriangle,2);//t->boundingBox[0][2];

      int vertexWithXMax = meshGeometry.getVertexIdInputTriangleWithMaxCoord(imap,idInputTriangle,0);//t->boundingBox[1][0];
      int vertexWithYMax = meshGeometry.getVertexIdInputTriangleWithMaxCoord(imap,idInputTriangle,1);//t->boundingBox[1][1];
      int vertexWithZMax = meshGeometry.getVertexIdInputTriangleWithMaxCoord(imap,idInputTriangle,2);//t->boundingBox[1][2];      

      xCellmin =  gridCellEachPointLevel2[imap][ vertexWithXMin ][0]  ;
      
      if ( gridCellEachPointLevel1[imap][ vertexWithXMin ][0] != ig ) //if this point is outside the grid
        xCellmin = 0;   

      
      yCellmin =  gridCellEachPointLevel2[imap][ vertexWithYMin ][1]  ;
      if ( gridCellEachPointLevel1[imap][ vertexWithYMin ][1] != jg ) //if this point is outside the grid
        yCellmin = 0;
      
      zCellmin =  gridCellEachPointLevel2[imap][ vertexWithZMin ][2]  ;
      if ( gridCellEachPointLevel1[imap][ vertexWithZMin ][2] != kg ) //if this point is outside the grid
        zCellmin = 0;


      xCellmax =  gridCellEachPointLevel2[imap][ vertexWithXMax ][0]  ;
      if ( gridCellEachPointLevel1[imap][ vertexWithXMax ][0] != ig ) //if this point is outside the grid
        xCellmax = gridSizeLevel2-1;

      yCellmax =  gridCellEachPointLevel2[imap][ vertexWithYMax ][1]  ;
      if ( gridCellEachPointLevel1[imap][ vertexWithYMax ][1] != jg ) //if this point is outside the grid
        yCellmax = gridSizeLevel2-1;

      zCellmax =  gridCellEachPointLevel2[imap][ vertexWithZMax ][2]  ;
      if ( gridCellEachPointLevel1[imap][ vertexWithZMax ][2] != kg ) //if this point is outside the grid
        zCellmax = gridSizeLevel2-1;


      


      assert(xCellmin>=0);
      assert(yCellmin>=0);
      assert(zCellmin>=0);
      assert(xCellmax<gridSizeLevel2);
      assert(yCellmax<gridSizeLevel2);
      assert(zCellmax<gridSizeLevel2);
      //insertTriangleInGrid(gridCells,t,imap, xCellmin,yCellmin,zCellmin,xCellmax,yCellmax, zCellmax);

      for(int x = xCellmin;x<=xCellmax;x++) {
		  	const int xBase = x*gridSize*gridSize;
		    for(int y = yCellmin;y<=yCellmax;y++) {
		    	const int yBase = y*gridSize;
		      for(int z = zCellmin;z<=zCellmax;z++) {
      			countTrianglesInEachGridTemp[xBase+yBase+z]++;
		      }
		    }
		  }
    }

    //cerr << "Counting..." << endl;

    int totalNumberTriangles = 0;
		for(int i=0;i<grid3;i++) {
			totalNumberTriangles += countTrianglesInEachGridTemp[i];
		}
		nestedGrid.triangles[imap] = new InputTriangle*[totalNumberTriangles];

		
		//make the start pointers point to the correct positions....
		int acccumTotalTriangles =0;
		for(int i=0;i<=grid3;i++) {
			nestedGrid.startPositionRaggedArray[imap][i] = (nestedGrid.triangles[imap] + acccumTotalTriangles);
			if (i<grid3)		
				acccumTotalTriangles += countTrianglesInEachGridTemp[i];
		}

		for(size_t i=0;i<grid3;i++) countTrianglesInEachGridTemp[i] = 0; //we will reuse this array to count the amount of triangles we've already inserted...
		
		//cerr << "Writting to memory.. " << endl;
		for(size_t triangle = 0;triangle<sz;triangle++) { 
      InputTriangle*t = trianglesInsert[imap][triangle];

      int idInputTriangle = t- &meshGeometry.inputTriangles[imap][0];

      int vertexWithXMin = meshGeometry.getVertexIdInputTriangleWithMinCoord(imap,idInputTriangle,0);//t->boundingBox[0][0];
      int vertexWithYMin = meshGeometry.getVertexIdInputTriangleWithMinCoord(imap,idInputTriangle,1);//t->boundingBox[0][1];
      int vertexWithZMin = meshGeometry.getVertexIdInputTriangleWithMinCoord(imap,idInputTriangle,2);//t->boundingBox[0][2];

      int vertexWithXMax = meshGeometry.getVertexIdInputTriangleWithMaxCoord(imap,idInputTriangle,0);//t->boundingBox[1][0];
      int vertexWithYMax = meshGeometry.getVertexIdInputTriangleWithMaxCoord(imap,idInputTriangle,1);//t->boundingBox[1][1];
      int vertexWithZMax = meshGeometry.getVertexIdInputTriangleWithMaxCoord(imap,idInputTriangle,2);//t->boundingBox[1][2];      

      xCellmin =  gridCellEachPointLevel2[imap][ vertexWithXMin ][0]  ;
      if ( gridCellEachPointLevel1[imap][ vertexWithXMin ][0] != ig ) //if this point is outside the grid
        xCellmin = 0;
      
      yCellmin =  gridCellEachPointLevel2[imap][ vertexWithYMin ][1]  ;
      if ( gridCellEachPointLevel1[imap][ vertexWithYMin ][1] != jg ) //if this point is outside the grid
        yCellmin = 0;
      
      zCellmin =  gridCellEachPointLevel2[imap][ vertexWithZMin ][2]  ;
      if ( gridCellEachPointLevel1[imap][ vertexWithZMin ][2] != kg ) //if this point is outside the grid
        zCellmin = 0;


      xCellmax =  gridCellEachPointLevel2[imap][ vertexWithXMax ][0]  ;
      if ( gridCellEachPointLevel1[imap][ vertexWithXMax ][0] != ig ) //if this point is outside the grid
        xCellmax = gridSizeLevel2-1;

      yCellmax =  gridCellEachPointLevel2[imap][ vertexWithYMax ][1]  ;
      if ( gridCellEachPointLevel1[imap][ vertexWithYMax ][1] != jg ) //if this point is outside the grid
        yCellmax = gridSizeLevel2-1;

      zCellmax =  gridCellEachPointLevel2[imap][ vertexWithZMax ][2]  ;
      if ( gridCellEachPointLevel1[imap][ vertexWithZMax ][2] != kg ) //if this point is outside the grid
        zCellmax = gridSizeLevel2-1;

      assert(xCellmin>=0);
      assert(yCellmin>=0);
      assert(zCellmin>=0);
      assert(xCellmax<gridSizeLevel2);
      assert(yCellmax<gridSizeLevel2);
      assert(zCellmax<gridSizeLevel2);

		  for(int x = xCellmin;x<=xCellmax;x++) {
		  	const int xBase = x*gridSize*gridSize;
		    for(int y = yCellmin;y<=yCellmax;y++) {
		    	const int yBase = y*gridSize;
		      for(int z = zCellmin;z<=zCellmax;z++) {
		      	InputTriangle** startPositionCell = nestedGrid.startPositionRaggedArray[imap][xBase+yBase+z];
		      	*(startPositionCell+countTrianglesInEachGridTemp[xBase+yBase+z]) = t;
      			countTrianglesInEachGridTemp[xBase+yBase+z]++;
		      }
		    }
		  }
		}
  }


}


void Nested3DGridWrapper::refineChildGrid(int ig,int jg,int kg) {
	int sizeNestedGrids[2];
	InputTriangle**trianglesInsert[2];
	sizeNestedGrids[0] = grid.numTrianglesInGridCell(0, gridSizeLevel1,ig,jg,kg);
	sizeNestedGrids[1] = grid.numTrianglesInGridCell(1, gridSizeLevel1,ig,jg,kg);
	trianglesInsert[0] = grid.getPointerStartListTriangles(0,gridSizeLevel1,ig,jg,kg);
	trianglesInsert[1] = grid.getPointerStartListTriangles(1,gridSizeLevel1,ig,jg,kg);



	//cerr << "Refining child grid.. " << ig << " " << jg << " " << kg << endl;
	//cerr << "Triangles: " << sizeNestedGrids[0] << " " << sizeNestedGrids[1] << endl;
	grid.childGrids[ig][jg][kg] = new Nested3DGrid();
  initNestedGridCells(trianglesInsert,sizeNestedGrids, *grid.childGrids[ig][jg][kg] ,ig,jg,kg);
}




void Nested3DGridWrapper::init(MeshIntersectionGeometry &meshGeometry, const int gridSizeLevel1, const int gridSizeLevel2, const long long prodThreshold) {
  timespec t0,t1,t2;

  meshGeometryPtr = &meshGeometry;

  this->gridSizeLevel1 = gridSizeLevel1;
  this->gridSizeLevel2 = gridSizeLevel2;
  this->gridSize2Levels = gridSizeLevel1*gridSizeLevel2;


  array<VertCoord,3> coordRange = meshGeometry.coordRangeMeshes();


  // Prevent any points exactly at the right edge of the grid.
  VertCoord rangec = max(max(coordRange[0], coordRange[1]),coordRange[2])
                     * slightlyMoreThanOne;

  cellScale2Levels = gridSize2Levels/rangec;
  cellWidth2Levels = rangec/gridSize2Levels;

  cellScaleLevel1 = gridSizeLevel1/rangec;
  cellWidthLevel1 = rangec/gridSizeLevel1;

  cellScaleLevel2 =  cellScaleLevel1*gridSizeLevel2;
  cellWidthLevel2 = cellWidthLevel1/gridSizeLevel2;


  assert(cellScale2Levels>0);
  assert(cellWidth2Levels>0);

  assert(cellScaleLevel1>0);
  assert(cellWidthLevel1>0);

  assert(cellScaleLevel2>0);
  assert(cellWidthLevel2>0);

  //cerr << "Resizing..." << endl;
  this->trianglesInGrid[0].resize(meshGeometryPtr->inputTriangles[0].size());
  //cerr << "resized..." << endl;

  const int numTriMesh0 = trianglesInGrid[0].size();
  for(int i=0;i<numTriMesh0;i++) trianglesInGrid[0][i] = &(meshGeometryPtr->inputTriangles[0][i]);

  this->trianglesInGrid[1].resize(meshGeometryPtr->inputTriangles[1].size());
  const int numTriMesh1 = trianglesInGrid[1].size();
  for(int i=0;i<numTriMesh1;i++) trianglesInGrid[1][i] = &(meshGeometryPtr->inputTriangles[1][i]);


  //Computing in what grid cell each point is...
  clock_gettime(CLOCK_REALTIME, &t0);
  computeGridCellWhereEachPointIs();  
  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to compute in what grid cell each point is: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n";
  

 
  clock_gettime(CLOCK_REALTIME, &t0);
  initGridCells();
  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to insert triangles into grid: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n";


  //Refine child grids...

  clock_gettime(CLOCK_REALTIME, &t0);
  clock_gettime(CLOCK_REALTIME, &t1);

  vector<array<int,3> > gridCellsToRefine;
  for (int ig=0; ig<gridSizeLevel1; ig++) 
    for(int jg=0;jg<gridSizeLevel1;jg++) 
      for(int kg=0;kg<gridSizeLevel1;kg++) {
        int sizeGridCellImap0 = grid.numTrianglesInGridCell(0, gridSizeLevel1,ig,jg,kg);
        int sizeGridCellImap1 = grid.numTrianglesInGridCell(1, gridSizeLevel1,ig,jg,kg);

        if(sizeGridCellImap0*sizeGridCellImap1 >= prodThreshold) { //criteria to refine grid cells..
          //refine this grid cell..
        	//cerr << "Will refine: " << ig << " " << jg << " " << kg << " ---> " << sizeGridCellImap0 << " " << sizeGridCellImap1 << endl;
          gridCellsToRefine.push_back({ig,jg,kg});
        }
      }

  clock_gettime(CLOCK_REALTIME, &t2);
  cerr << "Time to determine what cells to refine: " << convertTimeMsecs(diff(t1,t2))/1000 << "\n";
        
  int numCellsToRefine = gridCellsToRefine.size();
  cerr << "Number of cells to refine: " << numCellsToRefine << "\n";




  clock_gettime(CLOCK_REALTIME, &t1);

  #pragma omp parallel for schedule(dynamic)
  for(int i = 0; i < numCellsToRefine; ++i) {
    refineChildGrid(gridCellsToRefine[i][0],gridCellsToRefine[i][1],gridCellsToRefine[i][2]);
  }


  clock_gettime(CLOCK_REALTIME, &t2);
  cerr << "Time to refine cells: " << convertTimeMsecs(diff(t1,t2))/1000 << "\n";
  cerr << "Time (total) to refine cells: " << convertTimeMsecs(diff(t0,t2))/1000 << "\n";


/*  long long ctLevel1 = 0;
  long long ctLevel2 = 0;
  long long triangles2hash = 0;
  for (int ig=0; ig<gridSizeLevel1; ig++) 
    for(int jg=0;jg<gridSizeLevel1;jg++) 
      for(int kg=0;kg<gridSizeLevel1;kg++) {
        Nested3DGridCell &nestedGrid = grid.gridCells[ig][jg][kg];
        ctLevel1 += nestedGrid.triangles[0].size();
        if(nestedGrid.childGrid) {
          Nested3DGrid &grid = *nestedGrid.childGrid;
          long long mul = ig*gridSizeLevel1*gridSizeLevel1 + jg*gridSizeLevel1 + kg;
          for (int ig=0; ig<gridSizeLevel2; ig++) 
            for(int jg=0;jg<gridSizeLevel2;jg++) 
              for(int kg=0;kg<gridSizeLevel2;kg++) {
                Nested3DGridCell &nestedGrid = grid.gridCells[ig][jg][kg];
                ctLevel2 += nestedGrid.triangles[0].size();
                long long mul2 = ig*gridSizeLevel2*gridSizeLevel2 + jg*gridSizeLevel2 + kg;
                triangles2hash += nestedGrid.triangles[0].size()*mul2*mul;
              }
        }        
      }
  cerr << "Triangles in level 1: " << ctLevel1 << endl;
  cerr << "Triangles in level 2: " << ctLevel2 << endl;
  cerr << "Triangles level 2 hash: " << triangles2hash << endl;*/

}





void Nested3DGridWrapper::computeGridCellWhereEachPointIs() {
  for(int meshId=0;meshId<2;meshId++) { 

    int numPointsImap = meshGeometryPtr->getNumVertices(meshId);
    gridCellEachPointLevel1[meshId].resize(numPointsImap);
    gridCellEachPointLevel2[meshId].resize(numPointsImap);
  
    #pragma omp parallel
    {

        MeshIntersectionGeometry::TempVarsGetGridCellContainingVertex tempVars;

        int threadId = omp_get_thread_num();

          #pragma omp for schedule(dynamic,1000)
          for(int i=0;i<numPointsImap;i++) {
            const array<int,3> gridCell = meshGeometryPtr->getGridCellContainingVertex(meshId,i,cellScale2Levels,tempVars);

            gridCellEachPointLevel1[meshId][i][0] = gridCell[0]/gridSizeLevel2;
            gridCellEachPointLevel1[meshId][i][1] = gridCell[1]/gridSizeLevel2;
            gridCellEachPointLevel1[meshId][i][2] = gridCell[2]/gridSizeLevel2;

            gridCellEachPointLevel2[meshId][i][0] = gridCell[0]%gridSizeLevel2;
            gridCellEachPointLevel2[meshId][i][1] = gridCell[1]%gridSizeLevel2;
            gridCellEachPointLevel2[meshId][i][2] = gridCell[2]%gridSizeLevel2;

            
          }
    }    
  }
}