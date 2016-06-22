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



#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>
#include "common2.h"
#include <list>


using namespace std;



//Struct used (temporaly) to store labels of empty grid cells (that are completely inside objects)
//we could have added this information in the data structure below (and, thus, do not create another dedicated data structure) but
//we decided to create a separeta data structure in order to avoid adding more information to the uniform grid data structure
//(also, we do not need this information in all parts of our program...)
struct GridCellsLabels {
  vector<vector<vector<int > > > labels;
  vector<vector<vector<GridCellsLabels * > > > childGridLabels;

  bool has2ndLevel(int x,int y,int z) {
    return childGridLabels[x][y][z]!=NULL;
  }

  GridCellsLabels(int resolution) {
    childGridLabels = vector<vector<vector<GridCellsLabels * > > >(resolution,vector<vector<GridCellsLabels * > >(resolution, vector<GridCellsLabels *>(resolution,NULL)));
    labels = vector<vector<vector<int > > >(resolution,vector<vector<int > >(resolution, vector<int>(resolution,DONT_KNOW_ID)));
  }
  ~GridCellsLabels() {
    int gridSize = childGridLabels.size();
    for(int gx=0;gx<gridSize;gx++) 
      for(int gy=0;gy<gridSize;gy++)
        for(int gz=0;gz<gridSize;gz++)
          delete childGridLabels[gx][gy][gz];
  }
};




/********************************************************************************************************************
********************************************************************************************************************
********************************************************************************************************************
********************************************************************************************************************
********************************************************************************************************************/

// TODO: maybe use contin. fract. to simplify the quadtree coordinates
// TODO: use temporary variables allocated in the beginning of the problem...
//We consider that points in the edges of the quadtree belong to both quadrants (this is not a problem... some intersections may be computed twice)
//Maybe change cell_scale to allow multiplication instead of division...
//Maybe try to save memory with the unif grid boundin boxes that are stored in the grid...

//int rat2int2(const VertCoord &r) {
  //return  castIntWrap(big_int((numeratorWrap(r)/denominatorWrap(r))));
  //return (int)r.get_d();
//}

class Nested3DGrid;


class Nested3DGrid {
public:
	Triangle * *triangles[2]; //ragged array of triangles...

	//startPositionRaggedArray[i] points to the position of the triangles arrray where the linearized cell i starts (if i = gridSize^3 --> this position is one position after the end of the array)
	Triangle * **startPositionRaggedArray[2]; 

  Nested3DGrid ****childGrids;
  
  bool hasSecondLevel(int ig, int jg, int kg) const {
    return childGrids[ig][jg][kg] != NULL;
  }

  Nested3DGrid * getChildGrid(int ig, int jg, int kg) const {
  	return childGrids[ig][jg][kg];
  }

  Triangle ** getPointerStartListTriangles(const int imap,const int gridSize,int x,int y, int z) const {
		return startPositionRaggedArray[imap][gridSize*gridSize*x + y*gridSize + z];
	}

	int numTrianglesInGridCell(int imap, int gridSize,int x,int y, int z) const {
		return startPositionRaggedArray[imap][gridSize*gridSize*x + y*gridSize + z+1] - startPositionRaggedArray[imap][gridSize*gridSize*x + y*gridSize + z];
	}

	int numTriangles(int imap, int gridSize) const { //we need the grid size because we don't store the grid size in each nested grid... (to save memory)
		return startPositionRaggedArray[imap][gridSize*gridSize*gridSize] - startPositionRaggedArray[imap][0];
	}


	void deleteRaggedArrays();
  void deleteMemory(int gridSizeLevel1); //we do not have a "regular" destructor because we want to save memory by not storing the size of the grid...
};

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




struct Nested3DGridWrapper {
  void init(vector<Triangle *> trianglesInsert[2], const vector<Point> vertices[2], const int gridSizeLevel1, const int gridSizeLevel2, const Point &p0, const Point &p1,const long long prodThreshold);
  int gridSizeLevel1,gridSizeLevel2,gridSize2Levels; //gridSize2Levels is the product of the two grid sizes...
  Point box[2]; //3d box representing this uniform grid...
  VertCoord cellWidthLevel1,cellWidthLevel2,cellWidth2Levels;
  VertCoord cellScaleLevel1,cellScaleLevel2,cellScale2Levels;

  vector< array<int,3> > gridCellEachPointLevel1[2],gridCellEachPointLevel2[2];
  const vector<Point> *vertices[2];
  vector<Triangle *> *trianglesInGrid[2];

  void initGridCells(vector<Triangle *> trianglesInsert[2]);
  void initGridCells();
  Nested3DGrid grid;
  

  void computeGridCellWhereEachPointIs();

  void refineChildGrid(int ig,int jg,int kg);
	void initNestedGridCells(Triangle ** trianglesInsert[2],int numTrianglesInsert[2], Nested3DGrid &nestedGrid,int ig,int jg, int kg);

	~Nested3DGridWrapper() {
		grid.deleteMemory(gridSizeLevel1);
	}
  /*Triangle ** getPointerStartListTriangles(int ig, int jg, int kg,int idTriangle);
  int getNumTrianglesFirstLevelCell(int ig, int jg, int kg);
  Triangle *&getTriangleSecondLevelCell(Nested3DGrid &firstLevelGrid,int ig, int jg, int kg,int idTriangle);
  int getNumTriangleSecondLevelCell(Nested3DGrid &firstLevelGrid,int ig, int jg, int kg);
	*/

  


  //Returns the "global" cell where a point with a given coordinate is...
  //We call a global grid coordinate the coordinate supposing the grid has  resolution (gridSizeLevel1*gridSizeLevel2)^2
  //Thus, if the x global coordinate ix xg --> the coordinate of the point in the first level will be (xg/gridSizeLevel2) and in the second level it will be (xg%gridSizeLevel2)
  int x_global_cell_from_coord(const VertCoord &x, VertCoord &tempVar,big_int tempVarsInt[]) const {
    tempVar = x;
    tempVar -= box[0][0];
    tempVar *= cellScale2Levels;
    const int c  = convertToInt(tempVar,tempVarsInt);
    //const int c = rat2int2((x-box[0][0])*cellScale);
    //Removed this assertion because this will happen legally in the nested grids...
    //ASSERTT(c>=0&&c<gridSize,"Illegal c computed in x_cell_from_coord",cerr << PRINTC(x) << PRINTC(c)<< PRINTC(box) << PRINTN(cellScale));
    return c;
  }
  int x_cell_from_coord_level1(const VertCoord &x, VertCoord &tempVar,big_int tempVarsInt[]) const {
    return x_global_cell_from_coord(x, tempVar,tempVarsInt)/gridSizeLevel2;
  }
  int x_cell_from_coord_level2(const VertCoord &x, VertCoord &tempVar,big_int tempVarsInt[]) const {
    return x_global_cell_from_coord(x, tempVar,tempVarsInt)%gridSizeLevel2;
  }



  int y_global_cell_from_coord(const VertCoord &y, VertCoord &tempVar,big_int tempVarsInt[]) const{
    tempVar = y;
    tempVar -= box[0][1];
    tempVar *= cellScale2Levels;
    const int c  = convertToInt(tempVar,tempVarsInt);

    //const int c = rat2int2((y-box[0][1])*cellScale);
    // This can happen legally.
    //  ASSERTT(c>=0&&c<gridSize,"Illegal c computed in y_cell_from_coord",cerr << PRINTC(y) << PRINTC(c)<< PRINTC(box) << PRINTN(cellScale));
    return c;
  }
  int y_cell_from_coord_level1(const VertCoord &y, VertCoord &tempVar,big_int tempVarsInt[]) const {
    return y_global_cell_from_coord(y, tempVar,tempVarsInt)/gridSizeLevel2;
  }
  int y_cell_from_coord_level2(const VertCoord &y, VertCoord &tempVar,big_int tempVarsInt[]) const {
    return y_global_cell_from_coord(y, tempVar,tempVarsInt)%gridSizeLevel2;
  }

  int z_global_cell_from_coord(const VertCoord &z, VertCoord &tempVar,big_int tempVarsInt[]) const{
    tempVar = z;
    tempVar -= box[0][2];
    tempVar *= cellScale2Levels;
    const int c  = convertToInt(tempVar,tempVarsInt);

    //const int c = rat2int2((y-box[0][1])*cellScale);
    // This can happen legally.
    //  ASSERTT(c>=0&&c<gridSize,"Illegal c computed in y_cell_from_coord",cerr << PRINTC(y) << PRINTC(c)<< PRINTC(box) << PRINTN(cellScale));
    return c;
  }
  int z_cell_from_coord_level1(const VertCoord &z, VertCoord &tempVar,big_int tempVarsInt[]) const {
    return z_global_cell_from_coord(z, tempVar,tempVarsInt)/gridSizeLevel2;
  }
  int z_cell_from_coord_level2(const VertCoord &z, VertCoord &tempVar,big_int tempVarsInt[]) const {
    return z_global_cell_from_coord(z, tempVar,tempVarsInt)%gridSizeLevel2;
  }

   /*

  void insertTriangleInGrid(Triangle *t,const vector<Point> vertices[2],int imap,int xCellmin,int yCellmin,int zCellmin,int xCellmax,int yCellmax, int zCellmax);



	void initSeq(const vector<Triangle *> trianglesInsert[2], const vector<Point> vertices[2], const int gridSize, const Point &p0, const Point &p1, bool runInParallel);
	
	

	void make_grid();

	

  void printStats() {
    cerr << "Bounding box: " << box[0][0] << " " << box[0][1] << " -- " << box[1][0] << " " << box[1][1] << endl;
    cerr << "Grid size level 1: " << gridSizeLevel1 << endl;
    cerr << "Cell width level 1: " << cellWidthLevel1 << endl;
    cerr << "Cell scale level 1: " << cellScaleLevel1 << endl;
  }

  
	

	// Find the cell number from the internal coordinates.
	// A cell contains its left and lower edges, but not its right and upper edges.
	CellNo cell_from_coord(const Point &p, VertCoord &tempVar,big_int tempVarsInt[]) const{
    CellNo cell;
	  cell[0] = x_cell_from_coord(p[0],tempVar,tempVarsInt);
	  cell[1] = y_cell_from_coord(p[1],tempVar,tempVarsInt);
    cell[2] = z_cell_from_coord(p[2],tempVar,tempVarsInt);
	  return cell;
	}

	VertCoord x_coord_from_cell(const int c) const{   // return the coord of the left edge.
	  return c/cellScale+box[0][0];
	}

	VertCoord y_coord_from_cell(const int c) const{
	  return c/cellScale+box[0][1];
	}

  VertCoord z_coord_from_cell(const int c) const{
	  return c/cellScale+box[0][2];
	}

	VertCoord x_coord_from_cell(const int c, VertCoord &tempVar) const{   // return the coord of the left edge.
	  tempVar = c;
	  tempVar /= cellScale;
	  tempVar += box[0][0];
	  return tempVar;
	}

	VertCoord y_coord_from_cell(const int c, VertCoord &tempVar) const{
	  tempVar = c;
	  tempVar /= cellScale;
	  tempVar += box[0][1];
	  return tempVar;
	}

  VertCoord z_coord_from_cell(const int c, VertCoord &tempVar) const{
	  tempVar = c;
	  tempVar /= cellScale;
	  tempVar += box[0][2];
	  return tempVar;
	}*/

};









void Nested3DGridWrapper::initGridCells(vector<Triangle *> trianglesInsert[2]) {
  int xCellmin;
  int yCellmin;
  int zCellmin;

  int xCellmax;
  int yCellmax;
  int zCellmax; 
  
  int grid3 = gridSizeLevel1*gridSizeLevel1*gridSizeLevel1;
  int gridSize = gridSizeLevel1;
  vector<int> countTrianglesInEachGridTemp(grid3);

  
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

  for(int imap=0;imap<2;imap++) { // for each map...
    const size_t sz = trianglesInsert[imap].size();

    grid.startPositionRaggedArray[imap] = new Triangle **[grid3+1];
		for(int i=0;i<=grid3;i++) grid.startPositionRaggedArray[imap][i] = 0;

    for(size_t i=0;i<grid3;i++) countTrianglesInEachGridTemp[i] = 0;
    //we will performa first pass just to count the number of triangles we will insert in each grid cell...

    clock_gettime(CLOCK_REALTIME, &t0);
    for(size_t i = 0;i<sz;i++) { 
      Triangle *t = trianglesInsert[imap][i];
      xCellmin =  gridCellEachPointLevel1[imap][ t->boundingBox[0][0] ][0]  ;
      yCellmin =  gridCellEachPointLevel1[imap][ t->boundingBox[0][1] ][1]  ;
      zCellmin =  gridCellEachPointLevel1[imap][ t->boundingBox[0][2] ][2]  ;

      xCellmax =  gridCellEachPointLevel1[imap][ t->boundingBox[1][0] ][0]  ;
      yCellmax =  gridCellEachPointLevel1[imap][ t->boundingBox[1][1] ][1]  ;
      zCellmax =  gridCellEachPointLevel1[imap][ t->boundingBox[1][2] ][2]  ;

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
		grid.triangles[imap] = new Triangle *[totalNumberTriangles];


		
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
      Triangle *t = trianglesInsert[imap][i];
      xCellmin =  gridCellEachPointLevel1[imap][ t->boundingBox[0][0] ][0]  ;
      yCellmin =  gridCellEachPointLevel1[imap][ t->boundingBox[0][1] ][1]  ;
      zCellmin =  gridCellEachPointLevel1[imap][ t->boundingBox[0][2] ][2]  ;

      xCellmax =  gridCellEachPointLevel1[imap][ t->boundingBox[1][0] ][0]  ;
      yCellmax =  gridCellEachPointLevel1[imap][ t->boundingBox[1][1] ][1]  ;
      zCellmax =  gridCellEachPointLevel1[imap][ t->boundingBox[1][2] ][2]  ;

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
		      	Triangle ** startPositionCell = grid.startPositionRaggedArray[imap][xBase+yBase+z];
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

void Nested3DGridWrapper::initNestedGridCells(Triangle ** trianglesInsert[2],int numTrianglesInsert[2], Nested3DGrid &nestedGrid,int ig,int jg, int kg) {
  int xCellmin;
  int yCellmin;
  int zCellmin;

  int xCellmax;
  int yCellmax;
  int zCellmax; 

  //cerr << "initing.." << endl;
  int grid3 = gridSizeLevel2*gridSizeLevel2*gridSizeLevel2;
  int gridSize = gridSizeLevel2;
  vector<int> countTrianglesInEachGridTemp(grid3);

  //cerr << "Here..." << endl;
  nestedGrid.childGrids =  NULL; //we will use only up to two levels of the grid...

  //cerr << "initing.." << endl;

  //cerr << "Inserting triangles..." << endl;
  for(int imap=0;imap<2;imap++) { // for each map...
    const size_t sz = numTrianglesInsert[imap];

    //cerr << "Map, triangles: " << imap << " " << sz << endl;

    nestedGrid.startPositionRaggedArray[imap] = new Triangle **[grid3+1];
		for(int i=0;i<=grid3;i++) nestedGrid.startPositionRaggedArray[imap][i] = 0;

    for(size_t i=0;i<grid3;i++) countTrianglesInEachGridTemp[i] = 0;


   //cerr << "First pass.." << endl;
    //we will performa first pass just to count the number of triangles we will insert in each grid cell...
    for(size_t i = 0;i<sz;i++) { //for(const Chain &c : chains[imap]) {  // Iterate over edges in this map...
      Triangle *t = trianglesInsert[imap][i];
      xCellmin =  gridCellEachPointLevel2[imap][ t->boundingBox[0][0] ][0]  ;
      if ( gridCellEachPointLevel1[imap][ t->boundingBox[0][0] ][0] != ig ) //if this point is outside the grid
        xCellmin = 0;

      yCellmin =  gridCellEachPointLevel2[imap][ t->boundingBox[0][1] ][1]  ;
      if ( gridCellEachPointLevel1[imap][ t->boundingBox[0][1] ][1] != jg ) //if this point is outside the grid
        yCellmin = 0;

      zCellmin =  gridCellEachPointLevel2[imap][ t->boundingBox[0][2] ][2]  ;
      if ( gridCellEachPointLevel1[imap][ t->boundingBox[0][2] ][2] != kg ) //if this point is outside the grid
        zCellmin = 0;


      xCellmax =  gridCellEachPointLevel2[imap][ t->boundingBox[1][0] ][0]  ;
      if ( gridCellEachPointLevel1[imap][ t->boundingBox[1][0] ][0] != ig ) //if this point is outside the grid
        xCellmax = gridSizeLevel2-1;

      yCellmax =  gridCellEachPointLevel2[imap][ t->boundingBox[1][1] ][1]  ;
      if ( gridCellEachPointLevel1[imap][ t->boundingBox[1][1] ][1] != jg ) //if this point is outside the grid
        yCellmax = gridSizeLevel2-1;

      zCellmax =  gridCellEachPointLevel2[imap][ t->boundingBox[1][2] ][2]  ;
      if ( gridCellEachPointLevel1[imap][ t->boundingBox[1][2] ][2] != kg ) //if this point is outside the grid
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
		nestedGrid.triangles[imap] = new Triangle *[totalNumberTriangles];
		
		//make the start pointers point to the correct positions....
		int acccumTotalTriangles =0;
		for(int i=0;i<=grid3;i++) {
			nestedGrid.startPositionRaggedArray[imap][i] = (nestedGrid.triangles[imap] + acccumTotalTriangles);
			if (i<grid3)		
				acccumTotalTriangles += countTrianglesInEachGridTemp[i];
		}

		for(size_t i=0;i<grid3;i++) countTrianglesInEachGridTemp[i] = 0; //we will reuse this array to count the amount of triangles we've already inserted...
		
		//cerr << "Writting to memory.. " << endl;
		for(size_t i = 0;i<sz;i++) { 
      Triangle *t = trianglesInsert[imap][i];
      xCellmin =  gridCellEachPointLevel2[imap][ t->boundingBox[0][0] ][0]  ;
      if ( gridCellEachPointLevel1[imap][ t->boundingBox[0][0] ][0] != ig ) //if this point is outside the grid
        xCellmin = 0;

      yCellmin =  gridCellEachPointLevel2[imap][ t->boundingBox[0][1] ][1]  ;
      if ( gridCellEachPointLevel1[imap][ t->boundingBox[0][1] ][1] != jg ) //if this point is outside the grid
        yCellmin = 0;

      zCellmin =  gridCellEachPointLevel2[imap][ t->boundingBox[0][2] ][2]  ;
      if ( gridCellEachPointLevel1[imap][ t->boundingBox[0][2] ][2] != kg ) //if this point is outside the grid
        zCellmin = 0;


      xCellmax =  gridCellEachPointLevel2[imap][ t->boundingBox[1][0] ][0]  ;
      if ( gridCellEachPointLevel1[imap][ t->boundingBox[1][0] ][0] != ig ) //if this point is outside the grid
        xCellmax = gridSizeLevel2-1;

      yCellmax =  gridCellEachPointLevel2[imap][ t->boundingBox[1][1] ][1]  ;
      if ( gridCellEachPointLevel1[imap][ t->boundingBox[1][1] ][1] != jg ) //if this point is outside the grid
        yCellmax = gridSizeLevel2-1;

      zCellmax =  gridCellEachPointLevel2[imap][ t->boundingBox[1][2] ][2]  ;
      if ( gridCellEachPointLevel1[imap][ t->boundingBox[1][2] ][2] != kg ) //if this point is outside the grid
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
		      	Triangle ** startPositionCell = nestedGrid.startPositionRaggedArray[imap][xBase+yBase+z];
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
	Triangle **trianglesInsert[2];
	sizeNestedGrids[0] = grid.numTrianglesInGridCell(0, gridSizeLevel1,ig,jg,kg);
	sizeNestedGrids[1] = grid.numTrianglesInGridCell(1, gridSizeLevel1,ig,jg,kg);
	trianglesInsert[0] = grid.getPointerStartListTriangles(0,gridSizeLevel1,ig,jg,kg);
	trianglesInsert[1] = grid.getPointerStartListTriangles(1,gridSizeLevel1,ig,jg,kg);


	//cerr << "Refining child grid.. " << ig << " " << jg << " " << kg << endl;
	//cerr << "Triangles: " << sizeNestedGrids[0] << " " << sizeNestedGrids[1] << endl;
	grid.childGrids[ig][jg][kg] = new Nested3DGrid();
  initNestedGridCells(trianglesInsert,sizeNestedGrids, *grid.childGrids[ig][jg][kg] ,ig,jg,kg);
}
/*
void Nested3DGrid::refineChildGrids(const vector<Point> vertices[2], const long long prodThreshold, int sizeNested3DGrid,int maxRefineDepth) {
  if(maxRefineDepth<=0) return; //we will not refine further if the limit has veen achieved

	vector< array<int,3> > childToRefine;


  for(int xcell=0;xcell<gridSize;xcell++)
    for(int ycell=0;ycell<gridSize;ycell++)
      for(int zcell=0;zcell<gridSize;zcell++) {
        Nested3DGridCell &c = gridCells[xcell][ycell][zcell];
  			//sumSizes += c.e[0].size();
  			if (c.triangles[0].size() >= prodThreshold || c.triangles[1].size() >= prodThreshold) {
          array<int,3> cell{{xcell, ycell, zcell}};
  				childToRefine.push_back( cell );
  			}
      }



	//cerr << "Size child grid cells to refine: " << childToRefine.size() << endl;

	const int sz = childToRefine.size();

  cerr << "Number of cells to refine: " << sz << endl;
  
  
  #pragma omp parallel
  {
    Point p0;

    #pragma omp for schedule(dynamic,1)
  	for(int i=0;i<sz;i++) {
  		int xg = childToRefine[i][0];
  		int yg = childToRefine[i][1];
      int zg = childToRefine[i][2];

  		//cerr << i << endl;

  		Nested3DGridCell &c = gridCells[xg][yg][zg];
  		c.childGrid = new Nested3DGrid();
  		p0 = box[0];
  		p0[0] += cellWidth*xg;
  		p0[1] += cellWidth*yg;
      p0[2] += cellWidth*zg;
  		Point p1 = p0;
  		p1[0] += cellWidth;
  		p1[1] += cellWidth;
      p1[2] += cellWidth;
  		c.childGrid->initSeq(c.triangles, vertices,  sizeNested3DGrid,p0,p1,false );
      c.childGrid->refineChildGrids(vertices, prodThreshold, sizeNested3DGrid,maxRefineDepth-1);
  	}

  }
}



void Nested3DGrid::make_grid() {
  ASSERTC(gridSize>0, "Grid size must be positive.");

  // Prevent any points exactly at the right edge of the grid.
  VertCoord rangec = max(max(box[1][0]-box[0][0], box[1][1]-box[0][1]),box[1][2]-box[0][2])
                     * slightlyMoreThanOne;

  cellScale =  gridSize/rangec;
  cellWidth = rangec/gridSize;

  assert(cellScale>0);
  assert(cellWidth>0);

  gridCells.resize(gridSize);



  for (int ig=0; ig<gridSize; ig++) {
     gridCells[ig].resize(gridSize);
     for(int jg=0;jg<gridSize;jg++)
      gridCells[ig][jg].resize(gridSize);
  }

}






*/















void Nested3DGridWrapper::computeGridCellWhereEachPointIs() {
  for(int imap=0;imap<2;imap++) {
    const vector<Point> &points = *vertices[imap];

    int numPointsImap = points.size();
    gridCellEachPointLevel1[imap].resize(numPointsImap);
    gridCellEachPointLevel2[imap].resize(numPointsImap);

  
    #pragma omp parallel
    {

        VertCoord tempVar;
        big_int tempVarsInt[3];

          #pragma omp for schedule(dynamic,1000)
          for(int i=0;i<numPointsImap;i++) {
            tempVar = points[i][0];
            tempVar -= box[0][0];
            tempVar *= cellScale2Levels;
            const int x = convertToInt(tempVar,tempVarsInt);

            tempVar = points[i][1];
            tempVar -= box[0][1];
            tempVar *= cellScale2Levels;
            const int y = convertToInt(tempVar,tempVarsInt);

            tempVar = points[i][2];
            tempVar -= box[0][2];
            tempVar *= cellScale2Levels;
            const int z  = convertToInt(tempVar,tempVarsInt);

            gridCellEachPointLevel1[imap][i][0] = x/gridSizeLevel2;
            gridCellEachPointLevel1[imap][i][1] = y/gridSizeLevel2;
            gridCellEachPointLevel1[imap][i][2] = z/gridSizeLevel2;

            gridCellEachPointLevel2[imap][i][0] = x%gridSizeLevel2;
            gridCellEachPointLevel2[imap][i][1] = y%gridSizeLevel2;
            gridCellEachPointLevel2[imap][i][2] = z%gridSizeLevel2;
          }
    }    
  }
}

void Nested3DGridWrapper::init(vector<Triangle *> trianglesInsert[2], const vector<Point> vertices[2], const int gridSizeLevel1, const int gridSizeLevel2, const Point &p0, const Point &p1,const long long prodThreshold) {
  timespec t0,t1,t2;

  this->gridSizeLevel1 = gridSizeLevel1;
  this->gridSizeLevel2 = gridSizeLevel2;
  this->gridSize2Levels = gridSizeLevel1*gridSizeLevel2;

  box[0] = p0;
  box[1] = p1;

    // Prevent any points exactly at the right edge of the grid.
  VertCoord rangec = max(max(box[1][0]-box[0][0], box[1][1]-box[0][1]),box[1][2]-box[0][2])
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


  this->vertices[0] = &(vertices[0]);
  this->vertices[1] = &(vertices[1]);
  this->trianglesInGrid[0] = &(trianglesInsert[0]);
  this->trianglesInGrid[1] = &(trianglesInsert[1]);

  //Computing in what grid cell each point is...
  clock_gettime(CLOCK_REALTIME, &t0);
  computeGridCellWhereEachPointIs();  
  clock_gettime(CLOCK_REALTIME, &t1);
  cerr << "Time to compute in what grid cell each point is: " << convertTimeMsecs(diff(t0,t1))/1000 << "\n";
  

 
  clock_gettime(CLOCK_REALTIME, &t0);
  initGridCells(trianglesInsert);
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

        if(sizeGridCellImap0 >= prodThreshold || sizeGridCellImap1 >= prodThreshold) { //criteria to refine grid cells..
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