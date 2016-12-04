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

#ifndef NESTED_3D_GRID
#define NESTED_3D_GRID

#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <cstdlib>
#include <unordered_set>
#include <unordered_map>
#include "rationals.h"
#include "3d_objects.h"
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

  int get_global_x_coord_mesh_vertex(int iMesh,int vertexId) const {
    return gridCellEachPointLevel1[iMesh][vertexId][0]*gridSizeLevel2 + gridCellEachPointLevel2[iMesh][vertexId][0];
  }
  int get_global_y_coord_mesh_vertex(int iMesh,int vertexId) const {
    return gridCellEachPointLevel1[iMesh][vertexId][1]*gridSizeLevel2 + gridCellEachPointLevel2[iMesh][vertexId][1];
  }
  int get_global_z_coord_mesh_vertex(int iMesh,int vertexId) const {
    return gridCellEachPointLevel1[iMesh][vertexId][2]*gridSizeLevel2 + gridCellEachPointLevel2[iMesh][vertexId][2];
  }

  //Returns the "global" cell where a point with a given coordinate is...
  //We call a global grid coordinate the coordinate supposing the grid has  resolution (gridSizeLevel1*gridSizeLevel2)^3
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

   
};







#endif