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
#include <deque>
#include <set>
#include <map>
#include <string>
#include <cstdio>
#include <cmath>

using namespace std;

#ifndef FLOOD_FILL_SCANLINE
#define FLOOD_FILL_SCANLINE

#include "rationals.h"
#include "utils.h"
#include "nested3DGrid.h"
#include "3dGeometry.h"

struct ScanLineInterval {
	int x,y;
	int z0,z1;
};

#include <stack>
#include <array>
//given a seed cell (gx,gy,gz), set the connected component (considering the 4-neighborhood) of (gx,gy,gz) to id//
void setObjectInWhereEmptyCellIs2(int gx,int gy,int gz,int nx,int ny, int nz, int gridSize,int nestedGridSize,GridCellsLabels &cellsLabels,ObjectId id,int chunkStartX,int chunkStartY,int chunkStartZ,int chunkEndX,int chunkEndY,int chunkEndZ);

#endif