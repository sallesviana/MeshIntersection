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

#include "floodFillScanline.h"



#include <stack>
#include <array>
//given a seed cell (gx,gy,gz), set the connected component (considering the 4-neighborhood) of (gx,gy,gz) to id//
void setObjectInWhereEmptyCellIs2(int gx,int gy,int gz,int nx,int ny, int nz, int gridSize,int nestedGridSize,GridCellsLabels &cellsLabels,ObjectId id,int chunkStartX,int chunkStartY,int chunkStartZ,int chunkEndX,int chunkEndY,int chunkEndZ) {
	stack< array<int,3> > seedsLevel1;
	stack< pair<array<int,3>,array<int,3> > > seedsLevel2;

	std::vector<std::vector<std::vector<int> > > &idsFirstLevel = cellsLabels.labels;
	int oldId;
	if (nx<0) { //the cell is in the first level of the grid...
	  oldId = idsFirstLevel[gx][gy][gz];
	  seedsLevel1.push({gx,gy,gz});
	} else {
		oldId = cellsLabels.childGridLabels[gx][gy][gz]->labels[nx][ny][nz];
		seedsLevel2.push(pair<array<int,3>,array<int,3> >({gx,gy,gz},{nx,ny,nz}));
	}


	while(!seedsLevel1.empty() || !seedsLevel2.empty() ) {
		if(!seedsLevel1.empty()) {
			//cerr << "Getting top first level..." << endl;
			array<int,3> seed = seedsLevel1.top();
			seedsLevel1.pop();

			int x = seed[0];
			int y = seed[1];
			int z = seed[2];

			//cerr << "Processing: " << x <<  " " << y << " " << z << endl;
			if(idsFirstLevel[x][y][z] != oldId) continue;

			//assert(!cellsLabels.has2ndLevel(x,y,z));
			int i;
			for(i=z;i<chunkEndZ && idsFirstLevel[x][y][i]==oldId && !cellsLabels.has2ndLevel(x,y,i);i++) {
				idsFirstLevel[x][y][i] = id;
			}
			// position (i-1) is the last one that will be filled
			int endZ = i-1;
			for(i=z-1;i>=chunkStartZ && idsFirstLevel[x][y][i]==oldId  && !cellsLabels.has2ndLevel(x,y,i);i--) {
				idsFirstLevel[x][y][i] = id;
			}
			// (i+1) is the last z that will be filled (in reverse order..)
			int startZ = i+1;
			//so, we filled the interval [startZ, endZ] with id...
			//lets seed the neighbors...
			if (startZ!=chunkStartZ) { 
				//we have something before the start...  if we have a nested grid before the start, we need to add to the seeds the cells with id equal to the old one...
				//if we don't have a nested grid, this means that we stopped because we reached a cell that will not be labeled (so, we don't need to do anything.)
				if (cellsLabels.has2ndLevel(x,y,startZ-1)) {
					std::vector<std::vector<std::vector<int> > > &idsSecondLevel = cellsLabels.childGridLabels[x][y][startZ-1]->labels;
					for(int nx=0;nx<nestedGridSize;nx++)
						for(int ny=0;ny<nestedGridSize;ny++) {
							if(idsSecondLevel[nx][ny][nestedGridSize-1]==oldId) {
								seedsLevel2.push(  pair<array<int,3>,array<int,3> >({x,y,startZ-1},{nx,ny,nestedGridSize-1})  );
							}
						}
				}
			}
			if(endZ!=chunkEndZ-1) {
				//idem for the end...
				if (cellsLabels.has2ndLevel(x,y,endZ+1)) {
					std::vector<std::vector<std::vector<int> > > &idsSecondLevel = cellsLabels.childGridLabels[x][y][endZ+1]->labels;
					for(int nx=0;nx<nestedGridSize;nx++)
						for(int ny=0;ny<nestedGridSize;ny++) {
							if(idsSecondLevel[nx][ny][0]==oldId) {
								seedsLevel2.push(  pair<array<int,3>,array<int,3> >({x,y,endZ+1},{nx,ny,0})  );
							}
						}
				}
			}
			//now, lets seed the neighbors (in the 6 directions) of the filled cells (that, is, from [startZ] to [endZ])
			bool isertedInterval = false;
			if(y+1 < chunkEndY)
				for(int i=startZ;i<=endZ;i++) {
					//if its a new interval (whose seed was not inserted into the stack)
					if(cellsLabels.has2ndLevel(x,y+1,i)) {
						std::vector<std::vector<std::vector<int> > > &idsSecondLevel = cellsLabels.childGridLabels[x][y+1][i]->labels;
						
						for(int nx=0;nx<nestedGridSize;nx++) {
							isertedInterval = false;
							for(int nz=0;nz<nestedGridSize;nz++)
								if(idsSecondLevel[nx][0][nz]==oldId) {
									if(!isertedInterval) {	
										isertedInterval = true;
										seedsLevel2.push(pair<array<int,3>,array<int,3> >({x,y+1,i},{nx,0,nz}));									
									}									
								} else {
									isertedInterval = false;
								}
						}
						isertedInterval = false;
					} else {
						if( (idsFirstLevel[x][y+1][i]==oldId) ) {
							if(!isertedInterval) {
								seedsLevel1.push({x,y+1,i});
								isertedInterval = true;
							}
						} else {
							isertedInterval = false;
						}
					}
					
				}

			isertedInterval = false;		
			if(y-1 >= chunkStartY)
				for(int i=startZ;i<=endZ;i++) {
					//if its a new interval (whose seed was not inserted into the stack)
					if(cellsLabels.has2ndLevel(x,y-1,i)) {
						std::vector<std::vector<std::vector<int> > > &idsSecondLevel = cellsLabels.childGridLabels[x][y-1][i]->labels;
						
						for(int nx=0;nx<nestedGridSize;nx++) {
							isertedInterval = false;
							for(int nz=0;nz<nestedGridSize;nz++)
								if(idsSecondLevel[nx][nestedGridSize-1][nz]==oldId) {
									if(!isertedInterval) {	
										isertedInterval = true;
										seedsLevel2.push(pair<array<int,3>,array<int,3> >({x,y-1,i},{nx,nestedGridSize-1,nz}));									
									}									
								} else {
									isertedInterval = false;
								}
						}
						isertedInterval = false;
					} else {
						if( (idsFirstLevel[x][y-1][i]==oldId) ) {
							if(!isertedInterval) {
								seedsLevel1.push({x,y-1,i});
								isertedInterval = true;
							}
						} else {
							isertedInterval = false;
						}
					}
					
				}	

			isertedInterval = false;		
			if(x+1 < chunkEndX)
				for(int i=startZ;i<=endZ;i++) {
					//if its a new interval (whose seed was not inserted into the stack)
					if(cellsLabels.has2ndLevel(x+1,y,i)) {
						std::vector<std::vector<std::vector<int> > > &idsSecondLevel = cellsLabels.childGridLabels[x+1][y][i]->labels;
						
						for(int ny=0;ny<nestedGridSize;ny++) {
							isertedInterval = false;
							for(int nz=0;nz<nestedGridSize;nz++)
								if(idsSecondLevel[0][ny][nz]==oldId) {
									if(!isertedInterval) {	
										isertedInterval = true;
										seedsLevel2.push(pair<array<int,3>,array<int,3> >({x+1,y,i},{0,ny,nz}));									
									}									
								} else {
									isertedInterval = false;
								}
						}
						isertedInterval = false;
					} else {
						if( (idsFirstLevel[x+1][y][i]==oldId) ) {
							if(!isertedInterval) {
								seedsLevel1.push({x+1,y,i});
								isertedInterval = true;
							}
						} else {
							isertedInterval = false;
						}
					}
					
				}	




			isertedInterval = false;		
			if(x-1 >=chunkStartX)
				for(int i=startZ;i<=endZ;i++) {
					//if its a new interval (whose seed was not inserted into the stack)
					if(cellsLabels.has2ndLevel(x-1,y,i)) {
						std::vector<std::vector<std::vector<int> > > &idsSecondLevel = cellsLabels.childGridLabels[x-1][y][i]->labels;
						
						for(int ny=0;ny<nestedGridSize;ny++) {
							isertedInterval = false;
							for(int nz=0;nz<nestedGridSize;nz++)
								if(idsSecondLevel[nestedGridSize-1][ny][nz]==oldId) {
									if(!isertedInterval) {	
										isertedInterval = true;
										seedsLevel2.push(pair<array<int,3>,array<int,3> >({x-1,y,i},{nestedGridSize-1,ny,nz}));									
									}									
								} else {
									isertedInterval = false;
								}
						}
						isertedInterval = false;
					} else {
						if( (idsFirstLevel[x-1][y][i]==oldId) ) {
							if(!isertedInterval) {
								seedsLevel1.push({x-1,y,i});
								isertedInterval = true;
							}
						} else {
							isertedInterval = false;
						}
					}					
				}	
		} else {
			//cerr << "getting top second levle.." << endl;
			pair<array<int,3>,array<int,3> > &seed = seedsLevel2.top();
			seedsLevel2.pop();

			int x = seed.first[0];
			int y = seed.first[1];
			int z = seed.first[2];
			int nx = seed.second[0];
			int ny = seed.second[1];
			int nz = seed.second[2];

			//cerr << "Here..." << x << " " << y << " " << z << " " << nx << " " << ny << " " << nz << endl;
			assert(cellsLabels.has2ndLevel(x,y,z));
			if(cellsLabels.childGridLabels[x][y][z]->labels[nx][ny][nz] != oldId) continue;

			pair<int,int> endZ; //first is the first level z, second is the second level z...
			endZ.first = z;
			endZ.second = nz;

			//cerr << "Setting the line of cells in second level..." << endl;

			while(true) {
				cellsLabels.childGridLabels[x][y][endZ.first ]->labels[nx][ny][endZ.second] = id;
				if(endZ.second==nestedGridSize-1) { //we reached the end of a second level grid...
					if(endZ.first == chunkEndZ-1) {
						break;
					}
					if(!cellsLabels.has2ndLevel(x,y,endZ.first+1)) { // the next 1st level cell isn't nested... 
						break;
					}
					if(cellsLabels.childGridLabels[x][y][endZ.first+1]->labels[nx][ny][0]!=oldId) {
						break;
					} 
					endZ.first++;
					endZ.second=0;
				} else { //we havent reached the end of the nested grid yet... we just need to advance the nested cell z...
					if(cellsLabels.childGridLabels[x][y][endZ.first]->labels[nx][ny][endZ.second+1]!=oldId) {
						break;
					} 
					endZ.second++;
				}
			}
			//endZ will be the last nested cells that we set to id...

			pair<int,int> startZ; //first is the first level z, second is the second level z...
			startZ.first = z;
			startZ.second = nz;

			while(true) {
				cellsLabels.childGridLabels[x][y][startZ.first ]->labels[nx][ny][startZ.second] = id;
				if(startZ.second==0) { //we reached the end of a second level grid...
					if(startZ.first == chunkStartZ) {
						break;
					}
					if(!cellsLabels.has2ndLevel(x,y,startZ.first-1)) {
						break;
					}
					if(cellsLabels.childGridLabels[x][y][startZ.first-1]->labels[nx][ny][nestedGridSize-1]!=oldId) {
						break;
					} 
					startZ.first--;
					startZ.second=nestedGridSize-1;
				} else { //we havent reached the end of the nested grid yet... we just need to advance the nested cell z...
					if(cellsLabels.childGridLabels[x][y][startZ.first]->labels[nx][ny][startZ.second-1]!=oldId) {
						break;
					} 
					startZ.second--;
				}
			}
			//endZ will be the first nested cells that we set to id... that is, we've set the cells in [startZ,endZ]

			//cerr << "Line of cells set... now the next and previous in line.." << endl;

			//if next and second cells 
			if(startZ.second==0 && startZ.first!=chunkStartZ && !cellsLabels.has2ndLevel(x,y,startZ.first-1) && cellsLabels.labels[x][y][startZ.first-1]==oldId) {
				seedsLevel1.push({x,y,startZ.first-1});
			}
			if(endZ.second==nestedGridSize-1 && endZ.first!=chunkEndZ-1 && !cellsLabels.has2ndLevel(x,y,endZ.first+1) && cellsLabels.labels[x][y][endZ.first+1]==oldId) { // the next 1st level cell isn't nested... 
				seedsLevel1.push({x,y,endZ.first+1});
			}







			//cerr << "Setting around..." << endl;
			bool isertedInterval = false;
			if (ny+1<nestedGridSize) {
				//the y above is inside the same first level grid... all the grids until at least endZ have two levels!
				pair<int,int> currZ = startZ;
				while(  true  ) {
					if(currZ.first > endZ.first || (currZ.first==endZ.first && currZ.second > endZ.second)) break;
					if(  cellsLabels.childGridLabels[x][y][currZ.first]->labels[nx][ny+1][currZ.second]==oldId  ) {
						if(!isertedInterval) {
							assert(cellsLabels.has2ndLevel(x,y,currZ.first));
							seedsLevel2.push(pair<array<int,3>,array<int,3> >({x,y,currZ.first},{nx,ny+1,currZ.second}));
							isertedInterval= true;
						}
					} else isertedInterval= false;

					currZ.second++;
					if(currZ.second==nestedGridSize) {
						currZ.second = 0;
						currZ.first++;
					}
				}
			} else if(y+1 <chunkEndY){
				//the next neighbor y is actually in another first level grid cell...
				//we need to treat this...
				pair<int,int> currZ = startZ;

				bool insertedFirstLevel = false;
				bool insertedSecondLevel = false;
				while(true) {
					if(currZ.first > endZ.first || (currZ.first==endZ.first && currZ.second > endZ.second)) break;
					if(cellsLabels.has2ndLevel(x,y+1,currZ.first)) {
						insertedFirstLevel = false;
						if(  cellsLabels.childGridLabels[x][y+1][currZ.first]->labels[nx][0][currZ.second]==oldId  ) {
							if(!insertedSecondLevel) {
								assert(cellsLabels.has2ndLevel(x,y+1,currZ.first));
								seedsLevel2.push(pair<array<int,3>,array<int,3> >({x,y+1,currZ.first},{nx,0,currZ.second}));
								insertedSecondLevel= true;
							}
						} else insertedSecondLevel= false;

						currZ.second++;
						if(currZ.second==nestedGridSize) {
							currZ.second = 0;
							currZ.first++;
						}
					} else { //the top cell doesn't have a second level...
						insertedSecondLevel = false;
						if(  cellsLabels.labels[x][y+1][currZ.first]==oldId  ) {
							if(!insertedFirstLevel) {
								seedsLevel1.push({x,y+1,currZ.first});
								insertedFirstLevel= true;
							}
						} else insertedFirstLevel= false;

						currZ.first++;
						currZ.second=0;
					}					
				}
			}
			//cerr << "Finish" << endl << endl;

			//cerr << "Inesrting..." << endl;
			isertedInterval = false;
			if (ny-1>=0) {
				//the y below is inside the same first level grid... all the grids until at least endZ have two levels!
				pair<int,int> currZ = startZ;
				while(  true  ) {
					if(currZ.first > endZ.first || (currZ.first==endZ.first && currZ.second > endZ.second)) break;
					if(  cellsLabels.childGridLabels[x][y][currZ.first]->labels[nx][ny-1][currZ.second]==oldId  ) {
						if(!isertedInterval) {
							assert(cellsLabels.has2ndLevel(x,y,currZ.first));
							seedsLevel2.push(pair<array<int,3>,array<int,3> >({x,y,currZ.first},{nx,ny-1,currZ.second}));
							isertedInterval= true;
						}
					} else isertedInterval= false;

					currZ.second++;
					if(currZ.second==nestedGridSize) {
						currZ.second = 0;
						currZ.first++;
					}
				}
			} else if(y-1>=chunkStartY) {
				//the next neighbor y is actually in another first level grid cell...
				//we need to treat this...
				pair<int,int> currZ = startZ;

				bool insertedFirstLevel = false;
				bool insertedSecondLevel = false;
				while(true) {
					if(currZ.first > endZ.first || (currZ.first==endZ.first && currZ.second > endZ.second)) break;
					if(cellsLabels.has2ndLevel(x,y-1,currZ.first)) {
						insertedFirstLevel = false;
						if(  cellsLabels.childGridLabels[x][y-1][currZ.first]->labels[nx][nestedGridSize-1][currZ.second]==oldId  ) {
							if(!insertedSecondLevel) {
								assert(cellsLabels.has2ndLevel(x,y-1,currZ.first));
								seedsLevel2.push(pair<array<int,3>,array<int,3> >({x,y-1,currZ.first},{nx,nestedGridSize-1,currZ.second}));
								insertedSecondLevel= true;
							}
						} else insertedSecondLevel= false;

						currZ.second++;
						if(currZ.second==nestedGridSize) {
							currZ.second = 0;
							currZ.first++;
						}
					} else { //the bottom cell doesn't have a second level...
						insertedSecondLevel = false;
						if(  cellsLabels.labels[x][y-1][currZ.first]==oldId  ) {
							if(!insertedFirstLevel) {
								seedsLevel1.push({x,y-1,currZ.first});
								insertedFirstLevel= true;
							}
						} else insertedFirstLevel= false;

						currZ.first++;
						currZ.second=0;
					}					
				}
			}

			//Next and previous x directions...

			isertedInterval = false;
			if (nx+1<nestedGridSize) {
				//the y above is inside the same first level grid... all the grids until at least endZ have two levels!
				pair<int,int> currZ = startZ;
				while(  true  ) {
					if(currZ.first > endZ.first || (currZ.first==endZ.first && currZ.second > endZ.second)) break;
					if(  cellsLabels.childGridLabels[x][y][currZ.first]->labels[nx+1][ny][currZ.second]==oldId  ) {
						if(!isertedInterval) {
							assert(cellsLabels.has2ndLevel(x,y,currZ.first));
							seedsLevel2.push(pair<array<int,3>,array<int,3> >({x,y,currZ.first},{nx+1,ny,currZ.second}));
							isertedInterval= true;
						}
					} else isertedInterval= false;

					currZ.second++;
					if(currZ.second==nestedGridSize) {
						currZ.second = 0;
						currZ.first++;
					}
				}
			} else if(x+1 <chunkEndX){
				//the next neighbor y is actually in another first level grid cell...
				//we need to treat this...
				pair<int,int> currZ = startZ;

				bool insertedFirstLevel = false;
				bool insertedSecondLevel = false;
				while(true) {
					if(currZ.first > endZ.first || (currZ.first==endZ.first && currZ.second > endZ.second)) break;
					if(cellsLabels.has2ndLevel(x+1,y,currZ.first)) {
						insertedFirstLevel = false;
						if(  cellsLabels.childGridLabels[x+1][y][currZ.first]->labels[0][ny][currZ.second]==oldId  ) {
							if(!insertedSecondLevel) {
								assert(cellsLabels.has2ndLevel(x+1,y,currZ.first));
								seedsLevel2.push(pair<array<int,3>,array<int,3> >({x+1,y,currZ.first},{0,ny,currZ.second}));
								insertedSecondLevel= true;
							}
						} else insertedSecondLevel= false;

						currZ.second++;
						if(currZ.second==nestedGridSize) {
							currZ.second = 0;
							currZ.first++;
						}
					} else { //the top cell doesn't have a second level...
						insertedSecondLevel = false;
						if(  cellsLabels.labels[x+1][y][currZ.first]==oldId  ) {
							if(!insertedFirstLevel) {
								seedsLevel1.push({x+1,y,currZ.first});
								insertedFirstLevel= true;
							}
						} else insertedFirstLevel= false;

						currZ.first++;
						currZ.second=0;
					}					
				}
			}

			isertedInterval = false;
			if (nx-1>=0) {
				//the y below is inside the same first level grid... all the grids until at least endZ have two levels!
				pair<int,int> currZ = startZ;
				while(  true  ) {
					if(currZ.first > endZ.first || (currZ.first==endZ.first && currZ.second > endZ.second)) break;
					if(  cellsLabels.childGridLabels[x][y][currZ.first]->labels[nx-1][ny][currZ.second]==oldId  ) {
						if(!isertedInterval) {
							assert(cellsLabels.has2ndLevel(x,y,currZ.first));
							seedsLevel2.push(pair<array<int,3>,array<int,3> >({x,y,currZ.first},{nx-1,ny,currZ.second}));
							isertedInterval= true;
						}
					} else isertedInterval= false;

					currZ.second++;
					if(currZ.second==nestedGridSize) {
						currZ.second = 0;
						currZ.first++;
					}
				}
			} else if(x-1>=chunkStartX) {
				//the next neighbor y is actually in another first level grid cell...
				//we need to treat this...
				pair<int,int> currZ = startZ;

				bool insertedFirstLevel = false;
				bool insertedSecondLevel = false;
				while(true) {
					if(currZ.first > endZ.first || (currZ.first==endZ.first && currZ.second > endZ.second)) break;
					if(cellsLabels.has2ndLevel(x-1,y,currZ.first)) {
						insertedFirstLevel = false;
						if(  cellsLabels.childGridLabels[x-1][y][currZ.first]->labels[nestedGridSize-1][ny][currZ.second]==oldId  ) {
							if(!insertedSecondLevel) {
								assert(cellsLabels.has2ndLevel(x-1,y,currZ.first));
								seedsLevel2.push(pair<array<int,3>,array<int,3> >({x-1,y,currZ.first},{nestedGridSize-1,ny,currZ.second}));
								insertedSecondLevel= true;
							}
						} else insertedSecondLevel= false;

						currZ.second++;
						if(currZ.second==nestedGridSize) {
							currZ.second = 0;
							currZ.first++;
						}
					} else { //the bottom cell doesn't have a second level...
						insertedSecondLevel = false;
						if(  cellsLabels.labels[x-1][y][currZ.first]==oldId  ) {
							if(!insertedFirstLevel) {
								seedsLevel1.push({x-1,y,currZ.first});
								insertedFirstLevel= true;
							}
						} else insertedFirstLevel= false;

						currZ.first++;
						currZ.second=0;
					}					
				}
			}



			//now we just need to seed the grid lines around the line [startZ,endZ] (in the 6 directions...)
		}

		/*
		array<int,3> seed = seeds.top();
		seeds.pop();

		int x = seed[0];
		int y = seed[1];
		int z = seed[2];

		if(idsOfCellsTemp[x][y][z] != oldId) continue;

		//cout << x << " " << y << " " << z << endl;

		int i;
		for(i=z;i<gridSize && idsOfCellsTemp[x][y][i]==oldId;i++) {
			idsOfCellsTemp[x][y][i] = id;
		}
		// position (i-1) is the last one that will be filled
		int endZ = i-1;
		for(i=z-1;i>=0 && idsOfCellsTemp[x][y][i]==oldId;i--) {
			idsOfCellsTemp[x][y][i] = id;
		}
		// (i+1) is the last z that will be filled (in reverse order..)
		int startZ = i+1;
		//so, we filled the interval [startZ, endZ] with id...
		//lets seed the neighbors...

		//top and bottom y...
		bool isertedInterval = false;
		if(y+1 < gridSize)
			for(int i=startZ;i<=endZ;i++) {
				//if its a new interval (whose seed was not inserted into the stack)
				if( (idsOfCellsTemp[x][y+1][i]==oldId) ) {
					if(!isertedInterval) {
						seeds.push({x,y+1,i});
						isertedInterval = true;
					}
				} else {
					isertedInterval = false;
				}
			}

		isertedInterval = false;

		if(y-1>=0)
			for(int i=startZ;i<=endZ;i++) {
				//if its a new interval (whose seed was not inserted into the stack)
				if( (idsOfCellsTemp[x][y-1][i]==oldId)  ) {
					if(!isertedInterval) {
						seeds.push({x,y-1,i});
						isertedInterval = true;
					}
				} else {
					isertedInterval = false;
				}
			}

		//left and right x...
		isertedInterval = false;
		if(x+1 < gridSize)
			for(int i=startZ;i<=endZ;i++) {
				//if its a new interval (whose seed was not inserted into the stack)
				if( (idsOfCellsTemp[x+1][y][i]==oldId) ) {
					if(!isertedInterval) {
						seeds.push({x+1,y,i});
						isertedInterval = true;
					}
				} else {
					isertedInterval = false;
				}
			}

		isertedInterval = false;
		if(x-1>=0)
			for(int i=startZ;i<=endZ;i++) {
				//if its a new interval (whose seed was not inserted into the stack)
				if( (idsOfCellsTemp[x-1][y][i]==oldId)  ) {
					if(!isertedInterval) {
						seeds.push({x-1,y,i});
						isertedInterval = true;
					}
				} else {
					isertedInterval = false;
				}
			}
		*/

	}

} 