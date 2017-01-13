
array<VertCoord,3> MeshIntersectionGeometry::coordRangeMeshes() const {
  return {boundingBoxTwoMeshesTogetter[1][0]-boundingBoxTwoMeshesTogetter[0][0],boundingBoxTwoMeshesTogetter[1][1]-boundingBoxTwoMeshesTogetter[0][1],boundingBoxTwoMeshesTogetter[1][2]-boundingBoxTwoMeshesTogetter[0][2]};
} 

array<int,3> MeshIntersectionGeometry::getGridCellContainingVertex(const int meshId, const int vertexId, const VertCoord &cellScale, TempVarsGetGridCellContainingVertex &tempVars ) const {
	VertCoord &tempVar = tempVars.tempVertCoords;//tempVertCoords[poolToUse];
	const vector<Point> &points = verticesCoordinates[meshId];

  const Point &boundingBoxMin = boundingBoxTwoMeshesTogetter[0];

	tempVar = points[vertexId][0];
  tempVar -= boundingBoxMin[0];
  tempVar *= cellScale;
  const int x = convertToInt(tempVar,tempVars.tempVarsInt);

  tempVar = points[vertexId][1];
  tempVar -= boundingBoxMin[1];
  tempVar *= cellScale;
  const int y = convertToInt(tempVar,tempVars.tempVarsInt);

  tempVar = points[vertexId][2];
  tempVar -= boundingBoxMin[2];
  tempVar *= cellScale;
  const int z  = convertToInt(tempVar,tempVars.tempVarsInt);

  return {x,y,z};
}

int MeshIntersectionGeometry::getGridCellXContainingVertex(int meshId,const VertCoord &xCoord, const VertCoord &cellScale, TempVarsGetGridCellContainingVertex &tempVars ) const {
  VertCoord &tempVar = tempVars.tempVertCoords;
  const Point &boundingBoxMin = boundingBoxTwoMeshesTogetter[0];

  tempVar = xCoord;
  tempVar -= boundingBoxMin[0];
  tempVar *= cellScale;
  const int x = convertToInt(tempVar,tempVars.tempVarsInt);

  return x;
}
int MeshIntersectionGeometry::getGridCellYContainingVertex(int meshId,const VertCoord &yCoord, const VertCoord &cellScale, TempVarsGetGridCellContainingVertex &tempVars ) const {
  VertCoord &tempVar = tempVars.tempVertCoords;
  const Point &boundingBoxMin = boundingBoxTwoMeshesTogetter[0];

  tempVar = yCoord;
  tempVar -= boundingBoxMin[1];
  tempVar *= cellScale;
  const int y = convertToInt(tempVar,tempVars.tempVarsInt);

  return y;
}
int MeshIntersectionGeometry::getGridCellZContainingVertex(int meshId,const VertCoord &zCoord, const VertCoord &cellScale, TempVarsGetGridCellContainingVertex &tempVars ) const {
  VertCoord &tempVar = tempVars.tempVertCoords;
  const Point &boundingBoxMin = boundingBoxTwoMeshesTogetter[0];

  tempVar = zCoord;
  tempVar -= boundingBoxMin[2];
  tempVar *= cellScale;
  const int z  = convertToInt(tempVar,tempVars.tempVarsInt);

  return z;
}