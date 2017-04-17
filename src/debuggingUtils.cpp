void printVertexForDebugging(VertCoord v[3]) {
	cerr << "Vertex: [" << v[0].get_d() << " "<< v[1].get_d()<< " "  <<v[2].get_d() << "]\n";
}

void printVertexForDebugging(array<VertCoord,3> &v) {
	cerr << "Vertex: [" << v[0].get_d() << " "<< v[1].get_d()<< " "  <<v[2].get_d() << "]\n";
}

array<double,3> toDoublePoint(const Point &p) {
	array<double,3> pnew;
  for(int i=0;i<3;i++) {
    pnew[i] = p[i].get_d();
  }
  return pnew;
}


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



void storeOneTriangleAsGts(const string &path,array<double,3> a,array<double,3> b, array<double,3> c,bool reverseOrientation) {
	ofstream fout(path.c_str());
	int numVert = 3;
  int numEdges = 3;
  int numFaces = 1;

	fout << numVert << " " << numEdges << " " << numFaces << "\n";
  fout << a[0] << " " << a[1] << " " << a[2] << "\n";
  fout << b[0] << " " << b[1] << " " << b[2] << "\n";
  fout << c[0] << " " << c[1] << " " << c[2] << "\n";
  
  if(!reverseOrientation) {
  	fout << 1 << " " << 2 << "\n";
  	fout << 2 << " " << 3 << "\n";
  	fout << 3 << " " << 1 << "\n";
  	fout << "1 2 3"<< "\n";   	  	  	
  } else {
  	fout << 1 << " " << 3 << "\n";
  	fout << 3 << " " << 2 << "\n";
  	fout << 2 << " " << 1 << "\n";
  	fout << "1 2 3"<< "\n";   	
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


/*
  const Triangle *tA = intersectingTrianglesThatGeneratedEdges[i].first;
  const Triangle *tB = intersectingTrianglesThatGeneratedEdges[i].second;

*/
void printTriangleWithManyIntersection(const Triangle &t,int meshOfT,const vector<int> &edgesFromIntersection,const vector< pair<Triangle *,Triangle *> > &intersectingTrianglesThatGeneratedEdges) {
  cout << "Triangle with many intersection" << "\n";
  cout << "Num intersections: " << edgesFromIntersection.size() << "\n";
  cout << "Mesh of t: " << meshOfT << "\n";
  cout << t.p[0] << " " << t.p[1] << " " << t.p[2] << "\n";


  for(int edge:edgesFromIntersection) {
    const Triangle *ts;
    if(meshOfT==0) {
      ts = intersectingTrianglesThatGeneratedEdges[edge].second;
    } else {
      ts = intersectingTrianglesThatGeneratedEdges[edge].first;
    }
    cout << ts->p[0] << " " << ts->p[1] << " " << ts->p[2] << "\n";
  }

} 

void saveEdgesAsGTS(const vector< pair< array<VertCoord,3>,array<VertCoord,3> > > &edges,const vector<int> &edgesToSelect,const string &path) {
  vector<pair<array<double,3>,array<double,3>> > edgesToStore;
  for(int edgeId:edgesToSelect) {
    array<double,3> v0 = {edges[edgeId].first[0].get_d(),edges[edgeId].first[1].get_d(),edges[edgeId].first[2].get_d()};
    array<double,3> v1 = {edges[edgeId].second[0].get_d(),edges[edgeId].second[1].get_d(),edges[edgeId].second[2].get_d()};
    edgesToStore.push_back({v0,v1});
  }

  storeEdgesAsGts(path,edgesToStore );
}


void saveEdgesAsGTS(const MeshIntersectionGeometry &geom, const vector<pair<const Vertex *,const Vertex *>>  &edges,const string &path) {
  vector<pair<array<double,3>,array<double,3>> > edgesToStore;
  for(const pair<const Vertex *,const Vertex *> &edge:edges) {

    array<double,3> v0 = geom.getCoordinatesForDebugging(edge.first);
    array<double,3> v1 = geom.getCoordinatesForDebugging(edge.second);
    edgesToStore.push_back({v0,v1});
  }

  storeEdgesAsGts(path,edgesToStore );
}


void saveEdgesAsGTS(const vector< pair< int, int> > &edgesWithVertexIds, int meshWhereTriangleIs,  const string &path) {
  vector<pair<array<double,3>,array<double,3>> > edgesInThisTriangleDouble;
  for(const auto &e:edgesWithVertexIds) {
    const Point &v1 = *getPointFromVertexId(e.first,meshWhereTriangleIs);
    const Point &v2 = *getPointFromVertexId(e.second,meshWhereTriangleIs);
    
    edgesInThisTriangleDouble.push_back(make_pair(toDoublePoint(v1),toDoublePoint(v2)));
  }
  storeEdgesAsGts(path,edgesInThisTriangleDouble );
}
