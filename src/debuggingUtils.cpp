void printVertexForDebugging(VertCoord v[3]) {
	cerr << "Vertex: [" << v[0].get_d() << " "<< v[1].get_d()<< " "  <<v[2].get_d() << "]\n";
}

void printVertexForDebugging(array<VertCoord,3> &v) {
	cerr << "Vertex: [" << v[0].get_d() << " "<< v[1].get_d()<< " "  <<v[2].get_d() << "]\n";
}