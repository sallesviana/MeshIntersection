#ifndef SosPredicatesImpl_H
#define SosPredicatesImpl_H

#include <iostream>
#include <vector>
#include <array>

#include <gmpxx.h>


using namespace std;


typedef  mpq_class rational;
typedef  mpz_class big_int;

typedef rational VertCoord;
typedef array<VertCoord,3> Point;


class InputVertex {
	public:
	InputVertex(Point &p1): p(p1){};
	InputVertex() {}
	Point p;
};

class Triangle {
public:
	InputVertex points[3];
	const InputVertex * getInputVertex(int id) const {
		//cerr << points[id][0] << " " << "" << " " << endl;
		return &points[id];
	}
};

class VertexFromIntersection {
	public:
	InputVertex edge[2];
	Triangle triangle;
};

class SosPredicatesImpl {



public:

const Point& getCoordinates(const InputVertex &iv) const {
	return iv.p;
}


int orient2D_x0_000(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_y0_000(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_z0_000(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_x0_001(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_y0_001(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_z0_001(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_x0_010(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_y0_010(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_z0_010(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_x0_011(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_y0_011(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_z0_011(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_x0_100(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_y0_100(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_z0_100(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_x0_101(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_y0_101(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_z0_101(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_x0_110(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_y0_110(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_z0_110(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_x0_111(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_y0_111(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_z0_111(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
int orient2D_x0_000(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_000(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_000(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_001(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_001(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_001(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_010(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_010(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_010(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_011(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_011(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_011(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_100(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_100(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_100(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_101(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_101(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_101(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_110(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_110(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_110(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_111(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_111(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_111(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_000(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_000(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_000(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_001(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_001(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_001(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_010(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_010(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_010(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_011(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_011(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_011(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_100(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_100(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_100(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_101(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_101(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_101(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_110(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_110(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_110(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_111(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_111(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_111(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_000(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_000(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_000(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_001(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_001(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_001(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_010(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_010(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_010(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_011(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_011(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_011(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_100(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_100(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_100(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_101(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_101(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_101(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_110(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_110(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_110(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_x0_111(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_y0_111(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
int orient2D_z0_111(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;

};

#endif