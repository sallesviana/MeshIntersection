#include "sosPredicatesImpl.h"

/*
t00x -> 0, t00y -> 0, t00z -> 1, t01x -> 10, t01y -> 0, t01z -> 1, \
t02x -> 0, t02y -> 10, t02z -> 1, r00x -> 1, r00y -> 1, r00z -> 0, \
r01x -> 1, r01y -> 1, r01z -> h, t10x -> 0, t10y -> 0, t10z -> 1, \
t11x -> 10, t11y -> 0, t11z -> 1, t12x -> 0, t12y -> 10, t12z -> 1, \
r10x -> 2, r10y -> 1, r10z -> 0, r11x -> 1, r11y -> 1, r11z -> h, \
t20x -> 0, t20y -> 0, t20z -> 1, t21x -> 10, t21y -> 0, t21z -> 1, \
t22x -> 0, t22y -> 10, t22z -> 1, r20x -> 1, r20y -> 2, r20z -> 0, \
r21x -> 1, r21y -> 1, r21z -> 1
*/

int main() {
	Point t0 = {0,0,0};
	Point t1 = {10,0,0};
	Point t2 = {0,10,0};

	cerr << t0[0] << endl;
	cerr << t0[2] << endl;

	Triangle t;
	t.points[0] = t0;
	t.points[1] = t1;
	t.points[2] = t2;		

	Point r00 = {1,1,0};
	Point r01 = {1,1,1};
	Point r10 = {2,1,0};
	Point r11 = r01;
	Point r20 = {1,2,0};
	Point r21 = r01;

	VertexFromIntersection v0,v1,v2;
	v0.triangle = v1.triangle = v2.triangle = t;
	v0.edge[0] =  r00;
	v0.edge[1] =  r01;

	cerr << r01[0] << endl;

	v1.edge[0] =  r10;
	v1.edge[1] =  r11;

	v2.edge[0] =  r20;
	v2.edge[1] =  r21;		

	SosPredicatesImpl test;

	cout << test.orient2D_z0_000(v0, v1, v2) << endl;

}