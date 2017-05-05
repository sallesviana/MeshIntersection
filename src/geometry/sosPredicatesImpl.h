#include "../3dGeometry.h"

#ifndef SosPredicatesImpl_H
#define SosPredicatesImpl_H

using namespace std;


class SosPredicatesImpl {

public:
	SosPredicatesImpl(MeshIntersectionGeometry *geom,TempVarsSoSPredicatesImpl &tempVarsArg): tempVars(tempVarsArg), geometry(geom) {}

	int orientation1D(const InputVertex &v0, const InputVertex &v1,int whatAxisProjectTo) const;
	int orientation1D(const InputVertex &v0, const VertexFromIntersection &v1,int whatAxisProjectTo) const;
	int orientation1D(const VertexFromIntersection &v0, const VertexFromIntersection &v1,int whatAxisProjectTo) const;

	int orientation2D(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2,int whatPlaneProjectTo) const;
	int orientation2D(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2,int whatPlaneProjectTo) const;
	int orientation2D(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2,int whatPlaneProjectTo) const;
	int orientation2D(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2,int whatPlaneProjectTo) const;

	int orientation3D(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const InputVertex &v3) const;
	int orientation3D(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &v3) const;

private:
	MeshIntersectionGeometry *geometry;

	const Point& getCoordinates(const Vertex &iv) const {
		return geometry->getCoordinates(iv);
	}

	
	TempVarsSoSPredicatesImpl &tempVars;

	//given a vertex from intersection, returns the epsilon coefficient of a given coordinate
	//coord may be COORD_X,COORD_Y or COORD_Z
	//epsCoefficient may be 0, 1, 2 or 3
	const VertCoord&  getEpsCoefficientsVertexFromIntersection(const VertexFromIntersection &v0, int epsCoefficient,int coord) const;
	VertCoord getEpsCoefficientsVertexFromIntersectionFaster(const VertexFromIntersection &v0, int epsCoefficient,int coord) const;
	const VertCoord&  getEpsCoefficientsVertexFromIntersectionFaster2(const VertexFromIntersection &v0, int epsCoefficient,int coord) const;
	VertCoord getEpsCoefficientsVertexFromIntersectionOriginal(const VertexFromIntersection &v0, int epsCoefficient,int coord) const;





	int orientation1D_x(const InputVertex &v0, const InputVertex &v1) const;
	int orientation1D_x(const InputVertex &v0, const VertexFromIntersection &v1) const;
	int orientation1D_x(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const;

	int orientation2D_x0(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
	int orientation2D_x0(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
	int orientation2D_x0(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
	int orientation2D_x0(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;

	int orientation1D_y(const InputVertex &v0, const InputVertex &v1) const;
	int orientation1D_y(const InputVertex &v0, const VertexFromIntersection &v1) const;
	int orientation1D_y(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const;

	int orientation2D_y0(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
	int orientation2D_y0(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
	int orientation2D_y0(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
	int orientation2D_y0(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;

	int orientation1D_z(const InputVertex &v0, const InputVertex &v1) const;
	int orientation1D_z(const InputVertex &v0, const VertexFromIntersection &v1) const;
	int orientation1D_z(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const;

	int orientation2D_z0(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const;
	int orientation2D_z0(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const;
	int orientation2D_z0(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;
	int orientation2D_z0(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const;



	int orient1D_x_00(const InputVertex &v0, const InputVertex &v1) const;
	int orient1D_y_00(const InputVertex &v0, const InputVertex &v1) const;
	int orient1D_z_00(const InputVertex &v0, const InputVertex &v1) const;
	int orient1D_x_01(const InputVertex &v0, const InputVertex &v1) const;
	int orient1D_y_01(const InputVertex &v0, const InputVertex &v1) const;
	int orient1D_z_01(const InputVertex &v0, const InputVertex &v1) const;
	int orient1D_x_10(const InputVertex &v0, const InputVertex &v1) const;
	int orient1D_y_10(const InputVertex &v0, const InputVertex &v1) const;
	int orient1D_z_10(const InputVertex &v0, const InputVertex &v1) const;
	int orient1D_x_11(const InputVertex &v0, const InputVertex &v1) const;
	int orient1D_y_11(const InputVertex &v0, const InputVertex &v1) const;
	int orient1D_z_11(const InputVertex &v0, const InputVertex &v1) const;
	int orient1D_x_00(const InputVertex &v0, const VertexFromIntersection &v1) const;
	int orient1D_y_00(const InputVertex &v0, const VertexFromIntersection &v1) const;
	int orient1D_z_00(const InputVertex &v0, const VertexFromIntersection &v1) const;
	int orient1D_x_01(const InputVertex &v0, const VertexFromIntersection &v1) const;
	int orient1D_y_01(const InputVertex &v0, const VertexFromIntersection &v1) const;
	int orient1D_z_01(const InputVertex &v0, const VertexFromIntersection &v1) const;
	int orient1D_x_10(const InputVertex &v0, const VertexFromIntersection &v1) const;
	int orient1D_y_10(const InputVertex &v0, const VertexFromIntersection &v1) const;
	int orient1D_z_10(const InputVertex &v0, const VertexFromIntersection &v1) const;
	int orient1D_x_11(const InputVertex &v0, const VertexFromIntersection &v1) const;
	int orient1D_y_11(const InputVertex &v0, const VertexFromIntersection &v1) const;
	int orient1D_z_11(const InputVertex &v0, const VertexFromIntersection &v1) const;
	int orient1D_x_00(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const;
	int orient1D_y_00(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const;
	int orient1D_z_00(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const;
	int orient1D_x_01(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const;
	int orient1D_y_01(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const;
	int orient1D_z_01(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const;
	int orient1D_x_10(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const;
	int orient1D_y_10(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const;
	int orient1D_z_10(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const;
	int orient1D_x_11(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const;
	int orient1D_y_11(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const;
	int orient1D_z_11(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const;
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
	int orient3D_0000(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const InputVertex &v3) const;
	int orient3D_0001(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const InputVertex &v3) const;
	int orient3D_0010(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const InputVertex &v3) const;
	int orient3D_0011(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const InputVertex &v3) const;
	int orient3D_0100(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const InputVertex &v3) const;
	int orient3D_0101(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const InputVertex &v3) const;
	int orient3D_0110(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const InputVertex &v3) const;
	int orient3D_0111(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const InputVertex &v3) const;
	int orient3D_1000(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const InputVertex &v3) const;
	int orient3D_1001(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const InputVertex &v3) const;
	int orient3D_1010(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const InputVertex &v3) const;
	int orient3D_1011(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const InputVertex &v3) const;
	int orient3D_1100(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const InputVertex &v3) const;
	int orient3D_1101(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const InputVertex &v3) const;
	int orient3D_1110(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const InputVertex &v3) const;
	int orient3D_1111(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const InputVertex &v3) const;
	int orient3D_0000(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &v3) const;
	int orient3D_0001(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &v3) const;
	int orient3D_0010(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &v3) const;
	int orient3D_0011(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &v3) const;
	int orient3D_0100(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &v3) const;
	int orient3D_0101(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &v3) const;
	int orient3D_0110(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &v3) const;
	int orient3D_0111(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &v3) const;
	int orient3D_1000(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &v3) const;
	int orient3D_1001(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &v3) const;
	int orient3D_1010(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &v3) const;
	int orient3D_1011(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &v3) const;
	int orient3D_1100(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &v3) const;
	int orient3D_1101(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &v3) const;
	int orient3D_1110(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &v3) const;
	int orient3D_1111(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &v3) const;
};

#endif
