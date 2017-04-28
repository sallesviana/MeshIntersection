#include "sosPredicatesImpl.h"

int SosPredicatesImpl::orient2D_y0_010(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const{


/*****************************************************/
const VertCoord &iv0x = getCoordinates(v0)[0];
const VertCoord &iv0y = getCoordinates(v0)[1];
const VertCoord &iv0z = getCoordinates(v0)[2];

const VertCoord &iv1x = getCoordinates(v1)[0];
const VertCoord &iv1y = getCoordinates(v1)[1];
const VertCoord &iv1z = getCoordinates(v1)[2];

const VertCoord &v2t0x = getCoordinates( *(v2.triangle.getInputVertex(0)) )[0];
const VertCoord &v2t1x = getCoordinates( *(v2.triangle.getInputVertex(1)) )[0];
const VertCoord &v2t2x = getCoordinates( *(v2.triangle.getInputVertex(2)) )[0];
const VertCoord &v2r0x = getCoordinates(v2.edge[0])[0];
const VertCoord &v2r1x = getCoordinates(v2.edge[1])[0];
const VertCoord &v2t0y = getCoordinates( *(v2.triangle.getInputVertex(0)) )[1];
const VertCoord &v2t1y = getCoordinates( *(v2.triangle.getInputVertex(1)) )[1];
const VertCoord &v2t2y = getCoordinates( *(v2.triangle.getInputVertex(2)) )[1];
const VertCoord &v2r0y = getCoordinates(v2.edge[0])[1];
const VertCoord &v2r1y = getCoordinates(v2.edge[1])[1];
const VertCoord &v2t0z = getCoordinates( *(v2.triangle.getInputVertex(0)) )[2];
const VertCoord &v2t1z = getCoordinates( *(v2.triangle.getInputVertex(1)) )[2];
const VertCoord &v2t2z = getCoordinates( *(v2.triangle.getInputVertex(2)) )[2];
const VertCoord &v2r0z = getCoordinates(v2.edge[0])[2];
const VertCoord &v2r1z = getCoordinates(v2.edge[1])[2];

/*****************************************************/
VertCoord coeffV2ePower0x = (v2r0z*v2r1x*(-(v2t1y*v2t2x) + v2t0y*(-v2t1x + v2t2x) + v2t0x*(v2t1y - v2t2y) + v2t1x*v2t2y) + v2r1x*(v2t0z*v2t1y*v2t2x - v2t0y*v2t1z*v2t2x - v2t0z*v2t1x*v2t2y + v2t0x*v2t1z*v2t2y + v2t0y*v2t1x*v2t2z - v2t0x*v2t1y*v2t2z + v2r0y*(-(v2t0x*v2t1z) + v2t0z*(v2t1x - v2t2x) + v2t1z*v2t2x + v2t0x*v2t2z - v2t1x*v2t2z)) + v2r0x*(-(v2t0z*v2t1y*v2t2x) + v2t0y*v2t1z*v2t2x + v2t0z*v2t1x*v2t2y - v2t0x*v2t1z*v2t2y + v2r1z*(-(v2t0x*v2t1y) + v2t0y*(v2t1x - v2t2x) + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) - v2t0y*v2t1x*v2t2z + v2t0x*v2t1y*v2t2z + v2r1y*(-(v2t0z*v2t1x) + v2t0x*v2t1z + v2t0z*v2t2x - v2t1z*v2t2x - v2t0x*v2t2z + v2t1x*v2t2z)))/(v2r0y*v2t0z*v2t1x - v2r1y*v2t0z*v2t1x - v2r0x*v2t0z*v2t1y + v2r1x*v2t0z*v2t1y - v2r0y*v2t0x*v2t1z + v2r1y*v2t0x*v2t1z + v2r0x*v2t0y*v2t1z - v2r1x*v2t0y*v2t1z - v2r0y*v2t0z*v2t2x + v2r1y*v2t0z*v2t2x + v2r0y*v2t1z*v2t2x - v2r1y*v2t1z*v2t2x + v2r0x*v2t0z*v2t2y - v2r1x*v2t0z*v2t2y - v2r0x*v2t1z*v2t2y + v2r1x*v2t1z*v2t2y + v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + v2r0y*v2t0x*v2t2z - v2r1y*v2t0x*v2t2z - v2r0x*v2t0y*v2t2z + v2r1x*v2t0y*v2t2z - v2r0y*v2t1x*v2t2z + v2r1y*v2t1x*v2t2z + v2r0x*v2t1y*v2t2z - v2r1x*v2t1y*v2t2z);
VertCoord coeffV2ePower1x = (v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + (v2r0y - v2r1y)*(v2t0z*v2t1x - v2t0x*v2t1z - v2t0z*v2t2x + v2t1z*v2t2x + v2t0x*v2t2z - v2t1x*v2t2z))/(v2r0y*v2t0z*v2t1x - v2r1y*v2t0z*v2t1x - v2r0x*v2t0z*v2t1y + v2r1x*v2t0z*v2t1y - v2r0y*v2t0x*v2t1z + v2r1y*v2t0x*v2t1z + v2r0x*v2t0y*v2t1z - v2r1x*v2t0y*v2t1z - v2r0y*v2t0z*v2t2x + v2r1y*v2t0z*v2t2x + v2r0y*v2t1z*v2t2x - v2r1y*v2t1z*v2t2x + v2r0x*v2t0z*v2t2y - v2r1x*v2t0z*v2t2y - v2r0x*v2t1z*v2t2y + v2r1x*v2t1z*v2t2y + v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + v2r0y*v2t0x*v2t2z - v2r1y*v2t0x*v2t2z - v2r0x*v2t0y*v2t2z + v2r1x*v2t0y*v2t2z - v2r0y*v2t1x*v2t2z + v2r1y*v2t1x*v2t2z + v2r0x*v2t1y*v2t2z - v2r1x*v2t1y*v2t2z);
VertCoord coeffV2ePower2x = ((v2r0x - v2r1x)*(v2t0z*(v2t1x - v2t2x) + v2t1z*v2t2x - v2t1x*v2t2z + v2t0x*(-v2t1z + v2t2z)))/(-(v2r0y*v2t0z*v2t1x) + v2r1y*v2t0z*v2t1x + v2r0x*v2t0z*v2t1y - v2r1x*v2t0z*v2t1y + v2r0y*v2t0x*v2t1z - v2r1y*v2t0x*v2t1z - v2r0x*v2t0y*v2t1z + v2r1x*v2t0y*v2t1z + v2r0y*v2t0z*v2t2x - v2r1y*v2t0z*v2t2x - v2r0y*v2t1z*v2t2x + v2r1y*v2t1z*v2t2x - v2r0x*v2t0z*v2t2y + v2r1x*v2t0z*v2t2y + v2r0x*v2t1z*v2t2y - v2r1x*v2t1z*v2t2y + v2r0z*(-(v2t0x*v2t1y) + v2t0y*(v2t1x - v2t2x) + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r1z*(-(v2t0y*v2t1x) + v2t0x*v2t1y + v2t0y*v2t2x - v2t1y*v2t2x - v2t0x*v2t2y + v2t1x*v2t2y) - v2r0y*v2t0x*v2t2z + v2r1y*v2t0x*v2t2z + v2r0x*v2t0y*v2t2z - v2r1x*v2t0y*v2t2z + v2r0y*v2t1x*v2t2z - v2r1y*v2t1x*v2t2z - v2r0x*v2t1y*v2t2z + v2r1x*v2t1y*v2t2z);
VertCoord coeffV2ePower3x = ((v2r0x - v2r1x)*(v2t0y*(v2t1x - v2t2x) + v2t1y*v2t2x - v2t1x*v2t2y + v2t0x*(-v2t1y + v2t2y)))/(v2r0y*v2t0z*v2t1x - v2r1y*v2t0z*v2t1x - v2r0x*v2t0z*v2t1y + v2r1x*v2t0z*v2t1y - v2r0y*v2t0x*v2t1z + v2r1y*v2t0x*v2t1z + v2r0x*v2t0y*v2t1z - v2r1x*v2t0y*v2t1z - v2r0y*v2t0z*v2t2x + v2r1y*v2t0z*v2t2x + v2r0y*v2t1z*v2t2x - v2r1y*v2t1z*v2t2x + v2r0x*v2t0z*v2t2y - v2r1x*v2t0z*v2t2y - v2r0x*v2t1z*v2t2y + v2r1x*v2t1z*v2t2y + v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + v2r0y*v2t0x*v2t2z - v2r1y*v2t0x*v2t2z - v2r0x*v2t0y*v2t2z + v2r1x*v2t0y*v2t2z - v2r0y*v2t1x*v2t2z + v2r1y*v2t1x*v2t2z + v2r0x*v2t1y*v2t2z - v2r1x*v2t1y*v2t2z);
VertCoord coeffV2ePower0y = (v2r0z*v2r1y*(-(v2t1y*v2t2x) + v2t0y*(-v2t1x + v2t2x) + v2t0x*(v2t1y - v2t2y) + v2t1x*v2t2y) + v2r0y*(-(v2t0z*v2t1y*v2t2x) + v2t0y*v2t1z*v2t2x + v2t0z*v2t1x*v2t2y - v2t0x*v2t1z*v2t2y + v2r1z*(-(v2t0x*v2t1y) + v2t0y*(v2t1x - v2t2x) + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) - v2t0y*v2t1x*v2t2z + v2t0x*v2t1y*v2t2z + v2r1x*(v2t0z*v2t1y - v2t0y*v2t1z - v2t0z*v2t2y + v2t1z*v2t2y + v2t0y*v2t2z - v2t1y*v2t2z)) + v2r1y*(v2t0z*v2t1y*v2t2x - v2t0y*v2t1z*v2t2x - v2t0z*v2t1x*v2t2y + v2t0x*v2t1z*v2t2y + v2t0y*v2t1x*v2t2z - v2t0x*v2t1y*v2t2z + v2r0x*(v2t0y*v2t1z - v2t1z*v2t2y + v2t0z*(-v2t1y + v2t2y) - v2t0y*v2t2z + v2t1y*v2t2z)))/(v2r0y*v2t0z*v2t1x - v2r1y*v2t0z*v2t1x - v2r0x*v2t0z*v2t1y + v2r1x*v2t0z*v2t1y - v2r0y*v2t0x*v2t1z + v2r1y*v2t0x*v2t1z + v2r0x*v2t0y*v2t1z - v2r1x*v2t0y*v2t1z - v2r0y*v2t0z*v2t2x + v2r1y*v2t0z*v2t2x + v2r0y*v2t1z*v2t2x - v2r1y*v2t1z*v2t2x + v2r0x*v2t0z*v2t2y - v2r1x*v2t0z*v2t2y - v2r0x*v2t1z*v2t2y + v2r1x*v2t1z*v2t2y + v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + v2r0y*v2t0x*v2t2z - v2r1y*v2t0x*v2t2z - v2r0x*v2t0y*v2t2z + v2r1x*v2t0y*v2t2z - v2r0y*v2t1x*v2t2z + v2r1y*v2t1x*v2t2z + v2r0x*v2t1y*v2t2z - v2r1x*v2t1y*v2t2z);
VertCoord coeffV2ePower1y = ((v2r0y - v2r1y)*(v2t0z*(v2t1y - v2t2y) + v2t1z*v2t2y - v2t1y*v2t2z + v2t0y*(-v2t1z + v2t2z)))/(v2r0y*v2t0z*v2t1x - v2r1y*v2t0z*v2t1x - v2r0x*v2t0z*v2t1y + v2r1x*v2t0z*v2t1y - v2r0y*v2t0x*v2t1z + v2r1y*v2t0x*v2t1z + v2r0x*v2t0y*v2t1z - v2r1x*v2t0y*v2t1z - v2r0y*v2t0z*v2t2x + v2r1y*v2t0z*v2t2x + v2r0y*v2t1z*v2t2x - v2r1y*v2t1z*v2t2x + v2r0x*v2t0z*v2t2y - v2r1x*v2t0z*v2t2y - v2r0x*v2t1z*v2t2y + v2r1x*v2t1z*v2t2y + v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + v2r0y*v2t0x*v2t2z - v2r1y*v2t0x*v2t2z - v2r0x*v2t0y*v2t2z + v2r1x*v2t0y*v2t2z - v2r0y*v2t1x*v2t2z + v2r1y*v2t1x*v2t2z + v2r0x*v2t1y*v2t2z - v2r1x*v2t1y*v2t2z);
VertCoord coeffV2ePower2y = (v2r1z*(-(v2t0x*v2t1y) + v2t0y*(v2t1x - v2t2x) + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + (v2r0x - v2r1x)*(v2t0y*v2t1z - v2t1z*v2t2y + v2t0z*(-v2t1y + v2t2y) - v2t0y*v2t2z + v2t1y*v2t2z))/(v2r0y*v2t0z*v2t1x - v2r1y*v2t0z*v2t1x - v2r0x*v2t0z*v2t1y + v2r1x*v2t0z*v2t1y - v2r0y*v2t0x*v2t1z + v2r1y*v2t0x*v2t1z + v2r0x*v2t0y*v2t1z - v2r1x*v2t0y*v2t1z - v2r0y*v2t0z*v2t2x + v2r1y*v2t0z*v2t2x + v2r0y*v2t1z*v2t2x - v2r1y*v2t1z*v2t2x + v2r0x*v2t0z*v2t2y - v2r1x*v2t0z*v2t2y - v2r0x*v2t1z*v2t2y + v2r1x*v2t1z*v2t2y + v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + v2r0y*v2t0x*v2t2z - v2r1y*v2t0x*v2t2z - v2r0x*v2t0y*v2t2z + v2r1x*v2t0y*v2t2z - v2r0y*v2t1x*v2t2z + v2r1y*v2t1x*v2t2z + v2r0x*v2t1y*v2t2z - v2r1x*v2t1y*v2t2z);
VertCoord coeffV2ePower3y = ((v2r0y - v2r1y)*(v2t0y*(v2t1x - v2t2x) + v2t1y*v2t2x - v2t1x*v2t2y + v2t0x*(-v2t1y + v2t2y)))/(v2r0y*v2t0z*v2t1x - v2r1y*v2t0z*v2t1x - v2r0x*v2t0z*v2t1y + v2r1x*v2t0z*v2t1y - v2r0y*v2t0x*v2t1z + v2r1y*v2t0x*v2t1z + v2r0x*v2t0y*v2t1z - v2r1x*v2t0y*v2t1z - v2r0y*v2t0z*v2t2x + v2r1y*v2t0z*v2t2x + v2r0y*v2t1z*v2t2x - v2r1y*v2t1z*v2t2x + v2r0x*v2t0z*v2t2y - v2r1x*v2t0z*v2t2y - v2r0x*v2t1z*v2t2y + v2r1x*v2t1z*v2t2y + v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + v2r0y*v2t0x*v2t2z - v2r1y*v2t0x*v2t2z - v2r0x*v2t0y*v2t2z + v2r1x*v2t0y*v2t2z - v2r0y*v2t1x*v2t2z + v2r1y*v2t1x*v2t2z + v2r0x*v2t1y*v2t2z - v2r1x*v2t1y*v2t2z);
VertCoord coeffV2ePower0z = (v2r0z*(-(v2t0z*v2t1y*v2t2x) + v2t0y*v2t1z*v2t2x + v2t0z*v2t1x*v2t2y - v2t0x*v2t1z*v2t2y - v2t0y*v2t1x*v2t2z + v2t0x*v2t1y*v2t2z + v2r1y*(v2t0x*v2t1z - v2t1z*v2t2x + v2t0z*(-v2t1x + v2t2x) - v2t0x*v2t2z + v2t1x*v2t2z) + v2r1x*(v2t0z*v2t1y - v2t0y*v2t1z - v2t0z*v2t2y + v2t1z*v2t2y + v2t0y*v2t2z - v2t1y*v2t2z)) + v2r1z*(v2t0z*v2t1y*v2t2x - v2t0y*v2t1z*v2t2x - v2t0z*v2t1x*v2t2y + v2t0x*v2t1z*v2t2y + v2t0y*v2t1x*v2t2z - v2t0x*v2t1y*v2t2z + v2r0y*(-(v2t0x*v2t1z) + v2t0z*(v2t1x - v2t2x) + v2t1z*v2t2x + v2t0x*v2t2z - v2t1x*v2t2z) + v2r0x*(-(v2t0z*v2t1y) + v2t0y*v2t1z + v2t0z*v2t2y - v2t1z*v2t2y - v2t0y*v2t2z + v2t1y*v2t2z)))/(v2r0y*v2t0z*v2t1x - v2r1y*v2t0z*v2t1x - v2r0x*v2t0z*v2t1y + v2r1x*v2t0z*v2t1y - v2r0y*v2t0x*v2t1z + v2r1y*v2t0x*v2t1z + v2r0x*v2t0y*v2t1z - v2r1x*v2t0y*v2t1z - v2r0y*v2t0z*v2t2x + v2r1y*v2t0z*v2t2x + v2r0y*v2t1z*v2t2x - v2r1y*v2t1z*v2t2x + v2r0x*v2t0z*v2t2y - v2r1x*v2t0z*v2t2y - v2r0x*v2t1z*v2t2y + v2r1x*v2t1z*v2t2y + v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + v2r0y*v2t0x*v2t2z - v2r1y*v2t0x*v2t2z - v2r0x*v2t0y*v2t2z + v2r1x*v2t0y*v2t2z - v2r0y*v2t1x*v2t2z + v2r1y*v2t1x*v2t2z + v2r0x*v2t1y*v2t2z - v2r1x*v2t1y*v2t2z);
VertCoord coeffV2ePower1z = ((v2r0z - v2r1z)*(v2t0z*(v2t1y - v2t2y) + v2t1z*v2t2y - v2t1y*v2t2z + v2t0y*(-v2t1z + v2t2z)))/(v2r0y*v2t0z*v2t1x - v2r1y*v2t0z*v2t1x - v2r0x*v2t0z*v2t1y + v2r1x*v2t0z*v2t1y - v2r0y*v2t0x*v2t1z + v2r1y*v2t0x*v2t1z + v2r0x*v2t0y*v2t1z - v2r1x*v2t0y*v2t1z - v2r0y*v2t0z*v2t2x + v2r1y*v2t0z*v2t2x + v2r0y*v2t1z*v2t2x - v2r1y*v2t1z*v2t2x + v2r0x*v2t0z*v2t2y - v2r1x*v2t0z*v2t2y - v2r0x*v2t1z*v2t2y + v2r1x*v2t1z*v2t2y + v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + v2r0y*v2t0x*v2t2z - v2r1y*v2t0x*v2t2z - v2r0x*v2t0y*v2t2z + v2r1x*v2t0y*v2t2z - v2r0y*v2t1x*v2t2z + v2r1y*v2t1x*v2t2z + v2r0x*v2t1y*v2t2z - v2r1x*v2t1y*v2t2z);
VertCoord coeffV2ePower2z = ((v2r0z - v2r1z)*(v2t0z*(v2t1x - v2t2x) + v2t1z*v2t2x - v2t1x*v2t2z + v2t0x*(-v2t1z + v2t2z)))/(-(v2r0y*v2t0z*v2t1x) + v2r1y*v2t0z*v2t1x + v2r0x*v2t0z*v2t1y - v2r1x*v2t0z*v2t1y + v2r0y*v2t0x*v2t1z - v2r1y*v2t0x*v2t1z - v2r0x*v2t0y*v2t1z + v2r1x*v2t0y*v2t1z + v2r0y*v2t0z*v2t2x - v2r1y*v2t0z*v2t2x - v2r0y*v2t1z*v2t2x + v2r1y*v2t1z*v2t2x - v2r0x*v2t0z*v2t2y + v2r1x*v2t0z*v2t2y + v2r0x*v2t1z*v2t2y - v2r1x*v2t1z*v2t2y + v2r0z*(-(v2t0x*v2t1y) + v2t0y*(v2t1x - v2t2x) + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r1z*(-(v2t0y*v2t1x) + v2t0x*v2t1y + v2t0y*v2t2x - v2t1y*v2t2x - v2t0x*v2t2y + v2t1x*v2t2y) - v2r0y*v2t0x*v2t2z + v2r1y*v2t0x*v2t2z + v2r0x*v2t0y*v2t2z - v2r1x*v2t0y*v2t2z + v2r0y*v2t1x*v2t2z - v2r1y*v2t1x*v2t2z - v2r0x*v2t1y*v2t2z + v2r1x*v2t1y*v2t2z);
VertCoord coeffV2ePower3z = (v2r0y*(-(v2t0x*v2t1z) + v2t0z*(v2t1x - v2t2x) + v2t1z*v2t2x + v2t0x*v2t2z - v2t1x*v2t2z) + v2r1y*(v2t0x*v2t1z - v2t1z*v2t2x + v2t0z*(-v2t1x + v2t2x) - v2t0x*v2t2z + v2t1x*v2t2z) + (v2r0x - v2r1x)*(v2t0y*v2t1z - v2t1z*v2t2y + v2t0z*(-v2t1y + v2t2y) - v2t0y*v2t2z + v2t1y*v2t2z))/(v2r0y*v2t0z*v2t1x - v2r1y*v2t0z*v2t1x - v2r0x*v2t0z*v2t1y + v2r1x*v2t0z*v2t1y - v2r0y*v2t0x*v2t1z + v2r1y*v2t0x*v2t1z + v2r0x*v2t0y*v2t1z - v2r1x*v2t0y*v2t1z - v2r0y*v2t0z*v2t2x + v2r1y*v2t0z*v2t2x + v2r0y*v2t1z*v2t2x - v2r1y*v2t1z*v2t2x + v2r0x*v2t0z*v2t2y - v2r1x*v2t0z*v2t2y - v2r0x*v2t1z*v2t2y + v2r1x*v2t1z*v2t2y + v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + v2r0y*v2t0x*v2t2z - v2r1y*v2t0x*v2t2z - v2r0x*v2t0y*v2t2z + v2r1x*v2t0y*v2t2z - v2r0y*v2t1x*v2t2z + v2r1y*v2t1x*v2t2z + v2r0x*v2t1y*v2t2z - v2r1x*v2t1y*v2t2z);
/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_1 = coeffV2ePower0z*iv0x - coeffV2ePower0x*iv0z - coeffV2ePower0z*iv1x + iv0z*iv1x + coeffV2ePower0x*iv1z - iv0x*iv1z;
 
if(sgn(ans_1) != 0) return sgn(ans_1);

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_2 = -coeffV2ePower0z + coeffV2ePower1z*iv0x + iv0z - coeffV2ePower1x*iv0z - coeffV2ePower1z*iv1x + coeffV2ePower1x*iv1z;
 
if(sgn(ans_2) != 0) return sgn(ans_2);

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_3 = -coeffV2ePower1z + coeffV2ePower2z*iv0x - coeffV2ePower2x*iv0z - coeffV2ePower2z*iv1x + coeffV2ePower2x*iv1z;
 
if(sgn(ans_3) != 0) return sgn(ans_3);

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_4 = coeffV2ePower0x - coeffV2ePower2z - iv0x + coeffV2ePower3z*iv0x - coeffV2ePower3x*iv0z - coeffV2ePower3z*iv1x + coeffV2ePower3x*iv1z;
 
if(sgn(ans_4) != 0) return sgn(ans_4);

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_5 = coeffV2ePower1x - coeffV2ePower3z;
 
if(sgn(ans_5) != 0) return sgn(ans_5);

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_6 = coeffV2ePower2x;
 
if(sgn(ans_6) != 0) return sgn(ans_6);

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_7 = coeffV2ePower3x;
 
if(sgn(ans_7) != 0) return sgn(ans_7);

/*****************************************************/
/*****************************************************/


 return 0; 

 } 
