#include "sosPredicatesImpl.h"

int SosPredicatesImpl::orient2D_x0_101(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const{


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
cerr << "coeffV2ePower0x : "<<  coeffV2ePower0x<< endl;

VertCoord coeffV2ePower1x = ((v2r0x - v2r1x)*(v2t0z*(v2t1y - v2t2y) + v2t1z*v2t2y - v2t1y*v2t2z + v2t0y*(-v2t1z + v2t2z)))/(-(v2r0y*v2t0z*v2t1x) + v2r1y*v2t0z*v2t1x + v2r0x*v2t0z*v2t1y - v2r1x*v2t0z*v2t1y + v2r0y*v2t0x*v2t1z - v2r1y*v2t0x*v2t1z - v2r0x*v2t0y*v2t1z + v2r1x*v2t0y*v2t1z + v2r0y*v2t0z*v2t2x - v2r1y*v2t0z*v2t2x - v2r0y*v2t1z*v2t2x + v2r1y*v2t1z*v2t2x - v2r0x*v2t0z*v2t2y + v2r1x*v2t0z*v2t2y + v2r0x*v2t1z*v2t2y - v2r1x*v2t1z*v2t2y + v2r0z*(-(v2t0x*v2t1y) + v2t0y*(v2t1x - v2t2x) + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r1z*(-(v2t0y*v2t1x) + v2t0x*v2t1y + v2t0y*v2t2x - v2t1y*v2t2x - v2t0x*v2t2y + v2t1x*v2t2y) - v2r0y*v2t0x*v2t2z + v2r1y*v2t0x*v2t2z + v2r0x*v2t0y*v2t2z - v2r1x*v2t0y*v2t2z + v2r0y*v2t1x*v2t2z - v2r1y*v2t1x*v2t2z - v2r0x*v2t1y*v2t2z + v2r1x*v2t1y*v2t2z);
cerr << "coeffV2ePower1x : "<<  coeffV2ePower1x<< endl;

VertCoord coeffV2ePower2x = ((v2r0x - v2r1x)*(v2t0z*(v2t1x - v2t2x) + v2t1z*v2t2x - v2t1x*v2t2z + v2t0x*(-v2t1z + v2t2z)))/(v2r0y*v2t0z*v2t1x - v2r1y*v2t0z*v2t1x - v2r0x*v2t0z*v2t1y + v2r1x*v2t0z*v2t1y - v2r0y*v2t0x*v2t1z + v2r1y*v2t0x*v2t1z + v2r0x*v2t0y*v2t1z - v2r1x*v2t0y*v2t1z - v2r0y*v2t0z*v2t2x + v2r1y*v2t0z*v2t2x + v2r0y*v2t1z*v2t2x - v2r1y*v2t1z*v2t2x + v2r0x*v2t0z*v2t2y - v2r1x*v2t0z*v2t2y - v2r0x*v2t1z*v2t2y + v2r1x*v2t1z*v2t2y + v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + v2r0y*v2t0x*v2t2z - v2r1y*v2t0x*v2t2z - v2r0x*v2t0y*v2t2z + v2r1x*v2t0y*v2t2z - v2r0y*v2t1x*v2t2z + v2r1y*v2t1x*v2t2z + v2r0x*v2t1y*v2t2z - v2r1x*v2t1y*v2t2z);
cerr << "coeffV2ePower2x : "<<  coeffV2ePower2x<< endl;

VertCoord coeffV2ePower3x = -(((v2r0x - v2r1x)*(v2t0y*(v2t1x - v2t2x) + v2t1y*v2t2x - v2t1x*v2t2y + v2t0x*(-v2t1y + v2t2y)))/(v2r0y*v2t0z*v2t1x - v2r1y*v2t0z*v2t1x - v2r0x*v2t0z*v2t1y + v2r1x*v2t0z*v2t1y - v2r0y*v2t0x*v2t1z + v2r1y*v2t0x*v2t1z + v2r0x*v2t0y*v2t1z - v2r1x*v2t0y*v2t1z - v2r0y*v2t0z*v2t2x + v2r1y*v2t0z*v2t2x + v2r0y*v2t1z*v2t2x - v2r1y*v2t1z*v2t2x + v2r0x*v2t0z*v2t2y - v2r1x*v2t0z*v2t2y - v2r0x*v2t1z*v2t2y + v2r1x*v2t1z*v2t2y + v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + v2r0y*v2t0x*v2t2z - v2r1y*v2t0x*v2t2z - v2r0x*v2t0y*v2t2z + v2r1x*v2t0y*v2t2z - v2r0y*v2t1x*v2t2z + v2r1y*v2t1x*v2t2z + v2r0x*v2t1y*v2t2z - v2r1x*v2t1y*v2t2z));
cerr << "coeffV2ePower3x : "<<  coeffV2ePower3x<< endl;

VertCoord coeffV2ePower0y = (v2r0z*v2r1y*(-(v2t1y*v2t2x) + v2t0y*(-v2t1x + v2t2x) + v2t0x*(v2t1y - v2t2y) + v2t1x*v2t2y) + v2r0y*(-(v2t0z*v2t1y*v2t2x) + v2t0y*v2t1z*v2t2x + v2t0z*v2t1x*v2t2y - v2t0x*v2t1z*v2t2y + v2r1z*(-(v2t0x*v2t1y) + v2t0y*(v2t1x - v2t2x) + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) - v2t0y*v2t1x*v2t2z + v2t0x*v2t1y*v2t2z + v2r1x*(v2t0z*v2t1y - v2t0y*v2t1z - v2t0z*v2t2y + v2t1z*v2t2y + v2t0y*v2t2z - v2t1y*v2t2z)) + v2r1y*(v2t0z*v2t1y*v2t2x - v2t0y*v2t1z*v2t2x - v2t0z*v2t1x*v2t2y + v2t0x*v2t1z*v2t2y + v2t0y*v2t1x*v2t2z - v2t0x*v2t1y*v2t2z + v2r0x*(v2t0y*v2t1z - v2t1z*v2t2y + v2t0z*(-v2t1y + v2t2y) - v2t0y*v2t2z + v2t1y*v2t2z)))/(v2r0y*v2t0z*v2t1x - v2r1y*v2t0z*v2t1x - v2r0x*v2t0z*v2t1y + v2r1x*v2t0z*v2t1y - v2r0y*v2t0x*v2t1z + v2r1y*v2t0x*v2t1z + v2r0x*v2t0y*v2t1z - v2r1x*v2t0y*v2t1z - v2r0y*v2t0z*v2t2x + v2r1y*v2t0z*v2t2x + v2r0y*v2t1z*v2t2x - v2r1y*v2t1z*v2t2x + v2r0x*v2t0z*v2t2y - v2r1x*v2t0z*v2t2y - v2r0x*v2t1z*v2t2y + v2r1x*v2t1z*v2t2y + v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + v2r0y*v2t0x*v2t2z - v2r1y*v2t0x*v2t2z - v2r0x*v2t0y*v2t2z + v2r1x*v2t0y*v2t2z - v2r0y*v2t1x*v2t2z + v2r1y*v2t1x*v2t2z + v2r0x*v2t1y*v2t2z - v2r1x*v2t1y*v2t2z);
cerr << "coeffV2ePower0y : "<<  coeffV2ePower0y<< endl;

VertCoord coeffV2ePower1y = -(((v2r0y - v2r1y)*(v2t0z*(v2t1y - v2t2y) + v2t1z*v2t2y - v2t1y*v2t2z + v2t0y*(-v2t1z + v2t2z)))/(v2r0y*v2t0z*v2t1x - v2r1y*v2t0z*v2t1x - v2r0x*v2t0z*v2t1y + v2r1x*v2t0z*v2t1y - v2r0y*v2t0x*v2t1z + v2r1y*v2t0x*v2t1z + v2r0x*v2t0y*v2t1z - v2r1x*v2t0y*v2t1z - v2r0y*v2t0z*v2t2x + v2r1y*v2t0z*v2t2x + v2r0y*v2t1z*v2t2x - v2r1y*v2t1z*v2t2x + v2r0x*v2t0z*v2t2y - v2r1x*v2t0z*v2t2y - v2r0x*v2t1z*v2t2y + v2r1x*v2t1z*v2t2y + v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + v2r0y*v2t0x*v2t2z - v2r1y*v2t0x*v2t2z - v2r0x*v2t0y*v2t2z + v2r1x*v2t0y*v2t2z - v2r0y*v2t1x*v2t2z + v2r1y*v2t1x*v2t2z + v2r0x*v2t1y*v2t2z - v2r1x*v2t1y*v2t2z));
cerr << "coeffV2ePower1y : "<<  coeffV2ePower1y<< endl;

VertCoord coeffV2ePower2y = ((v2r0y - v2r1y)*(v2t0z*(v2t1x - v2t2x) + v2t1z*v2t2x - v2t1x*v2t2z + v2t0x*(-v2t1z + v2t2z)))/(v2r0y*v2t0z*v2t1x - v2r1y*v2t0z*v2t1x - v2r0x*v2t0z*v2t1y + v2r1x*v2t0z*v2t1y - v2r0y*v2t0x*v2t1z + v2r1y*v2t0x*v2t1z + v2r0x*v2t0y*v2t1z - v2r1x*v2t0y*v2t1z - v2r0y*v2t0z*v2t2x + v2r1y*v2t0z*v2t2x + v2r0y*v2t1z*v2t2x - v2r1y*v2t1z*v2t2x + v2r0x*v2t0z*v2t2y - v2r1x*v2t0z*v2t2y - v2r0x*v2t1z*v2t2y + v2r1x*v2t1z*v2t2y + v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + v2r0y*v2t0x*v2t2z - v2r1y*v2t0x*v2t2z - v2r0x*v2t0y*v2t2z + v2r1x*v2t0y*v2t2z - v2r0y*v2t1x*v2t2z + v2r1y*v2t1x*v2t2z + v2r0x*v2t1y*v2t2z - v2r1x*v2t1y*v2t2z);
cerr << "coeffV2ePower2y : "<<  coeffV2ePower2y<< endl;

VertCoord coeffV2ePower3y = -(((v2r0y - v2r1y)*(v2t0y*(v2t1x - v2t2x) + v2t1y*v2t2x - v2t1x*v2t2y + v2t0x*(-v2t1y + v2t2y)))/(v2r0y*v2t0z*v2t1x - v2r1y*v2t0z*v2t1x - v2r0x*v2t0z*v2t1y + v2r1x*v2t0z*v2t1y - v2r0y*v2t0x*v2t1z + v2r1y*v2t0x*v2t1z + v2r0x*v2t0y*v2t1z - v2r1x*v2t0y*v2t1z - v2r0y*v2t0z*v2t2x + v2r1y*v2t0z*v2t2x + v2r0y*v2t1z*v2t2x - v2r1y*v2t1z*v2t2x + v2r0x*v2t0z*v2t2y - v2r1x*v2t0z*v2t2y - v2r0x*v2t1z*v2t2y + v2r1x*v2t1z*v2t2y + v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + v2r0y*v2t0x*v2t2z - v2r1y*v2t0x*v2t2z - v2r0x*v2t0y*v2t2z + v2r1x*v2t0y*v2t2z - v2r0y*v2t1x*v2t2z + v2r1y*v2t1x*v2t2z + v2r0x*v2t1y*v2t2z - v2r1x*v2t1y*v2t2z));
cerr << "coeffV2ePower3y : "<<  coeffV2ePower3y<< endl;

VertCoord coeffV2ePower0z = (v2r0z*(-(v2t0z*v2t1y*v2t2x) + v2t0y*v2t1z*v2t2x + v2t0z*v2t1x*v2t2y - v2t0x*v2t1z*v2t2y - v2t0y*v2t1x*v2t2z + v2t0x*v2t1y*v2t2z + v2r1y*(v2t0x*v2t1z - v2t1z*v2t2x + v2t0z*(-v2t1x + v2t2x) - v2t0x*v2t2z + v2t1x*v2t2z) + v2r1x*(v2t0z*v2t1y - v2t0y*v2t1z - v2t0z*v2t2y + v2t1z*v2t2y + v2t0y*v2t2z - v2t1y*v2t2z)) + v2r1z*(v2t0z*v2t1y*v2t2x - v2t0y*v2t1z*v2t2x - v2t0z*v2t1x*v2t2y + v2t0x*v2t1z*v2t2y + v2t0y*v2t1x*v2t2z - v2t0x*v2t1y*v2t2z + v2r0y*(-(v2t0x*v2t1z) + v2t0z*(v2t1x - v2t2x) + v2t1z*v2t2x + v2t0x*v2t2z - v2t1x*v2t2z) + v2r0x*(-(v2t0z*v2t1y) + v2t0y*v2t1z + v2t0z*v2t2y - v2t1z*v2t2y - v2t0y*v2t2z + v2t1y*v2t2z)))/(v2r0y*v2t0z*v2t1x - v2r1y*v2t0z*v2t1x - v2r0x*v2t0z*v2t1y + v2r1x*v2t0z*v2t1y - v2r0y*v2t0x*v2t1z + v2r1y*v2t0x*v2t1z + v2r0x*v2t0y*v2t1z - v2r1x*v2t0y*v2t1z - v2r0y*v2t0z*v2t2x + v2r1y*v2t0z*v2t2x + v2r0y*v2t1z*v2t2x - v2r1y*v2t1z*v2t2x + v2r0x*v2t0z*v2t2y - v2r1x*v2t0z*v2t2y - v2r0x*v2t1z*v2t2y + v2r1x*v2t1z*v2t2y + v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + v2r0y*v2t0x*v2t2z - v2r1y*v2t0x*v2t2z - v2r0x*v2t0y*v2t2z + v2r1x*v2t0y*v2t2z - v2r0y*v2t1x*v2t2z + v2r1y*v2t1x*v2t2z + v2r0x*v2t1y*v2t2z - v2r1x*v2t1y*v2t2z);
cerr << "coeffV2ePower0z : "<<  coeffV2ePower0z<< endl;

VertCoord coeffV2ePower1z = ((v2r0z - v2r1z)*(v2t0z*(v2t1y - v2t2y) + v2t1z*v2t2y - v2t1y*v2t2z + v2t0y*(-v2t1z + v2t2z)))/(-(v2r0y*v2t0z*v2t1x) + v2r1y*v2t0z*v2t1x + v2r0x*v2t0z*v2t1y - v2r1x*v2t0z*v2t1y + v2r0y*v2t0x*v2t1z - v2r1y*v2t0x*v2t1z - v2r0x*v2t0y*v2t1z + v2r1x*v2t0y*v2t1z + v2r0y*v2t0z*v2t2x - v2r1y*v2t0z*v2t2x - v2r0y*v2t1z*v2t2x + v2r1y*v2t1z*v2t2x - v2r0x*v2t0z*v2t2y + v2r1x*v2t0z*v2t2y + v2r0x*v2t1z*v2t2y - v2r1x*v2t1z*v2t2y + v2r0z*(-(v2t0x*v2t1y) + v2t0y*(v2t1x - v2t2x) + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r1z*(-(v2t0y*v2t1x) + v2t0x*v2t1y + v2t0y*v2t2x - v2t1y*v2t2x - v2t0x*v2t2y + v2t1x*v2t2y) - v2r0y*v2t0x*v2t2z + v2r1y*v2t0x*v2t2z + v2r0x*v2t0y*v2t2z - v2r1x*v2t0y*v2t2z + v2r0y*v2t1x*v2t2z - v2r1y*v2t1x*v2t2z - v2r0x*v2t1y*v2t2z + v2r1x*v2t1y*v2t2z);
cerr << "coeffV2ePower1z : "<<  coeffV2ePower1z<< endl;

VertCoord coeffV2ePower2z = ((v2r0z - v2r1z)*(v2t0z*(v2t1x - v2t2x) + v2t1z*v2t2x - v2t1x*v2t2z + v2t0x*(-v2t1z + v2t2z)))/(v2r0y*v2t0z*v2t1x - v2r1y*v2t0z*v2t1x - v2r0x*v2t0z*v2t1y + v2r1x*v2t0z*v2t1y - v2r0y*v2t0x*v2t1z + v2r1y*v2t0x*v2t1z + v2r0x*v2t0y*v2t1z - v2r1x*v2t0y*v2t1z - v2r0y*v2t0z*v2t2x + v2r1y*v2t0z*v2t2x + v2r0y*v2t1z*v2t2x - v2r1y*v2t1z*v2t2x + v2r0x*v2t0z*v2t2y - v2r1x*v2t0z*v2t2y - v2r0x*v2t1z*v2t2y + v2r1x*v2t1z*v2t2y + v2r1z*(v2t0y*v2t1x - v2t0x*v2t1y - v2t0y*v2t2x + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r0z*(v2t0x*v2t1y - v2t1y*v2t2x + v2t0y*(-v2t1x + v2t2x) - v2t0x*v2t2y + v2t1x*v2t2y) + v2r0y*v2t0x*v2t2z - v2r1y*v2t0x*v2t2z - v2r0x*v2t0y*v2t2z + v2r1x*v2t0y*v2t2z - v2r0y*v2t1x*v2t2z + v2r1y*v2t1x*v2t2z + v2r0x*v2t1y*v2t2z - v2r1x*v2t1y*v2t2z);
cerr << "coeffV2ePower2z : "<<  coeffV2ePower2z<< endl;

VertCoord coeffV2ePower3z = ((v2r0z - v2r1z)*(v2t0y*(v2t1x - v2t2x) + v2t1y*v2t2x - v2t1x*v2t2y + v2t0x*(-v2t1y + v2t2y)))/(-(v2r0y*v2t0z*v2t1x) + v2r1y*v2t0z*v2t1x + v2r0x*v2t0z*v2t1y - v2r1x*v2t0z*v2t1y + v2r0y*v2t0x*v2t1z - v2r1y*v2t0x*v2t1z - v2r0x*v2t0y*v2t1z + v2r1x*v2t0y*v2t1z + v2r0y*v2t0z*v2t2x - v2r1y*v2t0z*v2t2x - v2r0y*v2t1z*v2t2x + v2r1y*v2t1z*v2t2x - v2r0x*v2t0z*v2t2y + v2r1x*v2t0z*v2t2y + v2r0x*v2t1z*v2t2y - v2r1x*v2t1z*v2t2y + v2r0z*(-(v2t0x*v2t1y) + v2t0y*(v2t1x - v2t2x) + v2t1y*v2t2x + v2t0x*v2t2y - v2t1x*v2t2y) + v2r1z*(-(v2t0y*v2t1x) + v2t0x*v2t1y + v2t0y*v2t2x - v2t1y*v2t2x - v2t0x*v2t2y + v2t1x*v2t2y) - v2r0y*v2t0x*v2t2z + v2r1y*v2t0x*v2t2z + v2r0x*v2t0y*v2t2z - v2r1x*v2t0y*v2t2z + v2r0y*v2t1x*v2t2z - v2r1y*v2t1x*v2t2z - v2r0x*v2t1y*v2t2z + v2r1x*v2t1y*v2t2z);
cerr << "coeffV2ePower3z : "<<  coeffV2ePower3z<< endl;

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_1 = -(coeffV2ePower0z*iv0y) + coeffV2ePower0y*iv0z + coeffV2ePower0z*iv1y - iv0z*iv1y - coeffV2ePower0y*iv1z + iv0y*iv1z;
cerr << "ans_1 : "<<  ans_1<< endl;

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_2 = -(coeffV2ePower1z*iv0y) + coeffV2ePower1y*iv0z + coeffV2ePower1z*iv1y - coeffV2ePower1y*iv1z;
cerr << "ans_2 : "<<  ans_2<< endl;

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_3 = -coeffV2ePower0z - coeffV2ePower2z*iv0y + coeffV2ePower2y*iv0z + coeffV2ePower2z*iv1y + iv1z - coeffV2ePower2y*iv1z;
cerr << "ans_3 : "<<  ans_3<< endl;

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_4 = coeffV2ePower0y - coeffV2ePower1z - coeffV2ePower3z*iv0y + coeffV2ePower3y*iv0z - iv1y + coeffV2ePower3z*iv1y - coeffV2ePower3y*iv1z;
cerr << "ans_4 : "<<  ans_4<< endl;

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_5 = coeffV2ePower1y - coeffV2ePower2z;
cerr << "ans_5 : "<<  ans_5<< endl;

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_6 = coeffV2ePower2y - coeffV2ePower3z;
cerr << "ans_6 : "<<  ans_6<< endl;

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_7 = coeffV2ePower3y;
cerr << "ans_7 : "<<  ans_7<< endl;

/*****************************************************/
/*****************************************************/


 return 0; 

 } 
