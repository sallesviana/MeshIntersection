#include "sosPredicatesImpl.h"

int SosPredicatesImpl::orient2D_z0_110(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const{


/*****************************************************/
const VertCoord &iv0x = getCoordinates(v0)[0];
const VertCoord &iv0y = getCoordinates(v0)[1];
const VertCoord &iv0z = getCoordinates(v0)[2];

const VertCoord &iv1x = getCoordinates(v1)[0];
const VertCoord &iv1y = getCoordinates(v1)[1];
const VertCoord &iv1z = getCoordinates(v1)[2];

const VertCoord &iv2x = getCoordinates(v2)[0];
const VertCoord &iv2y = getCoordinates(v2)[1];
const VertCoord &iv2z = getCoordinates(v2)[2];

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_1 = -(iv0y*(-iv0x + iv1x)) - iv0x*(iv0y - iv1y) + (iv0y - iv1y)*iv2x + (-iv0x + iv1x)*iv2y;
cerr << "ans_1 : "<<  ans_1<< endl;

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_2 = -iv0y + iv1y;
cerr << "ans_2 : "<<  ans_2<< endl;

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_3 = iv0x - iv1x;
cerr << "ans_3 : "<<  ans_3<< endl;

/*****************************************************/
/*****************************************************/


 return 0; 

 } 
