#include "sosPredicatesImpl.h"

int SosPredicatesImpl::orient1D_x_11(const InputVertex &v0, const InputVertex &v1) const{


/*****************************************************/
const VertCoord &iv0x = getCoordinates(v0)[0];
const VertCoord &iv0y = getCoordinates(v0)[1];
const VertCoord &iv0z = getCoordinates(v0)[2];

const VertCoord &iv1x = getCoordinates(v1)[0];
const VertCoord &iv1y = getCoordinates(v1)[1];
const VertCoord &iv1z = getCoordinates(v1)[2];

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_1 = -iv0x + iv1x;
 
if(sgn(ans_1) != 0) return sgn(ans_1);

/*****************************************************/
/*****************************************************/


 return 0; 

 } 
