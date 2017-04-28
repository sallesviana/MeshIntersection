#include "sosPredicatesImpl.h"

int SosPredicatesImpl::orient1D_x_01(const InputVertex &v0, const VertexFromIntersection &v1) const{


/*****************************************************/
const VertCoord &iv0x = getCoordinates(v0)[0];
const VertCoord &iv0y = getCoordinates(v0)[1];
const VertCoord &iv0z = getCoordinates(v0)[2];





/*****************************************************/
VertCoord coeffV0ePower0x = getEpsCoefficientsVertexFromIntersection(v0, 0, 0);
VertCoord coeffV0ePower1x = getEpsCoefficientsVertexFromIntersection(v0, 1, 0);
VertCoord coeffV0ePower2x = getEpsCoefficientsVertexFromIntersection(v0, 2, 0);
VertCoord coeffV0ePower3x = getEpsCoefficientsVertexFromIntersection(v0, 3, 0);
VertCoord coeffV0ePower0y = getEpsCoefficientsVertexFromIntersection(v0, 0, 1);
VertCoord coeffV0ePower1y = getEpsCoefficientsVertexFromIntersection(v0, 1, 1);
VertCoord coeffV0ePower2y = getEpsCoefficientsVertexFromIntersection(v0, 2, 1);
VertCoord coeffV0ePower3y = getEpsCoefficientsVertexFromIntersection(v0, 3, 1);
VertCoord coeffV0ePower0z = getEpsCoefficientsVertexFromIntersection(v0, 0, 2);
VertCoord coeffV0ePower1z = getEpsCoefficientsVertexFromIntersection(v0, 1, 2);
VertCoord coeffV0ePower2z = getEpsCoefficientsVertexFromIntersection(v0, 2, 2);
VertCoord coeffV0ePower3z = getEpsCoefficientsVertexFromIntersection(v0, 3, 2);
/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_1 = coeffV1ePower0x - iv0x;
 
if(sgn(ans_1) != 0) return sgn(ans_1);

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_2 = coeffV1ePower1x;
 
if(sgn(ans_2) != 0) return sgn(ans_2);

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_3 = coeffV1ePower2x;
 
if(sgn(ans_3) != 0) return sgn(ans_3);

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_4 = coeffV1ePower3x;
 
if(sgn(ans_4) != 0) return sgn(ans_4);

/*****************************************************/
/*****************************************************/


 return 0; 

 } 
