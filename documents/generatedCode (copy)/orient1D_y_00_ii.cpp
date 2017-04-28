#include "sosPredicatesImpl.h"

int SosPredicatesImpl::orient1D_y_00(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const{


/*****************************************************/








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


VertCoord coeffV1ePower0x = getEpsCoefficientsVertexFromIntersection(v1, 0, 0);
VertCoord coeffV1ePower1x = getEpsCoefficientsVertexFromIntersection(v1, 1, 0);
VertCoord coeffV1ePower2x = getEpsCoefficientsVertexFromIntersection(v1, 2, 0);
VertCoord coeffV1ePower3x = getEpsCoefficientsVertexFromIntersection(v1, 3, 0);
VertCoord coeffV1ePower0y = getEpsCoefficientsVertexFromIntersection(v1, 0, 1);
VertCoord coeffV1ePower1y = getEpsCoefficientsVertexFromIntersection(v1, 1, 1);
VertCoord coeffV1ePower2y = getEpsCoefficientsVertexFromIntersection(v1, 2, 1);
VertCoord coeffV1ePower3y = getEpsCoefficientsVertexFromIntersection(v1, 3, 1);
VertCoord coeffV1ePower0z = getEpsCoefficientsVertexFromIntersection(v1, 0, 2);
VertCoord coeffV1ePower1z = getEpsCoefficientsVertexFromIntersection(v1, 1, 2);
VertCoord coeffV1ePower2z = getEpsCoefficientsVertexFromIntersection(v1, 2, 2);
VertCoord coeffV1ePower3z = getEpsCoefficientsVertexFromIntersection(v1, 3, 2);
/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_1 = -coeffV0ePower0y + coeffV1ePower0y;
 
if(sgn(ans_1) != 0) return sgn(ans_1);

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_2 = -coeffV0ePower1y + coeffV1ePower1y;
 
if(sgn(ans_2) != 0) return sgn(ans_2);

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_3 = -coeffV0ePower2y + coeffV1ePower2y;
 
if(sgn(ans_3) != 0) return sgn(ans_3);

/*****************************************************/
/*****************************************************/
/*****************************************************/
/*****************************************************/
VertCoord ans_4 = -coeffV0ePower3y + coeffV1ePower3y;
 
if(sgn(ans_4) != 0) return sgn(ans_4);

/*****************************************************/
/*****************************************************/


 return 0; 

 } 
