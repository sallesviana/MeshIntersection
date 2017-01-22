/* Triangle/triangle intersection test routine,
 * by Tomas Moller, 1997.
 * See article "A Fast Triangle-Triangle Intersection Test",
 * Journal of Graphics Tools, 2(2), 1997
 * updated: 2001-06-20 (added line of intersection)
 *
 * int tri_tri_intersect(VertCoord V0[3],VertCoord V1[3],VertCoord V2[3],
 *                       VertCoord U0[3],VertCoord U1[3],VertCoord U2[3])
 *
 * parameters: vertices of triangle 1: V0,V1,V2
 *             vertices of triangle 2: U0,U1,U2
 * result    : returns 1 if the triangles intersect, otherwise 0
 *
 * Here is a version withouts divisions (a little faster)
 * int NoDivTriTriIsect(VertCoord V0[3],VertCoord V1[3],VertCoord V2[3],
 *                      VertCoord U0[3],VertCoord U1[3],VertCoord U2[3]);
 * 
 * This version computes the line of intersection as well (if they are not coplanar):
 * int tri_tri_intersect_with_isectline(VertCoord V0[3],VertCoord V1[3],VertCoord V2[3], 
 *				        VertCoord U0[3],VertCoord U1[3],VertCoord U2[3],int *coplanar,
 *				        VertCoord isectpt1[3],VertCoord isectpt2[3]);
 * coplanar returns whether the tris are coplanar
 * isectpt1, isectpt2 are the endpoints of the line of intersection
 */

#include <math.h>

#define FABS(x) ((VertCoord)fabs(x))        /* implement as is fastest on your machine */




/* some macros */
/*   
#define CROSS(dest,v1,v2)                      \
              dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
              dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
              dest[2]=v1[0]*v2[1]-v1[1]*v2[0];  */

#define CROSS(dest,v1,v2,temp)                      \
              dest[0]=v1[1]; \
              dest[0]*=v2[2]; \
              dest[2] =v1[2]; \
              dest[2] *= v2[1]; \
              dest[0] -= dest[2]; \
              dest[1] = v1[2]; \
              dest[1] *= v2[0]; \
              dest[2] = v1[0]; \
              dest[2] *= v2[2]; \
              dest[1] -= dest[2]; \
              dest[2] =  v1[0]; \
              dest[2] *= v2[1]; \
              temp = v1[1]; \
              temp *= v2[0]; \
              dest[2] -= temp;


//#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

void DOT(VertCoord &dest, const Point &v1, const Point &v2, VertCoord &tmp) {
  dest = v1[0];
  dest *= v2[0];
  tmp = v1[1];
  tmp *= v2[1];
  dest += tmp;
  tmp = v1[2];
  tmp *= v2[2];
  dest += tmp;  
}

void MinusDOT(VertCoord &dest, const Point &v1, const Point &v2, VertCoord &tmp) {
  dest = -v1[0];
  dest *= v2[0];
  tmp = v1[1];
  tmp *= v2[1];
  dest -= tmp;
  tmp = v1[2];
  tmp *= v2[2];
  dest -= tmp;  
}


//#define SUB(dest,v1,v2) dest[0]=v1[0]-v2[0]; dest[1]=v1[1]-v2[1]; dest[2]=v1[2]-v2[2]; 

//#define ADD(dest,v1,v2) dest[0]=v1[0]+v2[0]; dest[1]=v1[1]+v2[1]; dest[2]=v1[2]+v2[2]; 

//#define MULT(dest,v,factor) dest[0]=factor*v[0]; dest[1]=factor*v[1]; dest[2]=factor*v[2]; 




#define SUB(dest,v1,v2) dest[0]=v1[0]; dest[0]-=v2[0]; dest[1]=v1[1]; dest[1]-=v2[1]; dest[2]=v1[2]; dest[2]-=v2[2]; 

#define ADD(dest,v1,v2) dest[0]=v1[0]; dest[0]+=v2[0]; dest[1]=v1[1]; dest[1]+=v2[1]; dest[2]=v1[2]; dest[2]+=v2[2]; 

//#define MULT(dest,v,factor) dest[0]=factor; dest[0]*=v[0]; dest[1]=factor; dest[1]*=v[1]; dest[2]=factor; dest[2]*=v[2];


#define SET(dest,src) dest[0]=src[0]; dest[1]=src[1]; dest[2]=src[2]; 

void fabs(rational &dest, const rational &r) {
  dest = r;
  if (dest<0)
    dest *= -1;
}

/* sort so that a<=b */
#define SORT(a,b,temp)       \
             if(a>b)    \
             {          \
               temp=a;     \
               a=b;     \
               b=temp;     \
             }


/* this edge to edge test is based on Franlin Antonio's gem:
   "Faster Line Segment Intersection", in Graphics Gems III,
   pp. 199-202 */ 
/*#define EDGE_EDGE_TEST(V0,U0,U1)                      \
  Bx=U0[i0]-U1[i0];                                   \
  By=U0[i1]-U1[i1];                                   \
  Cx=V0[i0]-U0[i0];                                   \
  Cy=V0[i1]-U0[i1];                                   \
  f=Ay*Bx-Ax*By;                                      \
  d=By*Cx-Bx*Cy;                                      \
  if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))  \
  {                                                   \
    e=Ax*Cy-Ay*Cx;                                    \
    if(f>0)                                           \
    {                                                 \
      if(e>=0 && e<=f) return 1;                      \
    }                                                 \
    else                                              \
    {                                                 \
      if(e<=0 && e>=f) return 1;                      \
    }                                                 \
  }      */                          
#define EDGE_EDGE_TEST(V0,U0,U1)                      \
  Bx=U0[i0];                                    \
  Bx-=U1[i0];                                   \
  By=U0[i1];                                    \
  By-=U1[i1];                                   \
  Cx=V0[i0];                                    \
  Cx-=U0[i0];                                   \
  Cy=V0[i1];                                    \
  Cy-=U0[i1];                                   \
  f=Ay;                                   \
  f*=Bx;                                   \
  tmp = Ax;                                   \
  tmp *= By;                                   \
  f-= tmp;                                   \
  d = By;                                   \
  d*= Cx;                                   \
  tmp  = Bx;                                   \
  tmp *=Cy;                                   \
  d-= tmp;                                   \
  if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))  \
  {                                                   \
    e=Ax;                                   \
    e*=Cy;                                   \
    tmp = Ay;                                   \
    tmp *= Cx;                                   \
    e-= tmp;                                   \
    if(f>0)                                           \
    {                                                 \
      if(e>=0 && e<=f) return 1;                      \
    }                                                 \
    else                                              \
    {                                                 \
      if(e<=0 && e>=f) return 1;                      \
    }                                                 \
  }        


#define EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2,temp) \
{                                              \
  VertCoord &Ax=temp[0],&Ay=temp[1],&Bx=temp[2],&By=temp[3],&Cx=temp[4],&Cy=temp[5],&e=temp[6],&d=temp[7],&f=temp[8],&tmp=temp[9];               \
  Ax=V1[i0];  \
  Ax-=V0[i0];                            \
  Ay=V1[i1];    \
  Ay-=V0[i1];                            \
  /* test edge U0,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U0,U1);                    \
  /* test edge U1,U2 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U1,U2);                    \
  /* test edge U2,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U2,U0);                    \
}
/*
#define POINT_IN_TRI(V0,U0,U1,U2)           \
{                                           \
  VertCoord a,b,c,d0,d1,d2;                     \
  // is T1 completly inside T2?           \
  // check if V0 is inside tri(U0,U1,U2)  \
  a=U1[i1]-U0[i1];                          \
  b=-(U1[i0]-U0[i0]);                       \
  c=-a*U0[i0]-b*U0[i1];                     \
  d0=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U2[i1]-U1[i1];                          \
  b=-(U2[i0]-U1[i0]);                       \
  c=-a*U1[i0]-b*U1[i1];                     \
  d1=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U0[i1]-U2[i1];                          \
  b=-(U0[i0]-U2[i0]);                       \
  c=-a*U2[i0]-b*U2[i1];                     \
  d2=a*V0[i0]+b*V0[i1]+c;                   \
  if(d0*d1>0.0)                             \
  {                                         \
    if(d0*d2>0.0) return 1;                 \
  }                                         \
}*/

#define POINT_IN_TRI(V0,U0,U1,U2, tempRationals)       \
{                                           \
  VertCoord &a = tempRationals[0]; \
  VertCoord &b = tempRationals[1],&c=tempRationals[2],&d0=tempRationals[3],&d1=tempRationals[4],&d2=tempRationals[5],&tmp=tempRationals[6]; \
    /* is T1 completly inside T2? */          \
  /* check if V0 is inside tri(U0,U1,U2) */ \
  a=U1[i1]; \
  a-=U0[i1];                          \
  b=0; \
  b-=U1[i0]; \
  b+=U0[i0];                       \
  c = a; \
  c *= U0[i0]; \
  c *= -1; \
  d1 = b; \
  d1 *= U0[i1]; \
  c -= d1; \
            \
  d0 = a; \
  d0 *= V0[i0]; \
  d1 = b; \
  d1 *= V0[i1]; \
  d0 += d1; \
  d0 += c; \
            \
  a = U2[i1]; \
  a -= U1[i1]; \
  b = 0; \
  b -= U2[i0]; \
  b += U1[i0]; \
              \
  c = a; \
  c *= -1; \
  c *= U1[i0]; \
  tmp = b; \
  tmp *= U1[i1]; \
  c -= tmp; \
  d1 = a; \
  d1 *= V0[i0]; \
  tmp = b;   \
  tmp *= V0[i1]; \
  d1+= tmp; \
  d1 += c; \
               \
  a = U0[i1]; \
  a -= U2[i1]; \
  b = 0; \
  b -= U0[i0]; \
  b += U2[i0]; \
  c = 0; \
  c -= a; \
  c *= U2[i0]; \
  tmp = b; \
  tmp *= U2[i1]; \
  c -= tmp;   \
              \
  d2 = a; \
  d2 *= V0[i0]; \
  tmp= b; \
  tmp *=V0[i1]; \
  d2 += tmp;  \
  d2 += c;    \
             \
  tmp = d0; \
  tmp *= d1; \
  a = d0; \
  a *= d2; \
  if(tmp>0)                             \
  {                                         \
    if(a>0) return 1;                 \
  }                                         \
}

//Temp should hold at least 3 VertCoords
int coplanar_tri_tri(VertCoord N[3],VertCoord *V0,VertCoord *V1,VertCoord *V2,
                     VertCoord *U0,VertCoord *U1,VertCoord *U2,VertCoord *Temp)
{
   short i0,i1;
   /* first project onto an axis-aligned plane, that maximizes the area */
   /* of the triangles, compute indices: i0,i1. */
   //A[0]=fabs(N[0]);
   //A[1]=fabs(N[1]);
   //A[2]=fabs(N[2]);
   fabs(Temp[0],N[0]);
   fabs(Temp[1],N[1]);
   fabs(Temp[2],N[2]);
   if(Temp[0]>Temp[1])
   {
      if(Temp[0]>Temp[2])  
      {
          i0=1;      /* A[0] is greatest */
          i1=2;
      }
      else
      {
          i0=0;      /* A[2] is greatest */
          i1=1;
      }
   }
   else   /* A[0]<=A[1] */
   {
      if(Temp[2]>Temp[1])
      {
          i0=0;      /* A[2] is greatest */
          i1=1;                                           
      }
      else
      {
          i0=0;      /* A[1] is greatest */
          i1=2;
      }
    }            
    VertCoord *tempToPointInTri = Temp+3;   
    /* test all edges of triangle 1 against the edges of triangle 2 */
    EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2,tempToPointInTri);
    EDGE_AGAINST_TRI_EDGES(V1,V2,U0,U1,U2,tempToPointInTri);
    EDGE_AGAINST_TRI_EDGES(V2,V0,U0,U1,U2,tempToPointInTri);
                
    /* finally, test if tri1 is totally contained in tri2 or vice versa */
    
    POINT_IN_TRI(V0,U0,U1,U2,tempToPointInTri);
    POINT_IN_TRI(U0,V0,V1,V2,tempToPointInTri);

    return 0;
}



/* sort so that a<=b */
#define SORT2(a,b,smallest,temp)       \
             if(a>b)       \
             {             \
               temp=a;        \
               a=b;        \
               b=temp;        \
               smallest=1; \
             }             \
             else smallest=0;

//tmpVars should have at least 5 rationals
inline void isect2(Point &VTX0,Point & VTX1,Point & VTX2,const VertCoord &VV0, const VertCoord &VV1,const VertCoord &VV2,
	    const VertCoord &D0, const VertCoord &D1, const VertCoord &D2,VertCoord &isect0,VertCoord  &isect1,Point &isectpoint0,Point &isectpoint1, VertCoord *tmpVars) 
{
  //VertCoord tmp=D0/(D0-D1);          
  VertCoord &tmp = tmpVars[0];
  VertCoord &tmp2 = tmpVars[4];
  tmp = D0;
  tmp2 = D0;
  tmp2 -= D1;
  tmp /= tmp2;

  //VertCoord diff[3];
  VertCoord *diff = tmpVars+1;
  //*isect0=VV0+(VV1-VV0)*tmp;         
  isect0   =VV1; //+(VV1-VV0)*tmp; 
  isect0  -=VV0;
  isect0 *= tmp;
  isect0 += VV0;

  SUB(diff,VTX1,VTX0);              
  //MULT(diff,diff,tmp);     
  diff[0]*=tmp;
  diff[1]*=tmp;
  diff[2]*=tmp;

  ADD(isectpoint0,diff,VTX0);        
  //tmp=D0/(D0-D2);     
  tmp = D0;
  tmp2 = D0;
  tmp2 -= D2;
  tmp /= tmp2;

  //*isect1=VV0+(VV2-VV0)*tmp;       
  isect1 = VV2;
  isect1 -= VV0;
  isect1 *= tmp;
  isect1 += VV0;
  SUB(diff,VTX2,VTX0);                   
  //MULT(diff,diff,tmp);      
  diff[0]*=tmp;
  diff[1]*=tmp;
  diff[2]*=tmp;         
  ADD(isectpoint1,VTX0,diff);          
}




// VV0, VV1 and VV2 are the coordinates of Vert0,Vert1,... considering the largest coordinate (projected..)
inline int compute_intervals_isectline(Point &VERT0,Point & VERT1,Point &VERT2,
				       const VertCoord &VV0,const VertCoord &VV1,const VertCoord &VV2,const VertCoord &D0,const VertCoord &D1,const VertCoord &D2,
				       const int D0D1,const int D0D2,VertCoord &isect0,VertCoord & isect1,
				       Point &isectpoint0,Point & isectpoint1, pair<int,int> &edgeCreatedA0, pair<int,int> &edgeCreatedA1, VertCoord * tmpVars)
{

  //TODO: SoS when triangles touch...

  if(D0D1>0.0f)                                        
  {                                                    
    /* here we know that D0D2<=0.0 */                  
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */
    isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,isect0,isect1,isectpoint0,isectpoint1,tmpVars);

    edgeCreatedA0.first = 2;
    edgeCreatedA0.second = 0;

    edgeCreatedA1.first = 2;
    edgeCreatedA1.second = 1;
  } 
  else if(D0D2>0.0f)                                   
    {                                                   
    /* here we know that d0d1<=0.0 */             
    isect2(VERT1,VERT0,VERT2,VV1,VV0,VV2,D1,D0,D2,isect0,isect1,isectpoint0,isectpoint1,tmpVars);
    edgeCreatedA0.first = 1;
    edgeCreatedA0.second = 0;

    edgeCreatedA1.first = 1;
    edgeCreatedA1.second = 2;
  }                                                  
  else if(sgn(D1)*sgn(D2)>0.0f || D0!=0.0f)  //TODO: use signal here instead of multiplication... 
  {                                   
    /* here we know that d0d1<=0.0 or that D0!=0.0 */
    isect2(VERT0,VERT1,VERT2,VV0,VV1,VV2,D0,D1,D2,isect0,isect1,isectpoint0,isectpoint1,tmpVars);   

    edgeCreatedA0.first = 0;
    edgeCreatedA0.second = 1;

    edgeCreatedA1.first = 0;
    edgeCreatedA1.second = 2;
  }                                                  
  else if(D1!=0.0f)                                  
  {       
    assert(false); //one of the edges touch the plane?                                       
    isect2(VERT1,VERT0,VERT2,VV1,VV0,VV2,D1,D0,D2,isect0,isect1,isectpoint0,isectpoint1,tmpVars); 
  }                                         
  else if(D2!=0.0f)                                  
  {        
    assert(false);                                           
    isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,isect0,isect1,isectpoint0,isectpoint1,tmpVars);     
  }                                                 
  else                                               
  {                                                   
    /* triangles are coplanar */    
    return 1;
  }
  return 0;
}





int MeshIntersectionGeometry::intersectTwoTriangles(const InputTriangle &triMesh0,const InputTriangle &triMesh1,
				     Point &coordsPt1,VertexFromIntersection &vertexThatCreatedPt1, Point &coordsPt2,
             VertexFromIntersection &vertexThatCreatedPt2, TempVarsComputeIntersections &tempVars)
{
  Point &V0 = getCoordinates(*triMesh0.getInputVertex(0));
  Point &V1 = getCoordinates(*triMesh0.getInputVertex(1));
  Point &V2 = getCoordinates(*triMesh0.getInputVertex(2));

  Point &U0 = getCoordinates(*triMesh1.getInputVertex(0));
  Point &U1 = getCoordinates(*triMesh1.getInputVertex(1));
  Point &U2 = getCoordinates(*triMesh1.getInputVertex(2));

  
  //VertCoord du0du1,du0du2,dv0dv1,dv0dv2;
  // VertCoord &du0du1 = *(tempRationals+39);
  // VertCoord &du0du2 = *(tempRationals+40);
  // VertCoord &dv0dv1 = *(tempRationals+41);
  // VertCoord &dv0dv2 = *(tempRationals+42);
  int du0du1;
  int du0du2;
  int dv0dv1;
  int dv0dv2;

  short index;
  

  //until here we've used up to tempRationals+55 (next available is tempRationals+56)
  int smallest1,smallest2;
  
  /* compute plane equation of triangle(V0,V1,V2) */
  /*SUB(E1,V1,V0);
  SUB(E2,V2,V0);
  CROSS(N1,E1,E2,tmp);
  //d1=-DOT(N1,V0);
  DOT(d1,N1,V0,tmp);
  d1 *= -1;
  */

  const PlaneEquation &equationTri0 = getPlaneEquationInputTriangle(0, &triMesh0-&(inputTriangles[0][0]),tempVars.tempVarsComputeEquation);
  const Point &N1 = equationTri0.normal;
  const VertCoord &d1 = equationTri0.d;



  /* plane equation 1: N1.X+d1=0 */

  /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
  //du0=DOT(N1,U0)+d1;
  DOT(tempVars.du0,N1,U0,tempVars.tmp);
  tempVars.du0 += d1;

  //du1=DOT(N1,U1)+d1;
  DOT(tempVars.du1,N1,U1,tempVars.tmp);
  tempVars.du1 += d1;

  //du2=DOT(N1,U2)+d1;
  DOT(tempVars.du2,N1,U2,tempVars.tmp);
  tempVars.du2 += d1;

  /* coplanarity robustness check */
/*#if USE_EPSILON_TEST==TRUE
  if(fabs(du0)<EPSILON) du0=0.0;
  if(fabs(du1)<EPSILON) du1=0.0;
  if(fabs(du2)<EPSILON) du2=0.0;
#endif*/
  //du0du1=du0*du1;

  //TODO: SoS here...
  //du0du1=du0; //TODO: remove this multiplication! we only need the sign!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //du0du1*=du1;
  du0du1 = sgn(tempVars.du0)*sgn(tempVars.du1);
  //du0du2=du0*du2;
  //du0du2=du0;
  //du0du2*=du2;
  du0du2 = sgn(tempVars.du0)*sgn(tempVars.du2);

  if(du0du1>0 && du0du2>0) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute plane of triangle (U0,U1,U2) */
 /* SUB(E1,U1,U0);
  SUB(E2,U2,U0);
  CROSS(N2,E1,E2,tmp);
  //d2=-DOT(N2,U0);
  DOT(d2,N2,U0,tmp);
  d2 *= -1;*/
  const PlaneEquation &equationTri1 = getPlaneEquationInputTriangle(1, &triMesh1-&(inputTriangles[1][0]),tempVars.tempVarsComputeEquation);
  const Point &N2 = equationTri1.normal;
  const VertCoord &d2 = equationTri1.d;

  /* plane equation 2: N2.X+d2=0 */

  /* put V0,V1,V2 into plane equation 2 */
  //dv0=DOT(N2,V0)+d2;
  DOT(tempVars.dv0,N2,V0,tempVars.tmp);
  tempVars.dv0 += d2;
  //dv1=DOT(N2,V1)+d2;
  DOT(tempVars.dv1,N2,V1,tempVars.tmp);
  tempVars.dv1 += d2;
  //dv2=DOT(N2,V2)+d2;
  DOT(tempVars.dv2,N2,V2,tempVars.tmp);
  tempVars.dv2 += d2;

/*
#if USE_EPSILON_TEST==TRUE
  if(fabs(dv0)<EPSILON) dv0=0.0;
  if(fabs(dv1)<EPSILON) dv1=0.0;
  if(fabs(dv2)<EPSILON) dv2=0.0;
#endif*/

  //dv0,dv1,dv2 are the distances from v0 to the plane 2..
  //du0,du1,du2 are the distances from u0, ... to plane 1 

  //dv0dv1=dv0*dv1;
  //dv0dv1=dv0;
  //dv0dv1*=dv1;
  dv0dv1 = sgn(tempVars.dv0)*sgn(tempVars.dv1);
  //dv0dv2=dv0*dv2;
  //dv0dv2=dv0;
  //dv0dv2*=dv2;
  dv0dv2 = sgn(tempVars.dv0)*sgn(tempVars.dv2);
        
  if(dv0dv1>0 && dv0dv2>0) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */



  /* compute direction of intersection line */
  // L = tD + O , for some point O on the line...
  // D = N1 x N2 (cross)
  CROSS(tempVars.D,N1,N2,tempVars.tmp);

  /* compute an index to the largest component of D */
  //max=fabs(D[0]);
  fabs(tempVars.max,tempVars.D[0]);
  index=0;
  //b=fabs(D[1]);
  //c=fabs(D[2]);
  fabs(tempVars.b,tempVars.D[1]);
  fabs(tempVars.c,tempVars.D[2]);
  if(tempVars.b>tempVars.max) tempVars.max=tempVars.b,index=1;
  if(tempVars.c>tempVars.max) tempVars.max=tempVars.c,index=2;

  /* this is the simplified projection onto L*/

  //TODO: no need to copy here: use reference...
  pair<int,int> edgeCreatedA1;
  pair<int,int> edgeCreatedA2;
  /* compute interval for triangle 1 */
  int coplanar=compute_intervals_isectline(V0,V1,V2,V0[index],V1[index],V2[index],tempVars.dv0,tempVars.dv1,tempVars.dv2,
				       dv0dv1,dv0dv2,tempVars.isect1[0],tempVars.isect1[1],tempVars.isectpointA1,tempVars.isectpointA2,edgeCreatedA1,edgeCreatedA2,tempVars.tempRationals);

  //SoS: will never happen...
  if(coplanar) return 0;     


  pair<int,int> edgeCreatedB1;
  pair<int,int> edgeCreatedB2;
  /* compute interval for triangle 2 */
  compute_intervals_isectline(U0,U1,U2,U0[index],U1[index],U2[index],tempVars.du0,tempVars.du1,tempVars.du2,
			      du0du1,du0du2,tempVars.isect2[0],tempVars.isect2[1],tempVars.isectpointB1,tempVars.isectpointB2,edgeCreatedB1,edgeCreatedB2,tempVars.tempRationals);

  SORT2(tempVars.isect1[0],tempVars.isect1[1],smallest1,tempVars.tmp); //smallest1 == 0 iff isect1[0] < isect1[1]
  SORT2(tempVars.isect2[0],tempVars.isect2[1],smallest2,tempVars.tmp);

  if(tempVars.isect1[1]<tempVars.isect2[0] || tempVars.isect2[1]<tempVars.isect1[0]) return 0;

  /* at this point, we know that the triangles intersect */

  if(tempVars.isect2[0]<tempVars.isect1[0])
  {
    if(smallest1==0) { 
      SET(coordsPt1,tempVars.isectpointA1); 
      //the first vertex is vertex A1
      //A1 is the intersection between edgeCreatedA1 and the plane of the second triangle
      vertexThatCreatedPt1 = VertexFromIntersection(*triMesh0.getInputVertex(edgeCreatedA1.first), *triMesh0.getInputVertex(edgeCreatedA1.second), 
                                                      triMesh1);

    } // if(isect1[0] < isect1[1]) copy the point isectpointA1 to the output isectpoint1..
    else { 
      SET(coordsPt1,tempVars.isectpointA2); 

      vertexThatCreatedPt1 = VertexFromIntersection(*triMesh0.getInputVertex(edgeCreatedA2.first), *triMesh0.getInputVertex(edgeCreatedA2.second), 
                                                      triMesh1);
    }

    if(tempVars.isect2[1]<tempVars.isect1[1])
    {
      if(smallest2==0) { 
        SET(coordsPt2,tempVars.isectpointB2); 

        vertexThatCreatedPt2 = VertexFromIntersection(*triMesh1.getInputVertex(edgeCreatedB2.first), *triMesh1.getInputVertex(edgeCreatedB2.second), 
                                                      triMesh0);
      }
      else { 
        SET(coordsPt2,tempVars.isectpointB1); 

        vertexThatCreatedPt2 = VertexFromIntersection(*triMesh1.getInputVertex(edgeCreatedB1.first), *triMesh1.getInputVertex(edgeCreatedB1.second), 
                                                      triMesh0);
      }
    }
    else
    {
      if(smallest1==0) { 
        SET(coordsPt2,tempVars.isectpointA2); 

        vertexThatCreatedPt2 = VertexFromIntersection(*triMesh0.getInputVertex(edgeCreatedA2.first), *triMesh0.getInputVertex(edgeCreatedA2.second), 
                                                      triMesh1);
      }
      else { 
        SET(coordsPt2,tempVars.isectpointA1);

        vertexThatCreatedPt2 = VertexFromIntersection(*triMesh0.getInputVertex(edgeCreatedA1.first), *triMesh0.getInputVertex(edgeCreatedA1.second), 
                                                      triMesh1);
      }
    }
  }
  else
  {
    if(smallest2==0) { 
      SET(coordsPt1,tempVars.isectpointB1); 

      vertexThatCreatedPt1 = VertexFromIntersection(*triMesh1.getInputVertex(edgeCreatedB1.first), *triMesh1.getInputVertex(edgeCreatedB1.second), 
                                                      triMesh0);
    }
    else { 
      SET(coordsPt1,tempVars.isectpointB2); 

      vertexThatCreatedPt1 = VertexFromIntersection(*triMesh1.getInputVertex(edgeCreatedB2.first), *triMesh1.getInputVertex(edgeCreatedB2.second), 
                                                      triMesh0);
    }

    if(tempVars.isect2[1]>tempVars.isect1[1])
    {
      if(smallest1==0) { 
        SET(coordsPt2,tempVars.isectpointA2); 

        vertexThatCreatedPt2 = VertexFromIntersection(*triMesh0.getInputVertex(edgeCreatedA2.first), *triMesh0.getInputVertex(edgeCreatedA2.second), 
                                                      triMesh1);
      }
      else { 
        SET(coordsPt2,tempVars.isectpointA1); 

        vertexThatCreatedPt2 = VertexFromIntersection(*triMesh0.getInputVertex(edgeCreatedA1.first), *triMesh0.getInputVertex(edgeCreatedA1.second), 
                                                      triMesh1);
      }      
    }
    else
    {
      if(smallest2==0) { 
        SET(coordsPt2,tempVars.isectpointB2); 

        vertexThatCreatedPt2 = VertexFromIntersection(*triMesh1.getInputVertex(edgeCreatedB2.first), *triMesh1.getInputVertex(edgeCreatedB2.second), 
                                                      triMesh0);
      }
      else { 
        SET(coordsPt2,tempVars.isectpointB1); 

        vertexThatCreatedPt2 = VertexFromIntersection(*triMesh1.getInputVertex(edgeCreatedB1.first), *triMesh1.getInputVertex(edgeCreatedB1.second), 
                                                      triMesh0);
      } 
    }
  }
  return 1;
}
