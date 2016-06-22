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

/* if USE_EPSILON_TEST is true then we do a check: 
         if |dv|<EPSILON then dv=0.0;
   else no check is done (which is less robust)
*/
#define USE_EPSILON_TEST FALSE  
#define EPSILON 0.000001


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

void DOT(VertCoord &dest, VertCoord v1[3], VertCoord v2[3], VertCoord &tmp) {
  dest = v1[0];
  dest *= v2[0];
  tmp = v1[1];
  tmp *= v2[1];
  dest += tmp;
  tmp = v1[2];
  tmp *= v2[2];
  dest += tmp;

  
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
int coplanar_tri_tri(VertCoord N[3],VertCoord V0[3],VertCoord V1[3],VertCoord V2[3],
                     VertCoord U0[3],VertCoord U1[3],VertCoord U2[3],VertCoord *Temp)
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
inline void isect2(VertCoord VTX0[3],VertCoord VTX1[3],VertCoord VTX2[3],const VertCoord &VV0, const VertCoord &VV1,const VertCoord &VV2,
	    const VertCoord &D0, const VertCoord &D1, const VertCoord &D2,VertCoord *isect0,VertCoord *isect1,VertCoord isectpoint0[3],VertCoord isectpoint1[3], VertCoord *tmpVars) 
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
  *isect0   =VV1; //+(VV1-VV0)*tmp; 
  *isect0  -=VV0;
  *isect0 *= tmp;
  *isect0 += VV0;

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
  *isect1 = VV2;
  *isect1 -= VV0;
  *isect1 *= tmp;
  *isect1 += VV0;
  SUB(diff,VTX2,VTX0);                   
  //MULT(diff,diff,tmp);      
  diff[0]*=tmp;
  diff[1]*=tmp;
  diff[2]*=tmp;         
  ADD(isectpoint1,VTX0,diff);          
}


#if 0
#define ISECT2(VTX0,VTX1,VTX2,VV0,VV1,VV2,D0,D1,D2,isect0,isect1,isectpoint0,isectpoint1) \
              tmp=D0/(D0-D1);                    \
              isect0=VV0+(VV1-VV0)*tmp;          \
	      SUB(diff,VTX1,VTX0);               \
	      MULT(diff,diff,tmp);               \
              ADD(isectpoint0,diff,VTX0);        \
              tmp=D0/(D0-D2);                    
/*              isect1=VV0+(VV2-VV0)*tmp;          \ */
/*              SUB(diff,VTX2,VTX0);               \     */
/*              MULT(diff,diff,tmp);               \   */
/*              ADD(isectpoint1,VTX0,diff);           */
#endif

inline int compute_intervals_isectline(VertCoord VERT0[3],VertCoord VERT1[3],VertCoord VERT2[3],
				       const VertCoord &VV0,const VertCoord &VV1,const VertCoord &VV2,const VertCoord &D0,const VertCoord &D1,const VertCoord &D2,
				       const VertCoord &D0D1,const VertCoord &D0D2,VertCoord *isect0,VertCoord *isect1,
				       VertCoord isectpoint0[3],VertCoord isectpoint1[3], VertCoord * tmpVars)
{
  if(D0D1>0.0f)                                        
  {                                                    
    /* here we know that D0D2<=0.0 */                  
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */
    isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,isect0,isect1,isectpoint0,isectpoint1,tmpVars);
  } 
  else if(D0D2>0.0f)                                   
    {                                                   
    /* here we know that d0d1<=0.0 */             
    isect2(VERT1,VERT0,VERT2,VV1,VV0,VV2,D1,D0,D2,isect0,isect1,isectpoint0,isectpoint1,tmpVars);
  }                                                  
  else if(D1*D2>0.0f || D0!=0.0f)   
  {                                   
    /* here we know that d0d1<=0.0 or that D0!=0.0 */
    isect2(VERT0,VERT1,VERT2,VV0,VV1,VV2,D0,D1,D2,isect0,isect1,isectpoint0,isectpoint1,tmpVars);   
  }                                                  
  else if(D1!=0.0f)                                  
  {                                               
    isect2(VERT1,VERT0,VERT2,VV1,VV0,VV2,D1,D0,D2,isect0,isect1,isectpoint0,isectpoint1,tmpVars); 
  }                                         
  else if(D2!=0.0f)                                  
  {                                                   
    isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,isect0,isect1,isectpoint0,isectpoint1,tmpVars);     
  }                                                 
  else                                               
  {                                                   
    /* triangles are coplanar */    
    return 1;
  }
  return 0;
}

#define COMPUTE_INTERVALS_ISECTLINE(VERT0,VERT1,VERT2,VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,isect0,isect1,isectpoint0,isectpoint1,tmpVars) \
  if(D0D1>0.0f)                                         \
  {                                                     \
    /* here we know that D0D2<=0.0 */                   \
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
    isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,&isect0,&isect1,isectpoint0,isectpoint1,tmpVars);          \
  }                                                     
#if 0
  else if(D0D2>0.0f)                                    \
  {                                                     \
    /* here we know that d0d1<=0.0 */                   \
    isect2(VERT1,VERT0,VERT2,VV1,VV0,VV2,D1,D0,D2,&isect0,&isect1,isectpoint0,isectpoint1);          \
  }                                                     \
  else if(D1*D2>0.0f || D0!=0.0f)                       \
  {                                                     \
    /* here we know that d0d1<=0.0 or that D0!=0.0 */   \
    isect2(VERT0,VERT1,VERT2,VV0,VV1,VV2,D0,D1,D2,&isect0,&isect1,isectpoint0,isectpoint1);          \
  }                                                     \
  else if(D1!=0.0f)                                     \
  {                                                     \
    isect2(VERT1,VERT0,VERT2,VV1,VV0,VV2,D1,D0,D2,&isect0,&isect1,isectpoint0,isectpoint1);          \
  }                                                     \
  else if(D2!=0.0f)                                     \
  {                                                     \
    isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,&isect0,&isect1,isectpoint0,isectpoint1);          \
  }                                                     \
  else                                                  \
  {                                                     \
    /* triangles are coplanar */                        \
    coplanar=1;                                         \
    return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);      \
  }
#endif

int tri_tri_intersect_with_isectline(VertCoord V0[3],VertCoord V1[3],VertCoord V2[3],
				     VertCoord U0[3],VertCoord U1[3],VertCoord U2[3],int *coplanar,
				     VertCoord isectpt1[3],VertCoord isectpt2[3], VertCoord *tempRationals)
{
  //VertCoord E1[3],E2[3];
  VertCoord *E1 = tempRationals;
  VertCoord *E2 = tempRationals+3;

  //VertCoord N1[3],N2[3],d1,d2;
  VertCoord *N1 = tempRationals+6;
  VertCoord *N2 = tempRationals+9;
  VertCoord &d1 = *(tempRationals+12);
  VertCoord &d2 = *(tempRationals+13);
  
  //VertCoord du0,du1,du2,dv0,dv1,dv2;
  VertCoord &du0 = *(tempRationals+14);
  VertCoord &du1 = *(tempRationals+15);
  VertCoord &du2 = *(tempRationals+16);
  VertCoord &dv0 = *(tempRationals+17);
  VertCoord &dv1 = *(tempRationals+18);
  VertCoord &dv2 = *(tempRationals+19);

  //VertCoord D[3];
  VertCoord *D = tempRationals+20;
  // VertCoord isect1[2], isect2[2];
  VertCoord *isect1 = tempRationals+23;
  VertCoord *isect2 = tempRationals+25;
  //VertCoord isectpointA1[3],isectpointA2[3];
  VertCoord *isectpointA1 = tempRationals+27;
  VertCoord *isectpointA2 = tempRationals+30;
  //VertCoord isectpointB1[3],isectpointB2[3];
  VertCoord *isectpointB1 = tempRationals+33;
  VertCoord *isectpointB2 = tempRationals+36;

  //VertCoord du0du1,du0du2,dv0dv1,dv0dv2;
  VertCoord &du0du1 = *(tempRationals+39);
  VertCoord &du0du2 = *(tempRationals+40);
  VertCoord &dv0dv1 = *(tempRationals+41);
  VertCoord &dv0dv2 = *(tempRationals+42);

  short index;
  //VertCoord vp0,vp1,vp2;
  VertCoord &vp0 = *(tempRationals+43);
  VertCoord &vp1 = *(tempRationals+44);
  VertCoord &vp2 = *(tempRationals+45);
  //VertCoord up0,up1,up2;
  VertCoord &up0 = *(tempRationals+46);
  VertCoord &up1 = *(tempRationals+47);
  VertCoord &up2 = *(tempRationals+48);
  //VertCoord b,c,max;
  VertCoord &b = *(tempRationals+49);
  VertCoord &c = *(tempRationals+50);
  VertCoord &max = *(tempRationals+51);
  //VertCoord tmp,diff[3];
  VertCoord &tmp = *(tempRationals+52);
  VertCoord *diff = (tempRationals+53);

  //until here we've used up to tempRationals+55 (next available is tempRationals+56)
  int smallest1,smallest2;
  
  /* compute plane equation of triangle(V0,V1,V2) */
  SUB(E1,V1,V0);
  SUB(E2,V2,V0);
  CROSS(N1,E1,E2,tmp);
  //d1=-DOT(N1,V0);
  DOT(d1,N1,V0,tmp);
  d1 *= -1;

  /* plane equation 1: N1.X+d1=0 */

  /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
  //du0=DOT(N1,U0)+d1;
  DOT(du0,N1,U0,tmp);
  du0 += d1;

  //du1=DOT(N1,U1)+d1;
  DOT(du1,N1,U1,tmp);
  du1 += d1;

  //du2=DOT(N1,U2)+d1;
  DOT(du2,N1,U2,tmp);
  du2 += d1;

  /* coplanarity robustness check */
/*#if USE_EPSILON_TEST==TRUE
  if(fabs(du0)<EPSILON) du0=0.0;
  if(fabs(du1)<EPSILON) du1=0.0;
  if(fabs(du2)<EPSILON) du2=0.0;
#endif*/
  //du0du1=du0*du1;
  du0du1=du0;
  du0du1*=du1;
  //du0du2=du0*du2;
  du0du2=du0;
  du0du2*=du2;

  if(du0du1>0 && du0du2>0) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute plane of triangle (U0,U1,U2) */
  SUB(E1,U1,U0);
  SUB(E2,U2,U0);
  CROSS(N2,E1,E2,tmp);
  //d2=-DOT(N2,U0);
  DOT(d2,N2,U0,tmp);
  d2 *= -1;
  /* plane equation 2: N2.X+d2=0 */

  /* put V0,V1,V2 into plane equation 2 */
  //dv0=DOT(N2,V0)+d2;
  DOT(dv0,N2,V0,tmp);
  dv0 += d2;
  //dv1=DOT(N2,V1)+d2;
  DOT(dv1,N2,V1,tmp);
  dv1 += d2;
  //dv2=DOT(N2,V2)+d2;
  DOT(dv2,N2,V2,tmp);
  dv2 += d2;

/*
#if USE_EPSILON_TEST==TRUE
  if(fabs(dv0)<EPSILON) dv0=0.0;
  if(fabs(dv1)<EPSILON) dv1=0.0;
  if(fabs(dv2)<EPSILON) dv2=0.0;
#endif*/

  //dv0dv1=dv0*dv1;
  dv0dv1=dv0;
  dv0dv1*=dv1;
  //dv0dv2=dv0*dv2;
  dv0dv2=dv0;
  dv0dv2*=dv2;
        
  if(dv0dv1>0 && dv0dv2>0) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute direction of intersection line */
  CROSS(D,N1,N2,tmp);

  /* compute and index to the largest component of D */
  //max=fabs(D[0]);
  fabs(max,D[0]);
  index=0;
  //b=fabs(D[1]);
  //c=fabs(D[2]);
  fabs(b,D[1]);
  fabs(c,D[2]);
  if(b>max) max=b,index=1;
  if(c>max) max=c,index=2;

  /* this is the simplified projection onto L*/
  vp0=V0[index];
  vp1=V1[index];
  vp2=V2[index];
  
  up0=U0[index];
  up1=U1[index];
  up2=U2[index];

  /* compute interval for triangle 1 */
  *coplanar=compute_intervals_isectline(V0,V1,V2,vp0,vp1,vp2,dv0,dv1,dv2,
				       dv0dv1,dv0dv2,&isect1[0],&isect1[1],isectpointA1,isectpointA2,tempRationals+54);
  if(*coplanar) return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2,tempRationals+60);     


  /* compute interval for triangle 2 */
  compute_intervals_isectline(U0,U1,U2,up0,up1,up2,du0,du1,du2,
			      du0du1,du0du2,&isect2[0],&isect2[1],isectpointB1,isectpointB2,tempRationals+54);

  SORT2(isect1[0],isect1[1],smallest1,tmp);
  SORT2(isect2[0],isect2[1],smallest2,tmp);

  if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return 0;

  /* at this point, we know that the triangles intersect */

  if(isect2[0]<isect1[0])
  {
    if(smallest1==0) { SET(isectpt1,isectpointA1); }
    else { SET(isectpt1,isectpointA2); }

    if(isect2[1]<isect1[1])
    {
      if(smallest2==0) { SET(isectpt2,isectpointB2); }
      else { SET(isectpt2,isectpointB1); }
    }
    else
    {
      if(smallest1==0) { SET(isectpt2,isectpointA2); }
      else { SET(isectpt2,isectpointA1); }
    }
  }
  else
  {
    if(smallest2==0) { SET(isectpt1,isectpointB1); }
    else { SET(isectpt1,isectpointB2); }

    if(isect2[1]>isect1[1])
    {
      if(smallest1==0) { SET(isectpt2,isectpointA2); }
      else { SET(isectpt2,isectpointA1); }      
    }
    else
    {
      if(smallest2==0) { SET(isectpt2,isectpointB2); }
      else { SET(isectpt2,isectpointB1); } 
    }
  }
  return 1;
}
