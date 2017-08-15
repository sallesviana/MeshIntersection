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
inline void isect2(const Point &VTX0,const Point & VTX1,const Point & VTX2,const VertCoord &VV0, const VertCoord &VV1,const VertCoord &VV2,
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
inline void compute_intervals_isectline(const Point &VERT0,const Point & VERT1,const Point &VERT2,
				       const VertCoord &VV0,const VertCoord &VV1,const VertCoord &VV2,const VertCoord &D0,const VertCoord &D1,const VertCoord &D2,
				       const int D0D1,const int D0D2,VertCoord &isect0,VertCoord & isect1,
				       Point &isectpoint0,Point & isectpoint1, pair<int,int> &edgeCreatedA0, pair<int,int> &edgeCreatedA1, 
               const Point &v0, const Point &v1, const Point &v2, const Point &u0, const Point &u1, const Point &u2, //for debugging purposes..
               VertCoord * tmpVars)
{

  //TODO: SoS when triangles touch...

  if(D0D1>0)                                        
  {                                                    
    /* here we know that D0D2<=0.0 */                  
    /* that is D0, D1 are on the same side, D2 on the other or on the plane */
    isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,isect0,isect1,isectpoint0,isectpoint1,tmpVars);

    edgeCreatedA0.first = 2;
    edgeCreatedA0.second = 0;

    edgeCreatedA1.first = 2;
    edgeCreatedA1.second = 1;
  } 
  else if(D0D2>0)                                   
    {                                                   
    /* here we know that d0d1<=0.0 */             
    isect2(VERT1,VERT0,VERT2,VV1,VV0,VV2,D1,D0,D2,isect0,isect1,isectpoint0,isectpoint1,tmpVars);
    edgeCreatedA0.first = 1;
    edgeCreatedA0.second = 0;

    edgeCreatedA1.first = 1;
    edgeCreatedA1.second = 2;
  }                                                  
  else if(sgn(D1)*sgn(D2)>0 || sgn(D0)!=0)  //TODO: use signal here instead of multiplication... 
  {                                   
    /* here we know that d0d1<=0.0 or that D0!=0.0 */
    isect2(VERT0,VERT1,VERT2,VV0,VV1,VV2,D0,D1,D2,isect0,isect1,isectpoint0,isectpoint1,tmpVars);   

    edgeCreatedA0.first = 0;
    edgeCreatedA0.second = 1;

    edgeCreatedA1.first = 0;
    edgeCreatedA1.second = 2;
  }                                                  
  else if(sgn(D1)!=0)                                  
  {       
    //one of the edges touch the plane? (this is a coincidence...)                                 
    isect2(VERT1,VERT0,VERT2,VV1,VV0,VV2,D1,D0,D2,isect0,isect1,isectpoint0,isectpoint1,tmpVars); 

    edgeCreatedA0.first = 1;
    edgeCreatedA0.second = 0;

    edgeCreatedA1.first = 1;
    edgeCreatedA1.second = 2;
  }                                         
  else if(sgn(D2)!=0)                                  
  {                                                 
    isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,isect0,isect1,isectpoint0,isectpoint1,tmpVars);    

    edgeCreatedA0.first = 2;
    edgeCreatedA0.second = 0;

    edgeCreatedA1.first = 2;
    edgeCreatedA1.second = 1; 
  }                                                 
  else                                               
  {     
    #pragma omp critical
    {
      cerr << "Assertion failure in tri-tri intersection" << endl;
      cerr << "D0,D1,D2 : " << endl;
      cerr << D0.get_d() << endl;
      cerr << D1.get_d() << endl;
      cerr << D2.get_d() << endl;
      cerr << "Triangles: " << endl;
      for(int i=0;i<3;i++) cerr << v0[i].get_d() << endl; cerr << endl;     
      for(int i=0;i<3;i++) cerr << v1[i].get_d() << endl; cerr << endl;           
      for(int i=0;i<3;i++) cerr << v2[i].get_d() << endl; cerr << endl;           

      for(int i=0;i<3;i++) cerr << u0[i].get_d() << endl; cerr << endl;           
      for(int i=0;i<3;i++) cerr << u1[i].get_d() << endl; cerr << endl;           
      for(int i=0;i<3;i++) cerr << u2[i].get_d() << endl; cerr << endl;                                                   
      cerr << endl;
    }
    /* triangles are coplanar */    
    assert(false); //we should've detected this before...
  }
  
}





