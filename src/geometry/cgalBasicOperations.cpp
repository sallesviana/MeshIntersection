#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Projection_traits_xz_3.h>
#include <CGAL/Projection_traits_yz_3.h>
typedef CGAL::Projection_traits_xy_3<Kernel>  Projection_xy;
typedef CGAL::Projection_traits_xz_3<Kernel>  Projection_xz;
typedef CGAL::Projection_traits_yz_3<Kernel>  Projection_yz;
typedef Projection_xy::Point_2 Point2_xy;
typedef Projection_xz::Point_2 Point2_xz;
typedef Projection_yz::Point_2 Point2_yz;


const CGAL::Interval_nt<false> getCoord(const Point_3 &p, int coord) {
  switch(coord) {
    case 0: return p.x().approx();
    case 1: return p.y().approx();
    case 2: return p.z().approx();
    default: assert(false);
  }
}
//projection onto z=0
int orientationCGAL(const Point_3 &p0,const Point_3 &p1,const Point_3 &p2,int whatPlaneProjectOnto) {
  int coordY = 1;
  int coordX = 0;
  if(whatPlaneProjectOnto==PLANE_X0) { //if the triangle is projected to X=0 --> we need to use coordinates y,z (instead of x,y)
    coordX = 1;
    coordY = 2;
  } else if(whatPlaneProjectOnto ==PLANE_Y0) { //if the triangle is projected to Y=0 --> we need to use coordinates z,x (instead of x,y)
    coordX = 2;
    coordY = 0;
  }
  //CGAL::Orientation cgalAns = CGAL::orientation(p0,p1,p2,Point2_xy);
  const CGAL::Interval_nt<false> p0x = getCoord(p0,coordX);
  const CGAL::Interval_nt<false> p0y = getCoord(p0,coordY);

  const CGAL::Interval_nt<false> p1x = getCoord(p1,coordX);
  const CGAL::Interval_nt<false> p1y = getCoord(p1,coordY);

  const CGAL::Interval_nt<false> px = getCoord(p2,coordX);
  const CGAL::Interval_nt<false> py = getCoord(p2,coordY);

  const CGAL::Interval_nt<false>  ans = (p1x-p0x)*(py-p0y) -  (p1y-p0y)*(px-p0x);

  if(ans<0) return -1;
  if(ans>0) return 1;
  return 0;

}

//projection onto z=0
int orientationCGAL(const Point_3 &p0,const Point_3 &p1,const Point_3 &p2) {
  //CGAL::Orientation cgalAns = CGAL::orientation(p0,p1,p2,Point2_xy);
  const CGAL::Interval_nt<false> p0x = p0.x().approx();
  const CGAL::Interval_nt<false> p0y = p0.y().approx();

  const CGAL::Interval_nt<false> p1x = p1.x().approx();
  const CGAL::Interval_nt<false> p1y = p1.y().approx();

  const CGAL::Interval_nt<false> px = p2.x().approx();
  const CGAL::Interval_nt<false> py = p2.y().approx();

  const CGAL::Interval_nt<false>  ans = (p1x-p0x)*(py-p0y) -  (p1y-p0y)*(px-p0x);
  
  if(ans<0) return -1;
  if(ans>0) return 1;
  return 0;

 // int ans = sgn( (p1[coordX]-p0[coordX])*(p[coordY]-p0[coordY]) -  (p1[coordY]-p0[coordY])*(p[coordX]-p0[coordX]) );

  //int ans = 0;//signDeterminant4(getCoordinates(p1),getCoordinates(p2),getCoordinates(p3),getCoordinates(v));
  //if(cgalAns==CGAL::NEGATIVE) ans = 1;
  //else if (cgalAns==CGAL::POSITIVE) ans = -1;
}


int isVertexConvexMainImplCGAL(const Point_3 &p0,const Point_3 &p1,const Point_3 &p2, int whatPlaneProjectOnto)  {
  int coordY = 1;
  int coordX = 0;
  if(whatPlaneProjectOnto==PLANE_X0) { //if the triangle is projected to X=0 --> we need to use coordinates y,z (instead of x,y)
    coordX = 1;
    coordY = 2;
  } else if(whatPlaneProjectOnto ==PLANE_Y0) { //if the triangle is projected to Y=0 --> we need to use coordinates z,x (instead of x,y)
    coordX = 2;
    coordY = 0;
  }
  //CGAL::Orientation cgalAns = CGAL::orientation(p0,p1,p2,Point2_xy);
  const CGAL::Interval_nt<false> p1x = getCoord(p0,coordX);
  const CGAL::Interval_nt<false> p1y = getCoord(p0,coordY);

  const CGAL::Interval_nt<false> p2x = getCoord(p1,coordX);
  const CGAL::Interval_nt<false> p2y = getCoord(p1,coordY);

  const CGAL::Interval_nt<false> p3x = getCoord(p2,coordX);
  const CGAL::Interval_nt<false> p3y = getCoord(p2,coordY);  

  /* //REMOVE....
  VertCoord twoArea = vertex1[coordX]*(vertex2[coordY]-vertex3[coordY]) +
                      vertex2[coordX]*(vertex3[coordY]-vertex1[coordY]) +
                      vertex3[coordX]*(vertex1[coordY]-vertex2[coordY]) ;

  listVerticesToProcess[vertexId].convex = twoArea<0;
  return listVerticesToProcess[vertexId].convex;
  */
 // VertCoord *tempCoords  = tempVars.tempCoords;
  
  //tempCoords[0] = vertex2[coordY]; //VertCoord twoArea = vertex1[coordX]*(vertex2[coordY]-vertex3[coordY]);
                                   // + vertex2[coordX]*(vertex3[coordY]-vertex1[coordY]) + vertex3[coordX]*(vertex1[coordY]-vertex2[coordY]);
  
  const CGAL::Interval_nt<false>  ans = ( p1x*(p2y-p3y) + p2x*(p3y-p1y) + p3x*(p1y-p2y) ) ;
  if(ans<0) return 1;
  if(ans>0) return -1;
  return 0;

  //return -sgn(tempCoords[0]); //if the sign is negative --> the vertex is convex...
}