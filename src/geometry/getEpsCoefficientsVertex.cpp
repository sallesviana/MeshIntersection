#include "sosPredicatesImpl.h"

/*
VertCoord SosPredicatesImpl::getEpsCoefficientsVertexFromIntersectionFaster2(const VertexFromIntersection &v0, int epsCoefficient,int coord) const {
	assert(coord<=2);
	assert(epsCoefficient<=3);

	int meshId = v0.getMeshOfTriangleDefiningVertex();

	

	const VertCoord &t0x = getCoordinates( *(v0.triangle.getInputVertex(0)) )[0];
	const VertCoord &t1x = getCoordinates( *(v0.triangle.getInputVertex(1)) )[0];
	const VertCoord &t2x = getCoordinates( *(v0.triangle.getInputVertex(2)) )[0];
	const VertCoord &r0x = getCoordinates(v0.edge[0])[0];
	const VertCoord &r1x = getCoordinates(v0.edge[1])[0];
	const VertCoord &t0y = getCoordinates( *(v0.triangle.getInputVertex(0)) )[1];
	const VertCoord &t1y = getCoordinates( *(v0.triangle.getInputVertex(1)) )[1];
	const VertCoord &t2y = getCoordinates( *(v0.triangle.getInputVertex(2)) )[1];
	const VertCoord &r0y = getCoordinates(v0.edge[0])[1];
	const VertCoord &r1y = getCoordinates(v0.edge[1])[1];
	const VertCoord &t0z = getCoordinates( *(v0.triangle.getInputVertex(0)) )[2];
	const VertCoord &t1z = getCoordinates( *(v0.triangle.getInputVertex(1)) )[2];
	const VertCoord &t2z = getCoordinates( *(v0.triangle.getInputVertex(2)) )[2];
	const VertCoord &r0z = getCoordinates(v0.edge[0])[2];
	const VertCoord &r1z = getCoordinates(v0.edge[1])[2];

	if(epsCoefficient==0) return getCoordinates(v0)[coord];

	const InputTriangle &triangle = v0.triangle;
	int triId = triangle.getIdForEps();
	MeshIntersectionGeometry::TEpsDeterminanCommonTerms &tCommonTerms = geometry->tEpsDeterminanCommonTerms[triId];


	if(tCommonTerms.init==false) {
		omp_set_lock(&tCommonTerms.lock);

		if(tCommonTerms.init==false) { //is it still false?
			//let's init
			//TODO: add syncrhonization...
			tCommonTerms.commonTerms[0] = (t0y*(t1x - t2x) + t1y*t2x - t1x*t2y + t0x*(t2y-t1y));
			tCommonTerms.commonTerms[1] = (t0z*(t1x - t2x) + t1z*t2x - t1x*t2z + t0x*(t2z-t1z));
			tCommonTerms.commonTerms[2] = (t0z*(t1y - t2y) + t1z*t2y - t1y*t2z + t0y*(t2z-t1z));
			tCommonTerms.init = true;
		}

		omp_unset_lock(&tCommonTerms.lock);
	}

	int vertexInterId = v0.getIdForEps();
  MeshIntersectionGeometry::EpsCoefficientsVertexIntersection &epsCoefficients = geometry->epsCoefficientsVertexIntersection[vertexInterId];



	if(epsCoefficients.init==false) {
		omp_set_lock(&epsCoefficients.lock);

		if(epsCoefficients.init==false)
		{
			rational tDenominatorValue = (r0z - r1z)*tCommonTerms.commonTerms[0] - (r0y - r1y)*tCommonTerms.commonTerms[1] + (r0x - r1x)*tCommonTerms.commonTerms[2];
		
			rational tEpsCoefficients[3];
			if(meshId==0) {
				
			  tEpsCoefficients[0] = (t0z* t1y - t0y *t1z - t0z *t2y + t1z* t2y + t0y* t2z - t1y* t2z)/tDenominatorValue;//(t0z*(t1y - t2y) + t1z* t2y - t1y* t2z + t0y* (-t1z + t2z))/tDenominatorValue;
			  tEpsCoefficients[1] = (-t0z *t1x + t0x *t1z + t0z *t2x - t1z *t2x - t0x *t2z + t1x *t2z)/tDenominatorValue;//(-t1z*t2x + t0z *(-t1x + t2x) + t0x *(t1z - t2z) + t1x *t2z)/tDenominatorValue;
			  tEpsCoefficients[2] = (t0y* t1x - t0x* t1y - t0y* t2x + t1y* t2x + t0x *t2y - t1x *t2y)/tDenominatorValue;//(t0y*(t1x - t2x) + t1y *t2x - t1x *t2y + t0x *(-t1y + t2y))/tDenominatorValue;
			} else {
				

				tEpsCoefficients[0] =  (-t1z*t2y + t0z *(-t1y + t2y) + t0y *(t1z - t2z) + t1y* t2z)/tDenominatorValue;
				tEpsCoefficients[1] =  (t0z *(t1x - t2x) + t1z *t2x - t1x *t2z + t0x *(-t1z + t2z))/tDenominatorValue;
				tEpsCoefficients[2] =  (-t1y *t2x + t0y *(-t1x + t2x) + t0x *(t1y - t2y) + t1x *t2y)/tDenominatorValue;
			}

			rational r1r0x = r1x-r0x;
			rational r1r0y = r1y-r0y;
			rational r1r0z = r1z-r0z;		


			for(int eps=0;eps<3;eps++) {
				if(meshId==0) { //the triangle is in mesh 0 --> the ray is in mesh 1 --> each coord of ray translated by eps, eps2,eps3
					epsCoefficients.epsCoefficients[eps][0] =  (eps==0?1:0) + tEpsCoefficients[eps]*r1r0x;
					epsCoefficients.epsCoefficients[eps][1] =  (eps==1?1:0) + tEpsCoefficients[eps]*r1r0y;
					epsCoefficients.epsCoefficients[eps][2] =  (eps==2?1:0) + tEpsCoefficients[eps]*r1r0z;
				} else { //the ray is in mesh 0 --> it does not have eps coefficients...
					epsCoefficients.epsCoefficients[eps][0] =  tEpsCoefficients[eps]*r1r0x;
					epsCoefficients.epsCoefficients[eps][1] =  tEpsCoefficients[eps]*r1r0y;
					epsCoefficients.epsCoefficients[eps][2] =  tEpsCoefficients[eps]*r1r0z;
				}
				
			}

			epsCoefficients.init=true;
		}
		omp_unset_lock(&epsCoefficients.lock);

	}

	return epsCoefficients.epsCoefficients[epsCoefficient-1][coord];
}
*/

const VertCoord&  SosPredicatesImpl::getEpsCoefficientsVertexFromIntersectionFaster2(const VertexFromIntersection &v0, int epsCoefficient,int coord) const {
	assert(coord<=2);
	assert(epsCoefficient<=3);

	int meshId = v0.getMeshOfTriangleDefiningVertex();

	

	const VertCoord &t0x = getCoordinates( *(v0.triangle.getInputVertex(0)) )[0];
	const VertCoord &t1x = getCoordinates( *(v0.triangle.getInputVertex(1)) )[0];
	const VertCoord &t2x = getCoordinates( *(v0.triangle.getInputVertex(2)) )[0];
	const VertCoord &r0x = getCoordinates(v0.edge[0])[0];
	const VertCoord &r1x = getCoordinates(v0.edge[1])[0];
	const VertCoord &t0y = getCoordinates( *(v0.triangle.getInputVertex(0)) )[1];
	const VertCoord &t1y = getCoordinates( *(v0.triangle.getInputVertex(1)) )[1];
	const VertCoord &t2y = getCoordinates( *(v0.triangle.getInputVertex(2)) )[1];
	const VertCoord &r0y = getCoordinates(v0.edge[0])[1];
	const VertCoord &r1y = getCoordinates(v0.edge[1])[1];
	const VertCoord &t0z = getCoordinates( *(v0.triangle.getInputVertex(0)) )[2];
	const VertCoord &t1z = getCoordinates( *(v0.triangle.getInputVertex(1)) )[2];
	const VertCoord &t2z = getCoordinates( *(v0.triangle.getInputVertex(2)) )[2];
	const VertCoord &r0z = getCoordinates(v0.edge[0])[2];
	const VertCoord &r1z = getCoordinates(v0.edge[1])[2];

	if(epsCoefficient==0) return getCoordinates(v0)[coord];

	const InputTriangle &triangle = v0.triangle;
	int triId = triangle.getIdForEps();
	if(triId >=  geometry->tEpsDeterminanCommonTerms.size()) {
		cerr << "Error with triId: " << triId << " "  <<  geometry->tEpsDeterminanCommonTerms.size() << endl;
	}
	assert(triId <  geometry->tEpsDeterminanCommonTerms.size());
	if(triId==-1) {
		cerr << "Error when processing a vertex from intersection... id = -1 "<< endl;
	}
	assert(triId!=-1);
	MeshIntersectionGeometry::TEpsDeterminanCommonTerms &tCommonTerms = geometry->tEpsDeterminanCommonTerms[triId];


	if(tCommonTerms.init==false) {
		omp_set_lock(&tCommonTerms.lock);

		if(tCommonTerms.init==false) { //is it still false?
			//let's init
			//TODO: add syncrhonization...
			tempVars.t1xt2x = t1x; //(t1x - t2x)
			tempVars.t1xt2x -= t2x;

			tCommonTerms.commonTerms[0] = t0y;//(t0y*(tempVars.t1xt2x) + t1y*t2x - t1x*t2y + t0x*(t2y-t1y));
			tCommonTerms.commonTerms[0] *= tempVars.t1xt2x;
			tempVars.tmp = t1y;
			tempVars.tmp *=t2x;
			tCommonTerms.commonTerms[0] += tempVars.tmp;
			tempVars.tmp=t1x;
			tempVars.tmp*=t2y;
			tCommonTerms.commonTerms[0] -= tempVars.tmp;
			tempVars.tmp = t2y;
			tempVars.tmp -= t1y;
			tempVars.tmp *= t0x;
			tCommonTerms.commonTerms[0] += tempVars.tmp;


			tCommonTerms.commonTerms[1] = t0z;//(t0z*(tempVars.t1xt2x) + t1z*t2x - t1x*t2z + t0x*(t2z-t1z));
			tCommonTerms.commonTerms[1] *= tempVars.t1xt2x;
			tempVars.tmp = t1z;
			tempVars.tmp *=t2x;
			tCommonTerms.commonTerms[1] += tempVars.tmp;
			tempVars.tmp=t1x;
			tempVars.tmp*=t2z;
			tCommonTerms.commonTerms[1] -= tempVars.tmp;
			tempVars.tmp = t2z;
			tempVars.tmp -= t1z;
			tempVars.tmp *= t0x;
			tCommonTerms.commonTerms[1] += tempVars.tmp;

			//tCommonTerms.commonTerms[2] = (t0z*(t1y - t2y) + t1z*t2y - t1y*t2z + t0y*(t2z-t1z));

			tCommonTerms.commonTerms[2] = t1y;
			tCommonTerms.commonTerms[2] -= t2y;
			tCommonTerms.commonTerms[2] *= t0z;
			tempVars.tmp = t1z;
			tempVars.tmp *=t2y;
			tCommonTerms.commonTerms[2] += tempVars.tmp;
			tempVars.tmp=t1y;
			tempVars.tmp*=t2z;
			tCommonTerms.commonTerms[2] -= tempVars.tmp;
			tempVars.tmp = t2z;
			tempVars.tmp -= t1z;
			tempVars.tmp *= t0y;
			tCommonTerms.commonTerms[2] += tempVars.tmp;


			tCommonTerms.init = true;
		}

		omp_unset_lock(&tCommonTerms.lock);
	}

	int vertexInterId = v0.getIdForEps();
  MeshIntersectionGeometry::EpsCoefficientsVertexIntersection &epsCoefficients = geometry->epsCoefficientsVertexIntersection[vertexInterId];



	if(epsCoefficients.init==false) {
		omp_set_lock(&epsCoefficients.lock);

		if(epsCoefficients.init==false)
		{
			tempVars.tDenominatorValue = r0z;
			tempVars.tDenominatorValue -= r1z;
			tempVars.tDenominatorValue *= tCommonTerms.commonTerms[0];

			tempVars.tmp = r0y;
			tempVars.tmp -= r1y;
			tempVars.tmp *= tCommonTerms.commonTerms[1];
			tempVars.tDenominatorValue-= tempVars.tmp;


			tempVars.tmp = r0x;
			tempVars.tmp -= r1x;
			tempVars.tmp *= tCommonTerms.commonTerms[2];
			tempVars.tDenominatorValue+= tempVars.tmp;

			//rational tDenominatorValue = (r0z - r1z)*tCommonTerms.commonTerms[0] - (r0y - r1y)*tCommonTerms.commonTerms[1] + (r0x - r1x)*tCommonTerms.commonTerms[2];
		
			rational tEpsCoefficients[3];
			if(meshId==0) {
				/*(t0z (t1y - t2y) + t1z t2y - t1y t2z + t0y (-t1z + t2z)) \[Epsilon] + 
				(-t1z t2x + t0z (-t1x + t2x) + t0x (t1z - t2z) + t1x t2z) \[Epsilon]^2 + 
				(t0y (t1x - t2x) + t1y t2x - t1x t2y + t0x (-t1y + t2y)) \[Epsilon]^3	
			  */
			  //tEpsCoefficients[0] = (t0z*(t1y - t2y) + t1z* t2y - t1y* t2z + t0y* (-t1z + t2z))/tempVars.tDenominatorValue;
				tEpsCoefficients[0] = t1y;
				tEpsCoefficients[0] -= t2y;
				tEpsCoefficients[0] *= t0z;

				tempVars.tmp = t1z;
				tempVars.tmp *= t2y;
				tEpsCoefficients[0]+= tempVars.tmp;

				tempVars.tmp = t1y;
				tempVars.tmp *= t2z;
				tEpsCoefficients[0]-= tempVars.tmp;

				tempVars.tmp = t2z;
				tempVars.tmp -= t1z;
				tempVars.tmp *= t0y;
				tEpsCoefficients[0]+= tempVars.tmp;

				tEpsCoefficients[0]/= tempVars.tDenominatorValue;


			  //tEpsCoefficients[1] = (-t1z*t2x + t0z *(-t1x + t2x) + t0x *(t1z - t2z) + t1x *t2z)/tempVars.tDenominatorValue;
				tEpsCoefficients[1] = t2x;
				tEpsCoefficients[1] -= t1x;
				tEpsCoefficients[1] *= t0z;

				tempVars.tmp = t2x;
				tempVars.tmp *= t1z;
				tEpsCoefficients[1]-= tempVars.tmp;

				tempVars.tmp = t1z;
				tempVars.tmp -= t2z;
				tempVars.tmp *= t0x;
				tEpsCoefficients[1]+= tempVars.tmp;

				tempVars.tmp = t2z;
				tempVars.tmp *= t1x;
				tEpsCoefficients[1]+= tempVars.tmp;

				tEpsCoefficients[1]/= tempVars.tDenominatorValue;



			  //tEpsCoefficients[2] = (t0y*(t1x - t2x) + t1y *t2x - t1x *t2y + t0x *(-t1y + t2y))/tempVars.tDenominatorValue;
				tEpsCoefficients[2] = t1x;
				tEpsCoefficients[2] -= t2x;
				tEpsCoefficients[2] *= t0y;

				tempVars.tmp = t1y;
				tempVars.tmp *= t2x;
				tEpsCoefficients[2]+= tempVars.tmp;

				tempVars.tmp = t1x;
				tempVars.tmp *= t2y;
				tEpsCoefficients[2]-= tempVars.tmp;

				tempVars.tmp = t2y;
				tempVars.tmp -= t1y;
				tempVars.tmp *= t0x;
				tEpsCoefficients[2]+= tempVars.tmp;

				tEpsCoefficients[2]/= tempVars.tDenominatorValue;
			} else {
				/* -t1z t2y + t0z (-t1y + t2y) + t0y (t1z - t2z) + t1y t2z
				// t0z (t1x - t2x) + t1z t2x - t1x t2z + t0x (-t1z + t2z)
				// -t1y t2x + t0y (-t1x + t2x) + t0x (t1y - t2y) + t1x t2y
				*/

				//tEpsCoefficients[0] =  (-t1z*t2y + t0z *(-t1y + t2y) + t0y *(t1z - t2z) + t1y* t2z)/tempVars.tDenominatorValue;
				tEpsCoefficients[0] = t2y;
				tEpsCoefficients[0] -= t1y;
				tEpsCoefficients[0] *= t0z;

				tempVars.tmp = t1z;
				tempVars.tmp *= t2y;
				tEpsCoefficients[0]-= tempVars.tmp;

				tempVars.tmp = t1z;
				tempVars.tmp -= t2z;
				tempVars.tmp *= t0y;
				tEpsCoefficients[0]+= tempVars.tmp;

				tempVars.tmp = t2z;
				tempVars.tmp *= t1y;
				tEpsCoefficients[0]+= tempVars.tmp;

				tEpsCoefficients[0]/= tempVars.tDenominatorValue;


				//tEpsCoefficients[1] =  (t0z *(t1x - t2x) + t1z *t2x - t1x *t2z + t0x *(-t1z + t2z))/tempVars.tDenominatorValue;
				tEpsCoefficients[1] = t1x;
				tEpsCoefficients[1] -= t2x;
				tEpsCoefficients[1] *= t0z;

				tempVars.tmp = t1z;
				tempVars.tmp *= t2x;
				tEpsCoefficients[1]+= tempVars.tmp;

				tempVars.tmp = t1x;
				tempVars.tmp *= t2z;
				tEpsCoefficients[1]-= tempVars.tmp;

				tempVars.tmp = t2z;
				tempVars.tmp -= t1z;
				tempVars.tmp *= t0x;
				tEpsCoefficients[1]+= tempVars.tmp;

				tEpsCoefficients[1]/= tempVars.tDenominatorValue;


				//tEpsCoefficients[2] =  (-t1y *t2x + t0y *(-t1x + t2x) + t0x *(t1y - t2y) + t1x *t2y)/tempVars.tDenominatorValue;
				tEpsCoefficients[2] = t2x;
				tEpsCoefficients[2] -= t1x;
				tEpsCoefficients[2] *= t0y;

				tempVars.tmp = t2x;
				tempVars.tmp *= t1y;
				tEpsCoefficients[2]-= tempVars.tmp;

				tempVars.tmp = t1y;
				tempVars.tmp -= t2y;
				tempVars.tmp *= t0x;
				tEpsCoefficients[2]+= tempVars.tmp;

				tempVars.tmp = t2y;
				tempVars.tmp *= t1x;
				tEpsCoefficients[2]+= tempVars.tmp;

				tEpsCoefficients[2]/= tempVars.tDenominatorValue;
			}

			tempVars.r1r0x = r1x;
			tempVars.r1r0x -= r0x;

			tempVars.r1r0y = r1y;
			tempVars.r1r0y -= r0y;

			tempVars.r1r0z = r1z;		
			tempVars.r1r0z -= r0z;		

			for(int eps=0;eps<3;eps++) {
				if(meshId==0) { //the triangle is in mesh 0 --> the ray is in mesh 1 --> each coord of ray translated by eps, eps2,eps3
					//epsCoefficients.epsCoefficients[eps][0] =  (eps==0?1:0) + tEpsCoefficients[eps]*tempVars.r1r0x;
					//epsCoefficients.epsCoefficients[eps][1] =  (eps==1?1:0) + tEpsCoefficients[eps]*tempVars.r1r0y;
					//epsCoefficients.epsCoefficients[eps][2] =  (eps==2?1:0) + tEpsCoefficients[eps]*tempVars.r1r0z;
					epsCoefficients.epsCoefficients[eps][0] =  tEpsCoefficients[eps];
					epsCoefficients.epsCoefficients[eps][0] *= tempVars.r1r0x;
					if(eps==0) epsCoefficients.epsCoefficients[eps][0] += 1;

					epsCoefficients.epsCoefficients[eps][1] =  tEpsCoefficients[eps];
					epsCoefficients.epsCoefficients[eps][1] *= tempVars.r1r0y;
					if(eps==1) epsCoefficients.epsCoefficients[eps][1] += 1;

					epsCoefficients.epsCoefficients[eps][2] =  tEpsCoefficients[eps];
					epsCoefficients.epsCoefficients[eps][2] *=  tempVars.r1r0z;
					if(eps==2) epsCoefficients.epsCoefficients[eps][2] += 1;
				} else { //the ray is in mesh 0 --> it does not have eps coefficients...
					//epsCoefficients.epsCoefficients[eps][0] =  tEpsCoefficients[eps]*tempVars.r1r0x;
					//epsCoefficients.epsCoefficients[eps][1] =  tEpsCoefficients[eps]*tempVars.r1r0y;
					//epsCoefficients.epsCoefficients[eps][2] =  tEpsCoefficients[eps]*tempVars.r1r0z;
					epsCoefficients.epsCoefficients[eps][0] =  tEpsCoefficients[eps];
					epsCoefficients.epsCoefficients[eps][1] =  tEpsCoefficients[eps];
					epsCoefficients.epsCoefficients[eps][2] =  tEpsCoefficients[eps];

					epsCoefficients.epsCoefficients[eps][0] *=  tempVars.r1r0x;
					epsCoefficients.epsCoefficients[eps][1] *=  tempVars.r1r0y;
					epsCoefficients.epsCoefficients[eps][2] *=  tempVars.r1r0z;
				}				
			}
			epsCoefficients.init=true;
		}
		omp_unset_lock(&epsCoefficients.lock);

	}

	return epsCoefficients.epsCoefficients[epsCoefficient-1][coord];
}

VertCoord SosPredicatesImpl::getEpsCoefficientsVertexFromIntersectionFaster(const VertexFromIntersection &v0, int epsCoefficient,int coord) const {
	assert(coord<=2);
	assert(epsCoefficient<=3);

	int meshId = v0.getMeshOfTriangleDefiningVertex();


	const VertCoord &v0t0x = getCoordinates( *(v0.triangle.getInputVertex(0)) )[0];
	const VertCoord &v0t1x = getCoordinates( *(v0.triangle.getInputVertex(1)) )[0];
	const VertCoord &v0t2x = getCoordinates( *(v0.triangle.getInputVertex(2)) )[0];
	const VertCoord &v0r0x = getCoordinates(v0.edge[0])[0];
	const VertCoord &v0r1x = getCoordinates(v0.edge[1])[0];
	const VertCoord &v0t0y = getCoordinates( *(v0.triangle.getInputVertex(0)) )[1];
	const VertCoord &v0t1y = getCoordinates( *(v0.triangle.getInputVertex(1)) )[1];
	const VertCoord &v0t2y = getCoordinates( *(v0.triangle.getInputVertex(2)) )[1];
	const VertCoord &v0r0y = getCoordinates(v0.edge[0])[1];
	const VertCoord &v0r1y = getCoordinates(v0.edge[1])[1];
	const VertCoord &v0t0z = getCoordinates( *(v0.triangle.getInputVertex(0)) )[2];
	const VertCoord &v0t1z = getCoordinates( *(v0.triangle.getInputVertex(1)) )[2];
	const VertCoord &v0t2z = getCoordinates( *(v0.triangle.getInputVertex(2)) )[2];
	const VertCoord &v0r0z = getCoordinates(v0.edge[0])[2];
	const VertCoord &v0r1z = getCoordinates(v0.edge[1])[2];

	if(epsCoefficient==0) return getCoordinates(v0)[coord];

	if(meshId==0) {
		if(coord==0) {
			switch(epsCoefficient) {
				case 1: return (v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + (v0r0y - v0r1y)*(v0t0z*v0t1x - v0t0x*v0t1z - v0t0z*v0t2x + v0t1z*v0t2x + v0t0x*v0t2z - v0t1x*v0t2z))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 2: return ((v0r0x - v0r1x)*(v0t0z*(v0t1x - v0t2x) + v0t1z*v0t2x - v0t1x*v0t2z + v0t0x*(-v0t1z + v0t2z)))/(-(v0r0y*v0t0z*v0t1x) + v0r1y*v0t0z*v0t1x + v0r0x*v0t0z*v0t1y - v0r1x*v0t0z*v0t1y + v0r0y*v0t0x*v0t1z - v0r1y*v0t0x*v0t1z - v0r0x*v0t0y*v0t1z + v0r1x*v0t0y*v0t1z + v0r0y*v0t0z*v0t2x - v0r1y*v0t0z*v0t2x - v0r0y*v0t1z*v0t2x + v0r1y*v0t1z*v0t2x - v0r0x*v0t0z*v0t2y + v0r1x*v0t0z*v0t2y + v0r0x*v0t1z*v0t2y - v0r1x*v0t1z*v0t2y + v0r0z*(-(v0t0x*v0t1y) + v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r1z*(-(v0t0y*v0t1x) + v0t0x*v0t1y + v0t0y*v0t2x - v0t1y*v0t2x - v0t0x*v0t2y + v0t1x*v0t2y) - v0r0y*v0t0x*v0t2z + v0r1y*v0t0x*v0t2z + v0r0x*v0t0y*v0t2z - v0r1x*v0t0y*v0t2z + v0r0y*v0t1x*v0t2z - v0r1y*v0t1x*v0t2z - v0r0x*v0t1y*v0t2z + v0r1x*v0t1y*v0t2z);
				case 3: return ((v0r0x - v0r1x)*(v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x - v0t1x*v0t2y + v0t0x*(-v0t1y + v0t2y)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
			}
		} else if (coord==1) {
			switch(epsCoefficient) {
				case 1: return  ((v0r0y - v0r1y)*(v0t0z*(v0t1y - v0t2y) + v0t1z*v0t2y - v0t1y*v0t2z + v0t0y*(-v0t1z + v0t2z)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 2: return  (v0r1z*(-(v0t0x*v0t1y) + v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + (v0r0x - v0r1x)*(v0t0y*v0t1z - v0t1z*v0t2y + v0t0z*(-v0t1y + v0t2y) - v0t0y*v0t2z + v0t1y*v0t2z))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 3: return  ((v0r0y - v0r1y)*(v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x - v0t1x*v0t2y + v0t0x*(-v0t1y + v0t2y)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
			}
		} else {
			assert(coord==2);	
			switch(epsCoefficient) {
				case 1: return  ((v0r0z - v0r1z)*(v0t0z*(v0t1y - v0t2y) + v0t1z*v0t2y - v0t1y*v0t2z + v0t0y*(-v0t1z + v0t2z)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 2: return  ((v0r0z - v0r1z)*(v0t0z*(v0t1x - v0t2x) + v0t1z*v0t2x - v0t1x*v0t2z + v0t0x*(-v0t1z + v0t2z)))/(-(v0r0y*v0t0z*v0t1x) + v0r1y*v0t0z*v0t1x + v0r0x*v0t0z*v0t1y - v0r1x*v0t0z*v0t1y + v0r0y*v0t0x*v0t1z - v0r1y*v0t0x*v0t1z - v0r0x*v0t0y*v0t1z + v0r1x*v0t0y*v0t1z + v0r0y*v0t0z*v0t2x - v0r1y*v0t0z*v0t2x - v0r0y*v0t1z*v0t2x + v0r1y*v0t1z*v0t2x - v0r0x*v0t0z*v0t2y + v0r1x*v0t0z*v0t2y + v0r0x*v0t1z*v0t2y - v0r1x*v0t1z*v0t2y + v0r0z*(-(v0t0x*v0t1y) + v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r1z*(-(v0t0y*v0t1x) + v0t0x*v0t1y + v0t0y*v0t2x - v0t1y*v0t2x - v0t0x*v0t2y + v0t1x*v0t2y) - v0r0y*v0t0x*v0t2z + v0r1y*v0t0x*v0t2z + v0r0x*v0t0y*v0t2z - v0r1x*v0t0y*v0t2z + v0r0y*v0t1x*v0t2z - v0r1y*v0t1x*v0t2z - v0r0x*v0t1y*v0t2z + v0r1x*v0t1y*v0t2z);
				case 3: return  (v0r0y*(-(v0t0x*v0t1z) + v0t0z*(v0t1x - v0t2x) + v0t1z*v0t2x + v0t0x*v0t2z - v0t1x*v0t2z) + v0r1y*(v0t0x*v0t1z - v0t1z*v0t2x + v0t0z*(-v0t1x + v0t2x) - v0t0x*v0t2z + v0t1x*v0t2z) + (v0r0x - v0r1x)*(v0t0y*v0t1z - v0t1z*v0t2y + v0t0z*(-v0t1y + v0t2y) - v0t0y*v0t2z + v0t1y*v0t2z))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
			}
		}
	} else {
		if(coord==0) {
			switch(epsCoefficient) {
				case 1: return ((v0r0x - v0r1x)*(v0t0z*(v0t1y - v0t2y) + v0t1z*v0t2y - v0t1y*v0t2z + v0t0y*(-v0t1z + v0t2z)))/(-(v0r0y*v0t0z*v0t1x) + v0r1y*v0t0z*v0t1x + v0r0x*v0t0z*v0t1y - v0r1x*v0t0z*v0t1y + v0r0y*v0t0x*v0t1z - v0r1y*v0t0x*v0t1z - v0r0x*v0t0y*v0t1z + v0r1x*v0t0y*v0t1z + v0r0y*v0t0z*v0t2x - v0r1y*v0t0z*v0t2x - v0r0y*v0t1z*v0t2x + v0r1y*v0t1z*v0t2x - v0r0x*v0t0z*v0t2y + v0r1x*v0t0z*v0t2y + v0r0x*v0t1z*v0t2y - v0r1x*v0t1z*v0t2y + v0r0z*(-(v0t0x*v0t1y) + v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r1z*(-(v0t0y*v0t1x) + v0t0x*v0t1y + v0t0y*v0t2x - v0t1y*v0t2x - v0t0x*v0t2y + v0t1x*v0t2y) - v0r0y*v0t0x*v0t2z + v0r1y*v0t0x*v0t2z + v0r0x*v0t0y*v0t2z - v0r1x*v0t0y*v0t2z + v0r0y*v0t1x*v0t2z - v0r1y*v0t1x*v0t2z - v0r0x*v0t1y*v0t2z + v0r1x*v0t1y*v0t2z);
				case 2: return  ((v0r0x - v0r1x)*(v0t0z*(v0t1x - v0t2x) + v0t1z*v0t2x - v0t1x*v0t2z + v0t0x*(-v0t1z + v0t2z)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 3: return  -(((v0r0x - v0r1x)*(v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x - v0t1x*v0t2y + v0t0x*(-v0t1y + v0t2y)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z));
			}
		} else if (coord==1) {	
			switch(epsCoefficient) {
				case 1: return -(((v0r0y - v0r1y)*(v0t0z*(v0t1y - v0t2y) + v0t1z*v0t2y - v0t1y*v0t2z + v0t0y*(-v0t1z + v0t2z)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z));
				case 2: return ((v0r0y - v0r1y)*(v0t0z*(v0t1x - v0t2x) + v0t1z*v0t2x - v0t1x*v0t2z + v0t0x*(-v0t1z + v0t2z)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 3: return -(((v0r0y - v0r1y)*(v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x - v0t1x*v0t2y + v0t0x*(-v0t1y + v0t2y)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z));
			}
		} else {
			assert(coord==2);	
			switch(epsCoefficient) {
				case 1: return ((v0r0z - v0r1z)*(v0t0z*(v0t1y - v0t2y) + v0t1z*v0t2y - v0t1y*v0t2z + v0t0y*(-v0t1z + v0t2z)))/(-(v0r0y*v0t0z*v0t1x) + v0r1y*v0t0z*v0t1x + v0r0x*v0t0z*v0t1y - v0r1x*v0t0z*v0t1y + v0r0y*v0t0x*v0t1z - v0r1y*v0t0x*v0t1z - v0r0x*v0t0y*v0t1z + v0r1x*v0t0y*v0t1z + v0r0y*v0t0z*v0t2x - v0r1y*v0t0z*v0t2x - v0r0y*v0t1z*v0t2x + v0r1y*v0t1z*v0t2x - v0r0x*v0t0z*v0t2y + v0r1x*v0t0z*v0t2y + v0r0x*v0t1z*v0t2y - v0r1x*v0t1z*v0t2y + v0r0z*(-(v0t0x*v0t1y) + v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r1z*(-(v0t0y*v0t1x) + v0t0x*v0t1y + v0t0y*v0t2x - v0t1y*v0t2x - v0t0x*v0t2y + v0t1x*v0t2y) - v0r0y*v0t0x*v0t2z + v0r1y*v0t0x*v0t2z + v0r0x*v0t0y*v0t2z - v0r1x*v0t0y*v0t2z + v0r0y*v0t1x*v0t2z - v0r1y*v0t1x*v0t2z - v0r0x*v0t1y*v0t2z + v0r1x*v0t1y*v0t2z);
				case 2: return ((v0r0z - v0r1z)*(v0t0z*(v0t1x - v0t2x) + v0t1z*v0t2x - v0t1x*v0t2z + v0t0x*(-v0t1z + v0t2z)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 3: return ((v0r0z - v0r1z)*(v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x - v0t1x*v0t2y + v0t0x*(-v0t1y + v0t2y)))/(-(v0r0y*v0t0z*v0t1x) + v0r1y*v0t0z*v0t1x + v0r0x*v0t0z*v0t1y - v0r1x*v0t0z*v0t1y + v0r0y*v0t0x*v0t1z - v0r1y*v0t0x*v0t1z - v0r0x*v0t0y*v0t1z + v0r1x*v0t0y*v0t1z + v0r0y*v0t0z*v0t2x - v0r1y*v0t0z*v0t2x - v0r0y*v0t1z*v0t2x + v0r1y*v0t1z*v0t2x - v0r0x*v0t0z*v0t2y + v0r1x*v0t0z*v0t2y + v0r0x*v0t1z*v0t2y - v0r1x*v0t1z*v0t2y + v0r0z*(-(v0t0x*v0t1y) + v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r1z*(-(v0t0y*v0t1x) + v0t0x*v0t1y + v0t0y*v0t2x - v0t1y*v0t2x - v0t0x*v0t2y + v0t1x*v0t2y) - v0r0y*v0t0x*v0t2z + v0r1y*v0t0x*v0t2z + v0r0x*v0t0y*v0t2z - v0r1x*v0t0y*v0t2z + v0r0y*v0t1x*v0t2z - v0r1y*v0t1x*v0t2z - v0r0x*v0t1y*v0t2z + v0r1x*v0t1y*v0t2z);
			}
		}
	}
}

const VertCoord& SosPredicatesImpl::getEpsCoefficientsVertexFromIntersection(const VertexFromIntersection &v0, int epsCoefficient,int coord) const {
	/*if(v0.getMeshOfTriangleDefiningVertex()==1) {
		assert(getEpsCoefficientsVertexFromIntersectionFaster2(v0,epsCoefficient,coord)==getEpsCoefficientsVertexFromIntersectionOriginal(v0,epsCoefficient,coord));
		//cerr << "When is 1: " << endl;
		//cerr << getEpsCoefficientsVertexFromIntersectionFaster2(v0,epsCoefficient,coord) << endl;
		//cerr << getEpsCoefficientsVertexFromIntersectionOriginal(v0,epsCoefficient,coord) << endl;
	} else {
		cerr << v0.getMeshOfTriangleDefiningVertex() << " " << epsCoefficient << " " << coord << endl;
		cerr << getEpsCoefficientsVertexFromIntersectionFaster2(v0,epsCoefficient,coord) << endl;
		cerr << getEpsCoefficientsVertexFromIntersectionOriginal(v0,epsCoefficient,coord) << endl;
		cerr << getEpsCoefficientsVertexFromIntersectionFaster2(v0,epsCoefficient,coord)-getEpsCoefficientsVertexFromIntersectionOriginal(v0,epsCoefficient,coord) << endl;
	}*/

/*	VertCoord ansFaster2 = getEpsCoefficientsVertexFromIntersectionFaster2(v0,epsCoefficient,coord);

	VertCoord ans = getEpsCoefficientsVertexFromIntersectionFaster(v0,epsCoefficient,coord);
	if(ansFaster2!=ans) {
		assert(false);

	}*/
	//if(epsCoefficient!=0) {
	//	cerr << v0.getMeshOfTriangleDefiningVertex() << " " << epsCoefficient << " " << coord << " " << ansFaster2-ans << endl;
	//}
	return getEpsCoefficientsVertexFromIntersectionFaster2(v0,epsCoefficient,coord);
}


//given a vertex from intersection, returns the epsilon coefficient of a given coordinate
//coord may be 0,1 or 2
//epsCoefficient may be 0, 1, 2 or 3
VertCoord SosPredicatesImpl::getEpsCoefficientsVertexFromIntersectionOriginal(const VertexFromIntersection &v0, int epsCoefficient,int coord) const {
	//return getEpsCoefficientsVertexFromIntersectionFaster(v0,epsCoefficient,coord);

	assert(coord<=2);
	assert(epsCoefficient<=3);

	int meshId = v0.getMeshOfTriangleDefiningVertex();

	const VertCoord &v0t0x = getCoordinates( *(v0.triangle.getInputVertex(0)) )[0];
	const VertCoord &v0t1x = getCoordinates( *(v0.triangle.getInputVertex(1)) )[0];
	const VertCoord &v0t2x = getCoordinates( *(v0.triangle.getInputVertex(2)) )[0];
	const VertCoord &v0r0x = getCoordinates(v0.edge[0])[0];
	const VertCoord &v0r1x = getCoordinates(v0.edge[1])[0];
	const VertCoord &v0t0y = getCoordinates( *(v0.triangle.getInputVertex(0)) )[1];
	const VertCoord &v0t1y = getCoordinates( *(v0.triangle.getInputVertex(1)) )[1];
	const VertCoord &v0t2y = getCoordinates( *(v0.triangle.getInputVertex(2)) )[1];
	const VertCoord &v0r0y = getCoordinates(v0.edge[0])[1];
	const VertCoord &v0r1y = getCoordinates(v0.edge[1])[1];
	const VertCoord &v0t0z = getCoordinates( *(v0.triangle.getInputVertex(0)) )[2];
	const VertCoord &v0t1z = getCoordinates( *(v0.triangle.getInputVertex(1)) )[2];
	const VertCoord &v0t2z = getCoordinates( *(v0.triangle.getInputVertex(2)) )[2];
	const VertCoord &v0r0z = getCoordinates(v0.edge[0])[2];
	const VertCoord &v0r1z = getCoordinates(v0.edge[1])[2];

	if(meshId==0) {
		if(coord==0) {
			switch(epsCoefficient) {
				case 0: return (v0r0z*v0r1x*(-(v0t1y*v0t2x) + v0t0y*(-v0t1x + v0t2x) + v0t0x*(v0t1y - v0t2y) + v0t1x*v0t2y) + v0r1x*(v0t0z*v0t1y*v0t2x - v0t0y*v0t1z*v0t2x - v0t0z*v0t1x*v0t2y + v0t0x*v0t1z*v0t2y + v0t0y*v0t1x*v0t2z - v0t0x*v0t1y*v0t2z + v0r0y*(-(v0t0x*v0t1z) + v0t0z*(v0t1x - v0t2x) + v0t1z*v0t2x + v0t0x*v0t2z - v0t1x*v0t2z)) + v0r0x*(-(v0t0z*v0t1y*v0t2x) + v0t0y*v0t1z*v0t2x + v0t0z*v0t1x*v0t2y - v0t0x*v0t1z*v0t2y + v0r1z*(-(v0t0x*v0t1y) + v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) - v0t0y*v0t1x*v0t2z + v0t0x*v0t1y*v0t2z + v0r1y*(-(v0t0z*v0t1x) + v0t0x*v0t1z + v0t0z*v0t2x - v0t1z*v0t2x - v0t0x*v0t2z + v0t1x*v0t2z)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 1: return (v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + (v0r0y - v0r1y)*(v0t0z*v0t1x - v0t0x*v0t1z - v0t0z*v0t2x + v0t1z*v0t2x + v0t0x*v0t2z - v0t1x*v0t2z))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 2: return ((v0r0x - v0r1x)*(v0t0z*(v0t1x - v0t2x) + v0t1z*v0t2x - v0t1x*v0t2z + v0t0x*(-v0t1z + v0t2z)))/(-(v0r0y*v0t0z*v0t1x) + v0r1y*v0t0z*v0t1x + v0r0x*v0t0z*v0t1y - v0r1x*v0t0z*v0t1y + v0r0y*v0t0x*v0t1z - v0r1y*v0t0x*v0t1z - v0r0x*v0t0y*v0t1z + v0r1x*v0t0y*v0t1z + v0r0y*v0t0z*v0t2x - v0r1y*v0t0z*v0t2x - v0r0y*v0t1z*v0t2x + v0r1y*v0t1z*v0t2x - v0r0x*v0t0z*v0t2y + v0r1x*v0t0z*v0t2y + v0r0x*v0t1z*v0t2y - v0r1x*v0t1z*v0t2y + v0r0z*(-(v0t0x*v0t1y) + v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r1z*(-(v0t0y*v0t1x) + v0t0x*v0t1y + v0t0y*v0t2x - v0t1y*v0t2x - v0t0x*v0t2y + v0t1x*v0t2y) - v0r0y*v0t0x*v0t2z + v0r1y*v0t0x*v0t2z + v0r0x*v0t0y*v0t2z - v0r1x*v0t0y*v0t2z + v0r0y*v0t1x*v0t2z - v0r1y*v0t1x*v0t2z - v0r0x*v0t1y*v0t2z + v0r1x*v0t1y*v0t2z);
				case 3: return ((v0r0x - v0r1x)*(v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x - v0t1x*v0t2y + v0t0x*(-v0t1y + v0t2y)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
			}
		} else if (coord==1) {
			switch(epsCoefficient) {
				case 0: return (v0r0z*v0r1y*(-(v0t1y*v0t2x) + v0t0y*(-v0t1x + v0t2x) + v0t0x*(v0t1y - v0t2y) + v0t1x*v0t2y) + v0r0y*(-(v0t0z*v0t1y*v0t2x) + v0t0y*v0t1z*v0t2x + v0t0z*v0t1x*v0t2y - v0t0x*v0t1z*v0t2y + v0r1z*(-(v0t0x*v0t1y) + v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) - v0t0y*v0t1x*v0t2z + v0t0x*v0t1y*v0t2z + v0r1x*(v0t0z*v0t1y - v0t0y*v0t1z - v0t0z*v0t2y + v0t1z*v0t2y + v0t0y*v0t2z - v0t1y*v0t2z)) + v0r1y*(v0t0z*v0t1y*v0t2x - v0t0y*v0t1z*v0t2x - v0t0z*v0t1x*v0t2y + v0t0x*v0t1z*v0t2y + v0t0y*v0t1x*v0t2z - v0t0x*v0t1y*v0t2z + v0r0x*(v0t0y*v0t1z - v0t1z*v0t2y + v0t0z*(-v0t1y + v0t2y) - v0t0y*v0t2z + v0t1y*v0t2z)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 1: return  ((v0r0y - v0r1y)*(v0t0z*(v0t1y - v0t2y) + v0t1z*v0t2y - v0t1y*v0t2z + v0t0y*(-v0t1z + v0t2z)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 2: return  (v0r1z*(-(v0t0x*v0t1y) + v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + (v0r0x - v0r1x)*(v0t0y*v0t1z - v0t1z*v0t2y + v0t0z*(-v0t1y + v0t2y) - v0t0y*v0t2z + v0t1y*v0t2z))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 3: return  ((v0r0y - v0r1y)*(v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x - v0t1x*v0t2y + v0t0x*(-v0t1y + v0t2y)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
			}
		} else {
			switch(epsCoefficient) {
  			case 0: return (v0r0z*(-(v0t0z*v0t1y*v0t2x) + v0t0y*v0t1z*v0t2x + v0t0z*v0t1x*v0t2y - v0t0x*v0t1z*v0t2y - v0t0y*v0t1x*v0t2z + v0t0x*v0t1y*v0t2z + v0r1y*(v0t0x*v0t1z - v0t1z*v0t2x + v0t0z*(-v0t1x + v0t2x) - v0t0x*v0t2z + v0t1x*v0t2z) + v0r1x*(v0t0z*v0t1y - v0t0y*v0t1z - v0t0z*v0t2y + v0t1z*v0t2y + v0t0y*v0t2z - v0t1y*v0t2z)) + v0r1z*(v0t0z*v0t1y*v0t2x - v0t0y*v0t1z*v0t2x - v0t0z*v0t1x*v0t2y + v0t0x*v0t1z*v0t2y + v0t0y*v0t1x*v0t2z - v0t0x*v0t1y*v0t2z + v0r0y*(-(v0t0x*v0t1z) + v0t0z*(v0t1x - v0t2x) + v0t1z*v0t2x + v0t0x*v0t2z - v0t1x*v0t2z) + v0r0x*(-(v0t0z*v0t1y) + v0t0y*v0t1z + v0t0z*v0t2y - v0t1z*v0t2y - v0t0y*v0t2z + v0t1y*v0t2z)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 1: return  ((v0r0z - v0r1z)*(v0t0z*(v0t1y - v0t2y) + v0t1z*v0t2y - v0t1y*v0t2z + v0t0y*(-v0t1z + v0t2z)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 2: return  ((v0r0z - v0r1z)*(v0t0z*(v0t1x - v0t2x) + v0t1z*v0t2x - v0t1x*v0t2z + v0t0x*(-v0t1z + v0t2z)))/(-(v0r0y*v0t0z*v0t1x) + v0r1y*v0t0z*v0t1x + v0r0x*v0t0z*v0t1y - v0r1x*v0t0z*v0t1y + v0r0y*v0t0x*v0t1z - v0r1y*v0t0x*v0t1z - v0r0x*v0t0y*v0t1z + v0r1x*v0t0y*v0t1z + v0r0y*v0t0z*v0t2x - v0r1y*v0t0z*v0t2x - v0r0y*v0t1z*v0t2x + v0r1y*v0t1z*v0t2x - v0r0x*v0t0z*v0t2y + v0r1x*v0t0z*v0t2y + v0r0x*v0t1z*v0t2y - v0r1x*v0t1z*v0t2y + v0r0z*(-(v0t0x*v0t1y) + v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r1z*(-(v0t0y*v0t1x) + v0t0x*v0t1y + v0t0y*v0t2x - v0t1y*v0t2x - v0t0x*v0t2y + v0t1x*v0t2y) - v0r0y*v0t0x*v0t2z + v0r1y*v0t0x*v0t2z + v0r0x*v0t0y*v0t2z - v0r1x*v0t0y*v0t2z + v0r0y*v0t1x*v0t2z - v0r1y*v0t1x*v0t2z - v0r0x*v0t1y*v0t2z + v0r1x*v0t1y*v0t2z);
				case 3: return  (v0r0y*(-(v0t0x*v0t1z) + v0t0z*(v0t1x - v0t2x) + v0t1z*v0t2x + v0t0x*v0t2z - v0t1x*v0t2z) + v0r1y*(v0t0x*v0t1z - v0t1z*v0t2x + v0t0z*(-v0t1x + v0t2x) - v0t0x*v0t2z + v0t1x*v0t2z) + (v0r0x - v0r1x)*(v0t0y*v0t1z - v0t1z*v0t2y + v0t0z*(-v0t1y + v0t2y) - v0t0y*v0t2z + v0t1y*v0t2z))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
			}
		}
	} else {
		if(coord==0) {
			switch(epsCoefficient) {
				case 0: return (v0r0z*v0r1x*(-(v0t1y*v0t2x) + v0t0y*(-v0t1x + v0t2x) + v0t0x*(v0t1y - v0t2y) + v0t1x*v0t2y) + v0r1x*(v0t0z*v0t1y*v0t2x - v0t0y*v0t1z*v0t2x - v0t0z*v0t1x*v0t2y + v0t0x*v0t1z*v0t2y + v0t0y*v0t1x*v0t2z - v0t0x*v0t1y*v0t2z + v0r0y*(-(v0t0x*v0t1z) + v0t0z*(v0t1x - v0t2x) + v0t1z*v0t2x + v0t0x*v0t2z - v0t1x*v0t2z)) + v0r0x*(-(v0t0z*v0t1y*v0t2x) + v0t0y*v0t1z*v0t2x + v0t0z*v0t1x*v0t2y - v0t0x*v0t1z*v0t2y + v0r1z*(-(v0t0x*v0t1y) + v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) - v0t0y*v0t1x*v0t2z + v0t0x*v0t1y*v0t2z + v0r1y*(-(v0t0z*v0t1x) + v0t0x*v0t1z + v0t0z*v0t2x - v0t1z*v0t2x - v0t0x*v0t2z + v0t1x*v0t2z)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 1: return ((v0r0x - v0r1x)*(v0t0z*(v0t1y - v0t2y) + v0t1z*v0t2y - v0t1y*v0t2z + v0t0y*(-v0t1z + v0t2z)))/(-(v0r0y*v0t0z*v0t1x) + v0r1y*v0t0z*v0t1x + v0r0x*v0t0z*v0t1y - v0r1x*v0t0z*v0t1y + v0r0y*v0t0x*v0t1z - v0r1y*v0t0x*v0t1z - v0r0x*v0t0y*v0t1z + v0r1x*v0t0y*v0t1z + v0r0y*v0t0z*v0t2x - v0r1y*v0t0z*v0t2x - v0r0y*v0t1z*v0t2x + v0r1y*v0t1z*v0t2x - v0r0x*v0t0z*v0t2y + v0r1x*v0t0z*v0t2y + v0r0x*v0t1z*v0t2y - v0r1x*v0t1z*v0t2y + v0r0z*(-(v0t0x*v0t1y) + v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r1z*(-(v0t0y*v0t1x) + v0t0x*v0t1y + v0t0y*v0t2x - v0t1y*v0t2x - v0t0x*v0t2y + v0t1x*v0t2y) - v0r0y*v0t0x*v0t2z + v0r1y*v0t0x*v0t2z + v0r0x*v0t0y*v0t2z - v0r1x*v0t0y*v0t2z + v0r0y*v0t1x*v0t2z - v0r1y*v0t1x*v0t2z - v0r0x*v0t1y*v0t2z + v0r1x*v0t1y*v0t2z);
				case 2: return  ((v0r0x - v0r1x)*(v0t0z*(v0t1x - v0t2x) + v0t1z*v0t2x - v0t1x*v0t2z + v0t0x*(-v0t1z + v0t2z)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 3: return  -(((v0r0x - v0r1x)*(v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x - v0t1x*v0t2y + v0t0x*(-v0t1y + v0t2y)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z));
			}
		} else if (coord==1) {	
			switch(epsCoefficient) {
				case 0: return (v0r0z*v0r1y*(-(v0t1y*v0t2x) + v0t0y*(-v0t1x + v0t2x) + v0t0x*(v0t1y - v0t2y) + v0t1x*v0t2y) + v0r0y*(-(v0t0z*v0t1y*v0t2x) + v0t0y*v0t1z*v0t2x + v0t0z*v0t1x*v0t2y - v0t0x*v0t1z*v0t2y + v0r1z*(-(v0t0x*v0t1y) + v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) - v0t0y*v0t1x*v0t2z + v0t0x*v0t1y*v0t2z + v0r1x*(v0t0z*v0t1y - v0t0y*v0t1z - v0t0z*v0t2y + v0t1z*v0t2y + v0t0y*v0t2z - v0t1y*v0t2z)) + v0r1y*(v0t0z*v0t1y*v0t2x - v0t0y*v0t1z*v0t2x - v0t0z*v0t1x*v0t2y + v0t0x*v0t1z*v0t2y + v0t0y*v0t1x*v0t2z - v0t0x*v0t1y*v0t2z + v0r0x*(v0t0y*v0t1z - v0t1z*v0t2y + v0t0z*(-v0t1y + v0t2y) - v0t0y*v0t2z + v0t1y*v0t2z)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 1: return -(((v0r0y - v0r1y)*(v0t0z*(v0t1y - v0t2y) + v0t1z*v0t2y - v0t1y*v0t2z + v0t0y*(-v0t1z + v0t2z)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z));
				case 2: return ((v0r0y - v0r1y)*(v0t0z*(v0t1x - v0t2x) + v0t1z*v0t2x - v0t1x*v0t2z + v0t0x*(-v0t1z + v0t2z)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 3: return -(((v0r0y - v0r1y)*(v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x - v0t1x*v0t2y + v0t0x*(-v0t1y + v0t2y)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z));
			}
		} else {
			assert(coord==2);	
			switch(epsCoefficient) {
				case 0: return (v0r0z*(-(v0t0z*v0t1y*v0t2x) + v0t0y*v0t1z*v0t2x + v0t0z*v0t1x*v0t2y - v0t0x*v0t1z*v0t2y - v0t0y*v0t1x*v0t2z + v0t0x*v0t1y*v0t2z + v0r1y*(v0t0x*v0t1z - v0t1z*v0t2x + v0t0z*(-v0t1x + v0t2x) - v0t0x*v0t2z + v0t1x*v0t2z) + v0r1x*(v0t0z*v0t1y - v0t0y*v0t1z - v0t0z*v0t2y + v0t1z*v0t2y + v0t0y*v0t2z - v0t1y*v0t2z)) + v0r1z*(v0t0z*v0t1y*v0t2x - v0t0y*v0t1z*v0t2x - v0t0z*v0t1x*v0t2y + v0t0x*v0t1z*v0t2y + v0t0y*v0t1x*v0t2z - v0t0x*v0t1y*v0t2z + v0r0y*(-(v0t0x*v0t1z) + v0t0z*(v0t1x - v0t2x) + v0t1z*v0t2x + v0t0x*v0t2z - v0t1x*v0t2z) + v0r0x*(-(v0t0z*v0t1y) + v0t0y*v0t1z + v0t0z*v0t2y - v0t1z*v0t2y - v0t0y*v0t2z + v0t1y*v0t2z)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 1: return ((v0r0z - v0r1z)*(v0t0z*(v0t1y - v0t2y) + v0t1z*v0t2y - v0t1y*v0t2z + v0t0y*(-v0t1z + v0t2z)))/(-(v0r0y*v0t0z*v0t1x) + v0r1y*v0t0z*v0t1x + v0r0x*v0t0z*v0t1y - v0r1x*v0t0z*v0t1y + v0r0y*v0t0x*v0t1z - v0r1y*v0t0x*v0t1z - v0r0x*v0t0y*v0t1z + v0r1x*v0t0y*v0t1z + v0r0y*v0t0z*v0t2x - v0r1y*v0t0z*v0t2x - v0r0y*v0t1z*v0t2x + v0r1y*v0t1z*v0t2x - v0r0x*v0t0z*v0t2y + v0r1x*v0t0z*v0t2y + v0r0x*v0t1z*v0t2y - v0r1x*v0t1z*v0t2y + v0r0z*(-(v0t0x*v0t1y) + v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r1z*(-(v0t0y*v0t1x) + v0t0x*v0t1y + v0t0y*v0t2x - v0t1y*v0t2x - v0t0x*v0t2y + v0t1x*v0t2y) - v0r0y*v0t0x*v0t2z + v0r1y*v0t0x*v0t2z + v0r0x*v0t0y*v0t2z - v0r1x*v0t0y*v0t2z + v0r0y*v0t1x*v0t2z - v0r1y*v0t1x*v0t2z - v0r0x*v0t1y*v0t2z + v0r1x*v0t1y*v0t2z);
				case 2: return ((v0r0z - v0r1z)*(v0t0z*(v0t1x - v0t2x) + v0t1z*v0t2x - v0t1x*v0t2z + v0t0x*(-v0t1z + v0t2z)))/(v0r0y*v0t0z*v0t1x - v0r1y*v0t0z*v0t1x - v0r0x*v0t0z*v0t1y + v0r1x*v0t0z*v0t1y - v0r0y*v0t0x*v0t1z + v0r1y*v0t0x*v0t1z + v0r0x*v0t0y*v0t1z - v0r1x*v0t0y*v0t1z - v0r0y*v0t0z*v0t2x + v0r1y*v0t0z*v0t2x + v0r0y*v0t1z*v0t2x - v0r1y*v0t1z*v0t2x + v0r0x*v0t0z*v0t2y - v0r1x*v0t0z*v0t2y - v0r0x*v0t1z*v0t2y + v0r1x*v0t1z*v0t2y + v0r1z*(v0t0y*v0t1x - v0t0x*v0t1y - v0t0y*v0t2x + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r0z*(v0t0x*v0t1y - v0t1y*v0t2x + v0t0y*(-v0t1x + v0t2x) - v0t0x*v0t2y + v0t1x*v0t2y) + v0r0y*v0t0x*v0t2z - v0r1y*v0t0x*v0t2z - v0r0x*v0t0y*v0t2z + v0r1x*v0t0y*v0t2z - v0r0y*v0t1x*v0t2z + v0r1y*v0t1x*v0t2z + v0r0x*v0t1y*v0t2z - v0r1x*v0t1y*v0t2z);
				case 3: return ((v0r0z - v0r1z)*(v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x - v0t1x*v0t2y + v0t0x*(-v0t1y + v0t2y)))/(-(v0r0y*v0t0z*v0t1x) + v0r1y*v0t0z*v0t1x + v0r0x*v0t0z*v0t1y - v0r1x*v0t0z*v0t1y + v0r0y*v0t0x*v0t1z - v0r1y*v0t0x*v0t1z - v0r0x*v0t0y*v0t1z + v0r1x*v0t0y*v0t1z + v0r0y*v0t0z*v0t2x - v0r1y*v0t0z*v0t2x - v0r0y*v0t1z*v0t2x + v0r1y*v0t1z*v0t2x - v0r0x*v0t0z*v0t2y + v0r1x*v0t0z*v0t2y + v0r0x*v0t1z*v0t2y - v0r1x*v0t1z*v0t2y + v0r0z*(-(v0t0x*v0t1y) + v0t0y*(v0t1x - v0t2x) + v0t1y*v0t2x + v0t0x*v0t2y - v0t1x*v0t2y) + v0r1z*(-(v0t0y*v0t1x) + v0t0x*v0t1y + v0t0y*v0t2x - v0t1y*v0t2x - v0t0x*v0t2y + v0t1x*v0t2y) - v0r0y*v0t0x*v0t2z + v0r1y*v0t0x*v0t2z + v0r0x*v0t0y*v0t2z - v0r1x*v0t0y*v0t2z + v0r0y*v0t1x*v0t2z - v0r1y*v0t1x*v0t2z - v0r0x*v0t1y*v0t2z + v0r1x*v0t1y*v0t2z);
			}
		}
	}
}
