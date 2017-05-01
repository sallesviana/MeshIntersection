#include "originalAlgFromMathematicaSosPredicatesImpl.h"


int OriginalAlgFromMathematicaSosPredicatesImpl::orientation1D_x(const InputVertex &v0, const InputVertex &v1) const {
	const int meshId0 = v0.getMeshId();
	const int meshId1 = v1.getMeshId();

	if(meshId0==0) {
		if(meshId1==0) {
			return orient1D_x_00(v0,v1);
		} else {
			return orient1D_x_01(v0,v1);
		}
	} else {
		if(meshId1==0) {
			return orient1D_x_10(v0,v1);
		} else {
			return orient1D_x_11(v0,v1);
		}
	}
}
int OriginalAlgFromMathematicaSosPredicatesImpl::orientation1D_x(const InputVertex &v0, const VertexFromIntersection &v1) const {
	const int meshId0 = v0.getMeshId();
	const int meshId1 = v1.getMeshOfTriangleDefiningVertex();

	if(meshId0==0) {
		if(meshId1==0) {
			return orient1D_x_00(v0,v1);
		} else {
			return orient1D_x_01(v0,v1);
		}
	} else {
		if(meshId1==0) {
			return orient1D_x_10(v0,v1);
		} else {
			return orient1D_x_11(v0,v1);
		}
	}
}
int OriginalAlgFromMathematicaSosPredicatesImpl::orientation1D_x(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const {
	const int meshId0 = v0.getMeshOfTriangleDefiningVertex();
	const int meshId1 = v1.getMeshOfTriangleDefiningVertex();

	if(meshId0==0) {
		if(meshId1==0) {
			return orient1D_x_00(v0,v1);
		} else {
			return orient1D_x_01(v0,v1);
		}
	} else {
		if(meshId1==0) {
			return orient1D_x_10(v0,v1);
		} else {
			return orient1D_x_11(v0,v1);
		}
	}
}





int OriginalAlgFromMathematicaSosPredicatesImpl::orientation1D_y(const InputVertex &v0, const InputVertex &v1) const {
	const int meshId0 = v0.getMeshId();
	const int meshId1 = v1.getMeshId();

	if(meshId0==0) {
		if(meshId1==0) {
			return orient1D_y_00(v0,v1);
		} else {
			return orient1D_y_01(v0,v1);
		}
	} else {
		if(meshId1==0) {
			return orient1D_y_10(v0,v1);
		} else {
			return orient1D_y_11(v0,v1);
		}
	}
}
int OriginalAlgFromMathematicaSosPredicatesImpl::orientation1D_y(const InputVertex &v0, const VertexFromIntersection &v1) const {
	const int meshId0 = v0.getMeshId();
	const int meshId1 = v1.getMeshOfTriangleDefiningVertex();

	if(meshId0==0) {
		if(meshId1==0) {
			return orient1D_y_00(v0,v1);
		} else {
			return orient1D_y_01(v0,v1);
		}
	} else {
		if(meshId1==0) {
			return orient1D_y_10(v0,v1);
		} else {
			return orient1D_y_11(v0,v1);
		}
	}
}
int OriginalAlgFromMathematicaSosPredicatesImpl::orientation1D_y(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const {
	const int meshId0 = v0.getMeshOfTriangleDefiningVertex();
	const int meshId1 = v1.getMeshOfTriangleDefiningVertex();

	if(meshId0==0) {
		if(meshId1==0) {
			return orient1D_y_00(v0,v1);
		} else {
			return orient1D_y_01(v0,v1);
		}
	} else {
		if(meshId1==0) {
			return orient1D_y_10(v0,v1);
		} else {
			return orient1D_y_11(v0,v1);
		}
	}
}


int OriginalAlgFromMathematicaSosPredicatesImpl::orientation1D_z(const InputVertex &v0, const InputVertex &v1) const {
	const int meshId0 = v0.getMeshId();
	const int meshId1 = v1.getMeshId();

	if(meshId0==0) {
		if(meshId1==0) {
			return orient1D_z_00(v0,v1);
		} else {
			return orient1D_z_01(v0,v1);
		}
	} else {
		if(meshId1==0) {
			return orient1D_z_10(v0,v1);
		} else {
			return orient1D_z_11(v0,v1);
		}
	}
}
int OriginalAlgFromMathematicaSosPredicatesImpl::orientation1D_z(const InputVertex &v0, const VertexFromIntersection &v1) const {
	const int meshId0 = v0.getMeshId();
	const int meshId1 = v1.getMeshOfTriangleDefiningVertex();

	if(meshId0==0) {
		if(meshId1==0) {
			return orient1D_z_00(v0,v1);
		} else {
			return orient1D_z_01(v0,v1);
		}
	} else {
		if(meshId1==0) {
			return orient1D_z_10(v0,v1);
		} else {
			return orient1D_z_11(v0,v1);
		}
	}
}
int OriginalAlgFromMathematicaSosPredicatesImpl::orientation1D_z(const VertexFromIntersection &v0, const VertexFromIntersection &v1) const {
	const int meshId0 = v0.getMeshOfTriangleDefiningVertex();
	const int meshId1 = v1.getMeshOfTriangleDefiningVertex();

	if(meshId0==0) {
		if(meshId1==0) {
			return orient1D_z_00(v0,v1);
		} else {
			return orient1D_z_01(v0,v1);
		}
	} else {
		if(meshId1==0) {
			return orient1D_z_10(v0,v1);
		} else {
			return orient1D_z_11(v0,v1);
		}
	}
}


/*****************************************************************************************************************************************************/
/*****************************************************************************************************************************************************/
/*****************************************************************************************************************************************************/
/*****************************************************************************************************************************************************/

int OriginalAlgFromMathematicaSosPredicatesImpl::orientation2D_x0(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const {
	const int meshId0 = v0.getMeshId();
	const int meshId1 = v1.getMeshId();
	const int meshId2 = v2.getMeshId();

	if(meshId0==0) {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_x0_000(v0,v1,v2);
			} else {
				return orient2D_x0_001(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_x0_010(v0,v1,v2);
			} else {
				return orient2D_x0_011(v0,v1,v2);
			}
		}
	} else {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_x0_100(v0,v1,v2);
			} else {
				return orient2D_x0_101(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_x0_110(v0,v1,v2);
			} else {
				return orient2D_x0_111(v0,v1,v2);
			}
		}
	}
}
int OriginalAlgFromMathematicaSosPredicatesImpl::orientation2D_x0(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const {
	const int meshId0 = v0.getMeshId();
	const int meshId1 = v1.getMeshId();
	const int meshId2 = v2.getMeshOfTriangleDefiningVertex();

	if(meshId0==0) {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_x0_000(v0,v1,v2);
			} else {
				return orient2D_x0_001(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_x0_010(v0,v1,v2);
			} else {
				return orient2D_x0_011(v0,v1,v2);
			}
		}
	} else {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_x0_100(v0,v1,v2);
			} else {
				return orient2D_x0_101(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_x0_110(v0,v1,v2);
			} else {
				return orient2D_x0_111(v0,v1,v2);
			}
		}
	}
}
int OriginalAlgFromMathematicaSosPredicatesImpl::orientation2D_x0(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const {
	const int meshId0 = v0.getMeshId();
	const int meshId1 = v1.getMeshOfTriangleDefiningVertex();
	const int meshId2 = v2.getMeshOfTriangleDefiningVertex();

	if(meshId0==0) {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_x0_000(v0,v1,v2);
			} else {
				return orient2D_x0_001(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_x0_010(v0,v1,v2);
			} else {
				return orient2D_x0_011(v0,v1,v2);
			}
		}
	} else {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_x0_100(v0,v1,v2);
			} else {
				return orient2D_x0_101(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_x0_110(v0,v1,v2);
			} else {
				return orient2D_x0_111(v0,v1,v2);
			}
		}
	}
}
int OriginalAlgFromMathematicaSosPredicatesImpl::orientation2D_x0(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const {
	const int meshId0 = v0.getMeshOfTriangleDefiningVertex();
	const int meshId1 = v1.getMeshOfTriangleDefiningVertex();
	const int meshId2 = v2.getMeshOfTriangleDefiningVertex();

	if(meshId0==0) {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_x0_000(v0,v1,v2);
			} else {
				return orient2D_x0_001(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_x0_010(v0,v1,v2);
			} else {
				return orient2D_x0_011(v0,v1,v2);
			}
		}
	} else {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_x0_100(v0,v1,v2);
			} else {
				return orient2D_x0_101(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_x0_110(v0,v1,v2);
			} else {
				return orient2D_x0_111(v0,v1,v2);
			}
		}
	}
}



int OriginalAlgFromMathematicaSosPredicatesImpl::orientation2D_y0(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const {
	const int meshId0 = v0.getMeshId();
	const int meshId1 = v1.getMeshId();
	const int meshId2 = v2.getMeshId();

	if(meshId0==0) {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_y0_000(v0,v1,v2);
			} else {
				return orient2D_y0_001(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_y0_010(v0,v1,v2);
			} else {
				return orient2D_y0_011(v0,v1,v2);
			}
		}
	} else {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_y0_100(v0,v1,v2);
			} else {
				return orient2D_y0_101(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_y0_110(v0,v1,v2);
			} else {
				return orient2D_y0_111(v0,v1,v2);
			}
		}
	}
}
int OriginalAlgFromMathematicaSosPredicatesImpl::orientation2D_y0(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const {
	const int meshId0 = v0.getMeshId();
	const int meshId1 = v1.getMeshId();
	const int meshId2 = v2.getMeshOfTriangleDefiningVertex();

	if(meshId0==0) {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_y0_000(v0,v1,v2);
			} else {
				return orient2D_y0_001(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_y0_010(v0,v1,v2);
			} else {
				return orient2D_y0_011(v0,v1,v2);
			}
		}
	} else {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_y0_100(v0,v1,v2);
			} else {
				return orient2D_y0_101(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_y0_110(v0,v1,v2);
			} else {
				return orient2D_y0_111(v0,v1,v2);
			}
		}
	}
}
int OriginalAlgFromMathematicaSosPredicatesImpl::orientation2D_y0(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const {
	const int meshId0 = v0.getMeshId();
	const int meshId1 = v1.getMeshOfTriangleDefiningVertex();
	const int meshId2 = v2.getMeshOfTriangleDefiningVertex();

	if(meshId0==0) {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_y0_000(v0,v1,v2);
			} else {
				return orient2D_y0_001(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_y0_010(v0,v1,v2);
			} else {
				return orient2D_y0_011(v0,v1,v2);
			}
		}
	} else {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_y0_100(v0,v1,v2);
			} else {
				return orient2D_y0_101(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_y0_110(v0,v1,v2);
			} else {
				return orient2D_y0_111(v0,v1,v2);
			}
		}
	}
}
int OriginalAlgFromMathematicaSosPredicatesImpl::orientation2D_y0(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const {
	const int meshId0 = v0.getMeshOfTriangleDefiningVertex();
	const int meshId1 = v1.getMeshOfTriangleDefiningVertex();
	const int meshId2 = v2.getMeshOfTriangleDefiningVertex();

	if(meshId0==0) {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_y0_000(v0,v1,v2);
			} else {
				return orient2D_y0_001(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_y0_010(v0,v1,v2);
			} else {
				return orient2D_y0_011(v0,v1,v2);
			}
		}
	} else {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_y0_100(v0,v1,v2);
			} else {
				return orient2D_y0_101(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_y0_110(v0,v1,v2);
			} else {
				return orient2D_y0_111(v0,v1,v2);
			}
		}
	}
}


int OriginalAlgFromMathematicaSosPredicatesImpl::orientation2D_z0(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2) const {
	const int meshId0 = v0.getMeshId();
	const int meshId1 = v1.getMeshId();
	const int meshId2 = v2.getMeshId();

	if(meshId0==0) {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_z0_000(v0,v1,v2);
			} else {
				return orient2D_z0_001(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_z0_010(v0,v1,v2);
			} else {
				return orient2D_z0_011(v0,v1,v2);
			}
		}
	} else {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_z0_100(v0,v1,v2);
			} else {
				return orient2D_z0_101(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_z0_110(v0,v1,v2);
			} else {
				return orient2D_z0_111(v0,v1,v2);
			}
		}
	}
}
int OriginalAlgFromMathematicaSosPredicatesImpl::orientation2D_z0(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2) const {
	const int meshId0 = v0.getMeshId();
	const int meshId1 = v1.getMeshId();
	const int meshId2 = v2.getMeshOfTriangleDefiningVertex();

	if(meshId0==0) {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_z0_000(v0,v1,v2);
			} else {
				return orient2D_z0_001(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_z0_010(v0,v1,v2);
			} else {
				return orient2D_z0_011(v0,v1,v2);
			}
		}
	} else {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_z0_100(v0,v1,v2);
			} else {
				return orient2D_z0_101(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_z0_110(v0,v1,v2);
			} else {
				return orient2D_z0_111(v0,v1,v2);
			}
		}
	}
}
int OriginalAlgFromMathematicaSosPredicatesImpl::orientation2D_z0(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const {
	const int meshId0 = v0.getMeshId();
	const int meshId1 = v1.getMeshOfTriangleDefiningVertex();
	const int meshId2 = v2.getMeshOfTriangleDefiningVertex();

	if(meshId0==0) {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_z0_000(v0,v1,v2);
			} else {
				return orient2D_z0_001(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_z0_010(v0,v1,v2);
			} else {
				return orient2D_z0_011(v0,v1,v2);
			}
		}
	} else {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_z0_100(v0,v1,v2);
			} else {
				return orient2D_z0_101(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_z0_110(v0,v1,v2);
			} else {
				return orient2D_z0_111(v0,v1,v2);
			}
		}
	}
}
int OriginalAlgFromMathematicaSosPredicatesImpl::orientation2D_z0(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2) const {
	const int meshId0 = v0.getMeshOfTriangleDefiningVertex();
	const int meshId1 = v1.getMeshOfTriangleDefiningVertex();
	const int meshId2 = v2.getMeshOfTriangleDefiningVertex();

	if(meshId0==0) {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_z0_000(v0,v1,v2);
			} else {
				return orient2D_z0_001(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_z0_010(v0,v1,v2);
			} else {
				return orient2D_z0_011(v0,v1,v2);
			}
		}
	} else {
		if(meshId1==0) {
			if(meshId2==0) {
				return orient2D_z0_100(v0,v1,v2);
			} else {
				return orient2D_z0_101(v0,v1,v2);
			}
		} else {
			if(meshId2==0) {
				return orient2D_z0_110(v0,v1,v2);
			} else {
				return orient2D_z0_111(v0,v1,v2);
			}
		}
	}
}


/*****************************************************************************************************************************************************/
/*****************************************************************************************************************************************************/
/*****************************************************************************************************************************************************/
int OriginalAlgFromMathematicaSosPredicatesImpl::orientation1D(const InputVertex &v0, const InputVertex &v1,int whatAxisProjectTo) const {
	switch(whatAxisProjectTo) {
		case 0:  return orientation1D_x(v0,v1); break;
		case 1:  return orientation1D_y(v0,v1); break;
		case 2:  return orientation1D_z(v0,v1); break;
	}
	assert(false);
}

int OriginalAlgFromMathematicaSosPredicatesImpl::orientation1D(const InputVertex &v0, const VertexFromIntersection &v1,int whatAxisProjectTo) const {
	switch(whatAxisProjectTo) {
		case 0:  return orientation1D_x(v0,v1); break;
		case 1:  return orientation1D_y(v0,v1); break;
		case 2:  return orientation1D_z(v0,v1); break;
	}
	assert(false);
}

int OriginalAlgFromMathematicaSosPredicatesImpl::orientation1D(const VertexFromIntersection &v0, const VertexFromIntersection &v1,int whatAxisProjectTo) const {
	switch(whatAxisProjectTo) {
		case 0:  return orientation1D_x(v0,v1); break;
		case 1:  return orientation1D_y(v0,v1); break;
		case 2:  return orientation1D_z(v0,v1); break;
	}
	assert(false);
}


int OriginalAlgFromMathematicaSosPredicatesImpl::orientation2D(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2,int whatPlaneProjectTo) const{
	switch(whatPlaneProjectTo) {
		case 0:  return orientation2D_x0(v0,v1,v2); break;
		case 1:  return orientation2D_y0(v0,v1,v2); break;
		case 2:  return orientation2D_z0(v0,v1,v2); break;
	}
	assert(false);
}

int OriginalAlgFromMathematicaSosPredicatesImpl::orientation2D(const InputVertex &v0, const InputVertex &v1, const VertexFromIntersection &v2,int whatPlaneProjectTo) const{
	switch(whatPlaneProjectTo) {
		case 0:  return orientation2D_x0(v0,v1,v2); break;
		case 1:  return orientation2D_y0(v0,v1,v2); break;
		case 2:  return orientation2D_z0(v0,v1,v2); break;
	}
	assert(false);
}

int OriginalAlgFromMathematicaSosPredicatesImpl::orientation2D(const InputVertex &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2,int whatPlaneProjectTo) const{
	switch(whatPlaneProjectTo) {
		case 0:  return orientation2D_x0(v0,v1,v2); break;
		case 1:  return orientation2D_y0(v0,v1,v2); break;
		case 2:  return orientation2D_z0(v0,v1,v2); break;
	}
	assert(false);
}

int OriginalAlgFromMathematicaSosPredicatesImpl::orientation2D(const VertexFromIntersection &v0, const VertexFromIntersection &v1, const VertexFromIntersection &v2,int whatPlaneProjectTo) const{
	switch(whatPlaneProjectTo) {
		case 0:  return orientation2D_x0(v0,v1,v2); break;
		case 1:  return orientation2D_y0(v0,v1,v2); break;
		case 2:  return orientation2D_z0(v0,v1,v2); break;
	}
	assert(false);
}

/*****************************************************************************************************************************************************/

int OriginalAlgFromMathematicaSosPredicatesImpl::orientation3D(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const InputVertex &v3) const{
	const int meshId0 = v0.getMeshId();
	const int meshId1 = v1.getMeshId();
	const int meshId2 = v2.getMeshId();
	const int meshId3 = v3.getMeshId();

	if(meshId0==0) {
		if(meshId1==0) {
			if(meshId2==0) {
				if(meshId3==0) {
					return orient3D_0000(v0,v1,v2,v3);
				} else {
					return orient3D_0001(v0,v1,v2,v3);
				}
			} else { //v2 is 1...
				if(meshId3==0) {
					return orient3D_0010(v0,v1,v2,v3);
				} else {
					return orient3D_0011(v0,v1,v2,v3);
				}
			}
		} else { //v1 is 1
			if(meshId2==0) {
				if(meshId3==0) {
					return orient3D_0100(v0,v1,v2,v3);
				} else {
					return orient3D_0101(v0,v1,v2,v3);
				}
			} else { //v2 is 1...
				if(meshId3==0) {
					return orient3D_0110(v0,v1,v2,v3);
				} else {
					return orient3D_0111(v0,v1,v2,v3);
				}
			}
		}
	} else {
		if(meshId1==0) {
			if(meshId2==0) {
				if(meshId3==0) {
					return orient3D_1000(v0,v1,v2,v3);
				} else {
					return orient3D_1001(v0,v1,v2,v3);
				}
			} else { //v2 is 1...
				if(meshId3==0) {
					return orient3D_1010(v0,v1,v2,v3);
				} else {
					return orient3D_1011(v0,v1,v2,v3);
				}
			}
		} else { //v1 is 1
			if(meshId2==0) {
				if(meshId3==0) {
					return orient3D_1100(v0,v1,v2,v3);
				} else {
					return orient3D_1101(v0,v1,v2,v3);
				}
			} else { //v2 is 1...
				if(meshId3==0) {
					return orient3D_1110(v0,v1,v2,v3);
				} else {
					return orient3D_1111(v0,v1,v2,v3);
				}
			}
		}
	}
}

int OriginalAlgFromMathematicaSosPredicatesImpl::orientation3D(const InputVertex &v0, const InputVertex &v1, const InputVertex &v2, const VertexFromIntersection &v3) const{
	const int meshId0 = v0.getMeshId();
	const int meshId1 = v1.getMeshId();
	const int meshId2 = v2.getMeshId();
	const int meshId3 = v3.getMeshOfTriangleDefiningVertex();

	if(meshId0==0) {
		if(meshId1==0) {
			if(meshId2==0) {
				if(meshId3==0) {
					return orient3D_0000(v0,v1,v2,v3);
				} else {
					return orient3D_0001(v0,v1,v2,v3);
				}
			} else { //v2 is 1...
				if(meshId3==0) {
					return orient3D_0010(v0,v1,v2,v3);
				} else {
					return orient3D_0011(v0,v1,v2,v3);
				}
			}
		} else { //v1 is 1
			if(meshId2==0) {
				if(meshId3==0) {
					return orient3D_0100(v0,v1,v2,v3);
				} else {
					return orient3D_0101(v0,v1,v2,v3);
				}
			} else { //v2 is 1...
				if(meshId3==0) {
					return orient3D_0110(v0,v1,v2,v3);
				} else {
					return orient3D_0111(v0,v1,v2,v3);
				}
			}
		}
	} else {
		if(meshId1==0) {
			if(meshId2==0) {
				if(meshId3==0) {
					return orient3D_1000(v0,v1,v2,v3);
				} else {
					return orient3D_1001(v0,v1,v2,v3);
				}
			} else { //v2 is 1...
				if(meshId3==0) {
					return orient3D_1010(v0,v1,v2,v3);
				} else {
					return orient3D_1011(v0,v1,v2,v3);
				}
			}
		} else { //v1 is 1
			if(meshId2==0) {
				if(meshId3==0) {
					return orient3D_1100(v0,v1,v2,v3);
				} else {
					return orient3D_1101(v0,v1,v2,v3);
				}
			} else { //v2 is 1...
				if(meshId3==0) {
					return orient3D_1110(v0,v1,v2,v3);
				} else {
					return orient3D_1111(v0,v1,v2,v3);
				}
			}
		}
	}
}

