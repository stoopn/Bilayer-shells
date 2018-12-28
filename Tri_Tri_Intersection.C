/*
 *  Tri_Tri_Intersection.h
 *  newshell_prj
 *
 *  Created by Norbert Stoop on 11.05.09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


#include "Tri_Tri_Intersection.h"

/* Triangle/triangle intersection test routine,
 * by Tomas Moller, 1997.
 * See article "A Fast Triangle-Triangle Intersection Test",
 * Journal of Graphics Tools, 2(2), 1997
 * updated: 2001-06-20 (added line of intersection)
 *
 * int tri_tri_intersect(float V0[3],float V1[3],float V2[3],
 *                       float U0[3],float U1[3],float U2[3])
 *
 * parameters: vertices of triangle 1: V0,V1,V2
 *             vertices of triangle 2: U0,U1,U2
 * result    : returns 1 if the triangles intersect, otherwise 0
 *
 * Here is a version withouts divisions (a little faster)
 * int NoDivTriTriIsect(float V0[3],float V1[3],float V2[3],
 *                      float U0[3],float U1[3],float U2[3]);
 *
 * This version computes the line of intersection as well (if they are not coplanar):
 * int tri_tri_intersect_with_isectline(float V0[3],float V1[3],float V2[3],
 *               float U0[3],float U1[3],float U2[3],int *coplanar,
 *               float isectpt1[3],float isectpt2[3]);
 * coplanar returns whether the tris are coplanar
 * isectpt1, isectpt2 are the endpoints of the line of intersection
 */
// this edge to edge test is based on Franlin Antonio's gem:
// "Faster Line Segment Intersection", in Graphics Gems III, pp. 199-202


//The edge notation is as follos:
/*   
         2
       /   \
      /     \
     /       \
    / 2     1 \
   /           \
  /      0      \
 0 ------------- 1
 
We can find out which node is at edge e1 and e2: e1+e2 can be either 2 (=> node 0), 1 (node 1) or 3 (node 2).
*/
int NEXT[3]={1,2,0};
int twoedges2node[3] = {1,0,2};


bool edge_edge_test(RealVectorValue const& V0, RealVectorValue const& U0, RealVectorValue const& U1,
					unsigned int i0, unsigned int i1, double Ax, double Ay)
{
	real_type Bx=U0(i0)-U1(i0);
	real_type By=U0(i1)-U1(i1);
	real_type Cx=V0(i0)-U0(i0);
	real_type Cy=V0(i1)-U0(i1);
	real_type f=Ay*Bx-Ax*By;
	real_type d=By*Cx-Bx*Cy;
	if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))
	{
		real_type e=Ax*Cy-Ay*Cx;
		if(f>0)
		{
			if(e>=0 && e<=f) return true;
		}
		else
		{
			if(e<=0 && e>=f) return true;
		}
	}
	else
		return false;
}

bool edge_against_tri_edges(RealVectorValue const& V0, RealVectorValue const& V1,
							RealVectorValue const& U0, RealVectorValue const& U1, RealVectorValue const& U2,
							unsigned int i0, unsigned int i1)
{
	real_type Ax=V1(i0)-V0(i0);
	real_type Ay=V1(i1)-V0(i1);
	// test edge U0,U1 against V0,V1
	if ( edge_edge_test(V0,U0,U1,i0,i1,Ax,Ay) ) return true;
	// test edge U1,U2 against V0,V1
	if ( edge_edge_test(V0,U1,U2,i0,i1,Ax,Ay) ) return true;
	// test edge U2,U1 against V0,V1
	if ( edge_edge_test(V0,U2,U0,i0,i1,Ax,Ay) ) return true;
	return false;
}




bool point_in_tri(RealVectorValue const& V0,
				  RealVectorValue const& U0, RealVectorValue const& U1, RealVectorValue const& U2,
				  unsigned int i0, unsigned int i1)
{
	// is T1 completly inside T2?
	// check if V0 is inside tri(U0,U1,U2)
	real_type a=U1(i1)-U0(i1);
	real_type b=-(U1(i0)-U0(i0));
	real_type c=-a*U0(i0)-b*U0(i1);
	real_type d0=a*V0(i0)+b*V0(i1)+c;
	
	a=U2(i1)-U1(i1);
	b=-(U2(i0)-U1(i0));
	c=-a*U1(i0)-b*U1(i1);
	real_type d1=a*V0(i0)+b*V0(i1)+c;
	
	a=U0(i1)-U2(i1);
	b=-(U0(i0)-U2(i0));
	c=-a*U2(i0)-b*U2(i1);
	real_type d2=a*V0(i0)+b*V0(i1)+c;
	if(d0*d1>0.0)
	{
		if(d0*d2>0.0) return true;
	}
	return false;
}


bool coplanar_tri_tri(RealVectorValue const& N,
					  RealVectorValue const& V0, RealVectorValue const& V1, RealVectorValue const& V2,
					  RealVectorValue const& U0, RealVectorValue const& U1, RealVectorValue const& U2)
{
	
	unsigned int i0,i1;
	// first project onto an axis-aligned plane, that maximizes the area
	// of the triangles, compute indices: i0,i1.
	real_type A0=std::fabs(N(0));
	real_type A1=std::fabs(N(1));
	real_type A2=std::fabs(N(2));
	if(A0>A1)
	{
		if(A0>A2)
		{
			i0=1;      // A[0] is greatest
			i1=2;
		}
		else
		{
			i0=0;      // A[2] is greatest
			i1=1;
		}
	}
	else   // A[0]<=A[1]
	{
		if(A2>A1)
		{
			i0=0;      // A[2] is greatest
			i1=1;
		}
		else
		{
			i0=0;      // A[1] is greatest
			i1=2;
		}
	}
	
	// test all edges of triangle 1 against the edges of triangle 2
	if( edge_against_tri_edges(V0,V1,U0,U1,U2,i0,i1) ) return true;
	if( edge_against_tri_edges(V1,V2,U0,U1,U2,i0,i1) ) return true;
	if( edge_against_tri_edges(V2,V0,U0,U1,U2,i0,i1) ) return true;
	
	// finally, test if tri1 is totally contained in tri2 or vice versa
	if( point_in_tri(V0,U0,U1,U2,i0,i1) ) return true;
	if( point_in_tri(U0,V0,V1,V2,i0,i1) ) return true;
	
	return false;
}


void isect2(RealVectorValue const& VTX0, RealVectorValue const& VTX1, RealVectorValue const& VTX2,
			real_type VV0, real_type VV1, real_type VV2,
			real_type D0, real_type D1, real_type D2,
			real_type & isect0, real_type & isect1)
{
	real_type tmp;
//	RealVectorValue diff;
	tmp=D0/(D0-D1);
	isect0 = VV0+(VV1-VV0)*tmp;
//	diff = VTX1-VTX2;
//	diff *= tmp;
//	isectpoint0 = diff+VTX0;
	tmp=D0/(D0-D2);
	isect1 = VV0+(VV2-VV0)*tmp;
//	diff = VTX2-VTX0;
//	diff *= tmp;
//	isectpoint1 = VTX0 + diff;
}


bool compute_intervals_isectline(RealVectorValue const& VERT0, RealVectorValue const& VERT1, RealVectorValue const& VERT2,
								 real_type VV0, real_type VV1, real_type VV2, real_type D0, real_type D1, real_type D2,
								 real_type D0D1, real_type D0D2,
								 real_type & isect0, real_type & isect1, int& edge0, int& edge1)
{
	if(D0D1>0.0f)
	{
		// here we know that D0D2<=0.0
		// that is D0, D1 are on the same side, D2 on the other or on the plane
		isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,isect0,isect1);
		edge0 = 2;
		edge1 = 1;
	}
	else if(D0D2>0.0f)
	{
		// here we know that d0d1<=0.0
		isect2(VERT1,VERT0,VERT2,VV1,VV0,VV2,D1,D0,D2,isect0,isect1);
		edge0 = 0;
		edge1 = 1;
	}
	else if(D1*D2>0.0f || D0!=0.0f)
	{
		// here we know that d0d1<=0.0 or that D0!=0.0
		isect2(VERT0,VERT1,VERT2,VV0,VV1,VV2,D0,D1,D2,isect0,isect1);
		edge0 = 0;
		edge1 = 2;
	}
	else if(D1!=0.0f)
	{
		isect2(VERT1,VERT0,VERT2,VV1,VV0,VV2,D1,D0,D2,isect0,isect1);
		edge0 = 0;
		edge1 = 1;
	}
	else if(D2!=0.0f)
	{
		isect2(VERT2,VERT0,VERT1,VV2,VV0,VV1,D2,D0,D1,isect0,isect1);
		edge0 = 2;
		edge1 = 1;
	}
	else
	{
		// triangles are coplanar
		return true;
	}
	return false;
}

void line_line_cpp(RealVectorValue const & V0, RealVectorValue const & V1,  RealVectorValue const & U0, RealVectorValue const & U1, RealVectorValue& cpt, double& s1, double &s2)
{
	
	RealVectorValue   u = V1-V0;
	RealVectorValue   v = U1 - U0;
	RealVectorValue   w = V0 - U0;
	double  a = u*u;        // always >= 0
	double  b = u*v;
	double  c = v*v;        // always >= 0
	double  d = u*w;
	double  e = v*w;
	double  D = a*c - b*b;       // always >= 0
	double  sc, sN, sD = D;      // sc = sN / sD, default sD = D >= 0
	double  tc, tN, tD = D;      // tc = tN / tD, default tD = D >= 0
	
	// compute the line parameters of the two closest points
	if (D < EPSILON) { // the lines are almost parallel
		sN = 0.0;        // force using point P0 on segment S1
		sD = 1.0;        // to prevent possible division by 0.0 later
		tN = e;
		tD = c;
	}
	else {                // get the closest points on the infinite lines
		sN = (b*e - c*d);
		tN = (a*e - b*d);
		if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
			sN = 0.0;
			tN = e;
			tD = c;
		}
		else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
			sN = sD;
			tN = e + b;
			tD = c;
		}
	}
	
	if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
		tN = 0.0;
		// recompute sc for this edge
		if (-d < 0.0)
			sN = 0.0;
		else if (-d > a)
			sN = sD;
		else {
			sN = -d;
			sD = a;
		}
	}
	else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
		tN = tD;
		// recompute sc for this edge
		if ((-d + b) < 0.0)
			sN = 0;
		else if ((-d + b) > a)
			sN = sD;
		else {
			sN = (-d + b);
			sD = a;
		}
	}
	// finally do the division to get sc and tc
	sc = (fabs(sN) < EPSILON ? 0.0 : sN / sD);
	tc = (fabs(tN) < EPSILON ? 0.0 : tN / tD);
	
	// get the difference of the two closest points
	cpt = w + (sc * u) - (tc * v);  // = S1(sc) - S2(tc)
	
	s1 = sc;
	s2 = tc;
	//s1t=sc;
	//s2t=tc;
}




bool triangle_triangle_interval_overlap( RealVectorValue const & V0, RealVectorValue const & V1, RealVectorValue const & V2
										, RealVectorValue const & U0, RealVectorValue const & U1, RealVectorValue const & U2
										, bool & coplanar, RealVectorValue & cpt, int& iedge1, int& iedge2, int& inode1, int& inode2
										, bool epsilon_test)
{
	
	// compute plane equation of triangle(V0,V1,V2)
	RealVectorValue E1 = V1-V0;
	RealVectorValue E2 = V2-V0;
	RealVectorValue N1 = E1.cross(E2);
	real_type d1 = -N1*V0;
	// plane equation 1: N1.X+d1=0
	
	// put U0,U1,U2 into plane equation 1 to compute signed distances to the plane
	real_type du0 = N1*U0+d1;
	real_type du1 = N1*U1+d1;
	real_type du2 = N1*U2+d1;
	
	// coplanarity robustness check
	if (epsilon_test)
	{
		if( std::fabs(du0) < EPSILON ) du0=0.0;
		if( std::fabs(du1) < EPSILON ) du1=0.0;
		if( std::fabs(du2) < EPSILON ) du2=0.0;
	}
	
	real_type du0du1=du0*du1;
	real_type du0du2=du0*du2;
	
	if(du0du1>0.0f && du0du2>0.0f) // same sign on all of them + not equal 0 ?
		return false;                // no intersection occurs
	
	// compute plane of triangle (U0,U1,U2)
	E1 = U1-U0;
	E2 = U2,U0;
	RealVectorValue N2 = E1.cross(E2);
	real_type d2 = -N2*U0;
	// plane equation 2: N2.X+d2=0
	
	// put V0,V1,V2 into plane equation 2
	real_type dv0 = N2*V0+d2;
	real_type dv1 = N2*V1+d2;
	real_type dv2 = N2*V2+d2;
	
	// coplanarity robustness check
	if (epsilon_test)
	{
		if( std::fabs(dv0) < EPSILON ) dv0=0.0;
		if( std::fabs(dv1) < EPSILON ) dv1=0.0;
		if( std::fabs(dv2) < EPSILON ) dv2=0.0;
	}
	
	real_type dv0dv1 = dv0*dv1;
	real_type dv0dv2 = dv0*dv2;
	
	if(dv0dv1>0.0f && dv0dv2>0.0f) // same sign on all of them + not equal 0 ?
		return false;                // no intersection occurs
	
	// compute direction of intersection line
	RealVectorValue D = N1.cross(N2);
	
	// compute and index to the largest component of D
	real_type max=std::fabs(D(0));
	unsigned int index=0;
	real_type b=std::fabs(D(1));
	real_type c=std::fabs(D(2));
	if(b>max) max=b,index=1;
	if(c>max) max=c,index=2;
	
	// this is the simplified projection onto L
	real_type vp0 = V0(index);
	real_type vp1 = V1(index);
	real_type vp2 = V2(index);
	
	real_type up0 = U0(index);
	real_type up1 = U1(index);
	real_type up2 = U2(index);
	
	real_type isect1[2], isect2[2];
//	RealVectorValue isectpointA1, isectpointA2;
	
	// compute interval for triangle 1
	

	int edge1[2];  //Gives us the two edges which could be intersecting for triangle 1
	int edge2[2];  //Same for triangle 2
	coplanar=compute_intervals_isectline(V0,V1,V2,vp0,vp1,vp2,dv0,dv1,dv2, 
										 dv0dv1,dv0dv2,isect1[0],isect1[1],edge1[0],edge1[1]);
	
	if(coplanar) return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);
	
	//RealVectorValue isectpointB1, isectpointB2;
	
	// compute interval for triangle 2
	compute_intervals_isectline(U0,U1,U2,up0,up1,up2,du0,du1,du2,
								du0du1,du0du2,isect2[0],isect2[1],edge2[0],edge2[1]);
	
	// sort so that a<=b
	unsigned int smallest1 = 0;
	if (isect1[0] > isect1[1])
	{
		std::swap(isect1[0], isect1[1]);
		smallest1=1;
	}
	unsigned int smallest2 = 0;
	if (isect2[0] > isect2[1])
	{
		std::swap(isect2[0], isect2[1]);
		smallest2=1;
	}
	
	if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return false;
	
	// at this point, we know that the triangles intersect	
	
	const RealVectorValue* nodes1[3]; nodes1[0] = &V0; nodes1[1] = &V1; nodes1[2] = &V2;
	const RealVectorValue* nodes2[3]; nodes2[0] = &U0; nodes2[1] = &U1; nodes2[2] = &U2;

	
	//This array tells us which node is shared by two edges. We use the above notation for edge numbering
	
	inode1=-1; //Flag that tells us if it is a node-face or edge-edge intersection
	inode2=-1;
	
	if(isect2[0]<isect1[0])
	{
		//isectpt1 = smallest1==0 ? isectpointA1 : isectpointA2;
		iedge1 = smallest1==0 ? edge1[0] : edge1[1];
		if(isect2[1]<isect1[1])
		{
		//	isectpt2 = smallest2==0 ? isectpointB2 : isectpointB1;
			iedge2 = smallest2==0 ? edge2[0] : edge2[1];
			//Now lets find the closest point projection from iedge1 on iedge2
			double s1,s2;
			line_line_cpp(*nodes1[iedge1], *nodes1[NEXT[iedge1]], *nodes2[iedge2], *nodes2[NEXT[iedge2]], cpt,s1,s2);
		}
		else
		{
		//	isectpt2 = smallest1==0 ? isectpointA2 : isectpointA1;
			//This means that both intersection points are from triangle 1. Ie. a node from triangle 1 is penetrating the
			//face of triangle 2. Let's find out which node it is:
			inode1 = twoedges2node[edge1[0]+edge1[1]-1];
			//Find the closest point project. cpt is the point on triangle 2 closest to the inode1
			if (inode1==0)
				cpt = V0 - dv0*N2/N2.size();
			else if (inode1==1)
				cpt = V1 - dv1*N2/N2.size();
			else
				cpt = V2 - dv2*N2/N2.size();
 
		}
	}
	else
	{
		//isectpt1 = smallest2==0 ? isectpointB1 : isectpointB2;
		iedge2 = smallest2==0 ? edge2[0] : edge2[1];
		if(isect2[1]>isect1[1])
		{
		//	isectpt2 = smallest1==0 ? isectpointA2 : isectpointA1;
			iedge1 = smallest1==0 ? edge1[0] : edge1[1];
			//Find the closest point projection between the two edges:
			double s1,s2;
			line_line_cpp(*nodes1[iedge1], *nodes1[NEXT[iedge1]], *nodes2[iedge2], *nodes2[NEXT[iedge2]], cpt,s1,s2);
			
		}
		else
		{
		//	isectpt2 = smallest2==0 ? isectpointB2 : isectpointB1;
		    //This means a node from triangle 2 is penetrating the
			//face of triangle 1. Let's find out which node it is:
			inode2 = twoedges2node[edge2[0]+edge2[1]-1];
			//Find the closest point project. cpt is the point on triangle 1 closest to the inode2
			if (inode2==0)
				cpt = U0 - du0*N1/N1.size();
			else if (inode1==1)
				cpt = U1 - du1*N1/N1.size();
			else
				cpt = U2 - du2*N1/N1.size();
		}
	}
	return true;
}


void cpt_point_triangle(RealVectorValue& P, RealVectorValue const & V0, RealVectorValue const & V1, RealVectorValue const & V2, 
						RealVectorValue& CPT, double& v, double &w)
{
	RealVectorValue ab = V1-V0;
	RealVectorValue ac = V2-V0;
	RealVectorValue ap =  P-V0;
	double d1 = ab*ap;
	double d2 = ac*ap;
	if (d1<= 0.0f && d2 <= 0.0f)
	{
		v=0; w=0;
		CPT = V0;
		return;
	}
	RealVectorValue bp = P-V1;
	double d3 = ab*bp;
	double d4 = ac*bp;
	if (d3 >= 0.0f && d4 <= d3) {
		v=1; w=0;
		CPT = V1;
		return;
	}
	double vc = d1*d4 - d3*d2;
	if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f)
	{
		v = d1 / (d1-d3);
		w = 0;
		CPT = V0+v*ab;
		return;
	}
	RealVectorValue cp = P - V2;
	double d5 = ab*cp;
	double d6 = ac*cp;
	if (d6 >= 0.0f && d5 <= d6)
	{
		v=0; w=1;
		CPT = V2; 
		return;
	}
	double vb = d5*d2 - d1*d6;
	if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f)
	{
		w = d2/(d2-d6);
		v=0.;
		CPT = V0 + w*ac;
		return;
	}
	double va = d3*d6 - d5*d4;
	if (va <= 0.0f && (d4-d3) >= 0.0f && (d5 - d6) >= 0.0f)
	{
		w = (d4 - d3)/((d4-d3) + (d5-d6));
		v = 1.-w;
		CPT = V1 + w*(V2-V1);
		return;
	}
	double denom = 1.0 / (va + vb + vc);
	v = vb*denom;
	w = vc*denom;
	CPT = V0 + ab*v + ac*w;
	return;
}
