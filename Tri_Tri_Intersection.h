/*
 *  Tri_Tri_Intersection.h
 *  newshell_prj
 *
 *  Created by Norbert Stoop on 11.05.09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TRI_TRI_INTERSECTION_H
#define TRI_TRI_INTERSECTION_H

#include <vector_value.h>
#include <cmath>

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
#define EPSILON 1e-20
typedef double real_type;

bool edge_edge_test(RealVectorValue const& V0, RealVectorValue const& U0, RealVectorValue const& U1,
					unsigned int i0, unsigned int i1, double Ax, double Ay);


bool edge_against_tri_edges(RealVectorValue const& V0, RealVectorValue const& V1,
							RealVectorValue const& U0, RealVectorValue const& U1, RealVectorValue const& U2,
							unsigned int i0, unsigned int i1);

// this edge to edge test is based on Franlin Antonio's gem:
// "Faster Line Segment Intersection", in Graphics Gems III, pp. 199-202
bool edge_edge_test(RealVectorValue const& V0, RealVectorValue const& U0, RealVectorValue const& U1,
					unsigned int i0, unsigned int i1, double Ax, double Ay);

bool point_in_tri(RealVectorValue const& V0,
				  RealVectorValue const& U0, RealVectorValue const& U1, RealVectorValue const& U2,
				  unsigned int i0, unsigned int i1);


bool coplanar_tri_tri(RealVectorValue const& N,
					  RealVectorValue const& V0, RealVectorValue const& V1, RealVectorValue const& V2,
					  RealVectorValue const& U0, RealVectorValue const& U1, RealVectorValue const& U2);


void isect2(RealVectorValue const& VTX0, RealVectorValue const& VTX1, RealVectorValue const& VTX2,
			real_type VV0, real_type VV1, real_type VV2,
			real_type D0, real_type D1, real_type D2,
			real_type & isect0, real_type & isect1, RealVectorValue & isectpoint0, RealVectorValue & isectpoint1);


bool compute_intervals_isectline(RealVectorValue const& VERT0, RealVectorValue const& VERT1, RealVectorValue const& VERT2,
								 real_type VV0, real_type VV1, real_type VV2, real_type D0, real_type D1, real_type D2,
								 real_type D0D1, real_type D0D2,
								 real_type & isect0, real_type & isect1,
								 RealVectorValue & isectpoint0, RealVectorValue & isectpoint1);


bool triangle_triangle_interval_overlap( RealVectorValue const & V0, RealVectorValue const & V1, RealVectorValue const & V2
										, RealVectorValue const & U0, RealVectorValue const & U1, RealVectorValue const & U2
										, bool & coplanar        , RealVectorValue & cpt, int& iedge1, int& iedge2, int& inode1, int& inode2
										, bool epsilon_test = true );

void line_line_cpp(RealVectorValue const & V0, RealVectorValue const & V1,  RealVectorValue const & U0, RealVectorValue const & U1, 
				   RealVectorValue& cpt, double &s1, double &s2);

void cpt_point_triangle(RealVectorValue& P, RealVectorValue const & V0, RealVectorValue const & V1, RealVectorValue const & V2, 
						RealVectorValue& CPT, double& v, double &w);
#endif
