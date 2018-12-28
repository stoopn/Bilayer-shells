/*
 *  BoundingBox.h
 *  newshell_prj
 *
 *  Created by Norbert Stoop on 10.05.09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef BoundingBox_H
#define BoundingBox_H

#include <vector_value.h>

class BoundingBox 
{
public:
	double mmin[3];
	double mmax[3];
	void update_from_triangle(const RealVectorValue& V1, const RealVectorValue& V2, const RealVectorValue& V3, const double thickness)
	{
		mmin[0] = std::min(V1(0),std::min(V2(0),V3(0))) - thickness;
		mmin[1] = std::min(V1(1),std::min(V2(1),V3(1))) - thickness;
		mmin[2] = std::min(V1(2),std::min(V2(2),V3(2))) - thickness;
		
		mmax[0] = std::max(V1(0),std::max(V2(0),V3(0))) + thickness;
		mmax[1] = std::max(V1(1),std::max(V2(1),V3(1))) + thickness;
		mmax[2] = std::max(V1(2),std::max(V2(2),V3(2))) + thickness;
	}
	void update_from_sphere(const RealVectorValue& P, double radius)
	{
		mmin[0] = P(0)-radius;
		mmin[1] = P(1)-radius;
		mmin[2] = P(2)-radius;
		
		mmax[0] = P(0)+radius;
		mmax[1] = P(1)+radius;
		mmax[2] = P(2)+radius;
		
	}
	bool test_overlap(const BoundingBox& other)
	{
		bool result = (mmin[0] > other.mmax[0]);
		result |= (mmax[0] < other.mmin[0]);
		result |= (mmin[1] > other.mmax[1]);
		result |= (mmax[1] < other.mmin[1]);
		result |= (mmin[2] > other.mmax[2]);
		result |= (mmax[2] < other.mmin[2]);
		return (!result);
	}
	void print_info() 
	{
		std::cout<<"mmin: "<<mmin[0]<<", "<<mmin[1]<<", "<<mmin[2]<<"\n";
		std::cout<<"mmax: "<<mmax[0]<<", "<<mmax[1]<<", "<<mmax[2]<<"\n";

	}
};


#endif
