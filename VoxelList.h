/*
 *  VoxelList.h
 *  newshell_prj
 *
 *  Created by Norbert Stoop on 10.05.09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef VOXELLIST_H
#define VOXELLIST_H

#include "shellelement.h"
#include "BoundingBox.h"
#include "Tri_Tri_Intersection.h"



class SphereNode
{
public:
  SphereNode() {
	node=NULL; 
	VoxelIndices.clear(); 
	is_ghost=true; 
	normalct=1; 
	normal.zero(); 
	exclude=false; 
  }
	Node* node;
	BoundingBox Bbox;
	std::vector<int> VoxelIndices;
	RealVectorValue current_position;
	RealVectorValue current_velocity;
	bool is_ghost;
	RealVectorValue normal;
	int normalct;
	bool exclude;
};

class TriContact
{
public:
	const ShellCacheObj* obj1;
	const ShellCacheObj* obj2;
	
	double v1,w1,v2,w2,delta;
	
	RealVectorValue dist;
	double dist2;
	
	void print_info() {
		std::cout<<"Contact id1:"<<obj1->elem->id()<<", id2: "<<obj2->elem->id()<<", d "<<dist<<", d^2: "<<dist2
		<<", (v1,w1,v2,w2)={"<<v1<<w1<<v2<<w2<<"}"<<std::endl;
	}
};





class VoxelCell
{
public:
	VoxelCell() {
		
	}
	~VoxelCell() { 
		ShellObjects.clear(); 
	}
	void clear_objects() {
		ShellObjects.clear(); 
	}
	void add(ShellCacheObj* obj) { 
		ShellObjects.push_back(obj); 
	}	
	int n_objects() {
		return ShellObjects.size();
	}
	ShellCacheObj* get_object(int i) {
		return ShellObjects[i]; 
	}
private:
	std::vector<ShellCacheObj*> ShellObjects;
};




class VoxelList {
	
public:
	VoxelList();

	SphereNode* get_node(int i) { return &obj_voxel_map[i]; }
	
	void clear();
	void add(ShellCacheObj* selem);
	void update(SphereNode* obj);
	void initialize(RealVectorValue P1, RealVectorValue P2, int Nv, double thickness, int Nelements, int Nnodes, std::vector<ShellCacheObj> *scl);
	void update_tri_tri_contacts();
	void update_surf_surf_contacts();
	
	
	int get_n_voxels() { return Voxels.size(); }
	int get_n_nodes() { return obj_voxel_map.size(); }
	int get_n_elements() { return n_elements; }
	
	void print_info();
	void clear_contacts();
	

	
	//check if a ridge is formed
	bool ridge_check(const ShellCacheObj* obj1, const ShellCacheObj* obj2);
	
	
	bool edge_edge_test(ShellCacheObj* obj1, ShellCacheObj* obj2);
	bool tri_tri_test(ShellCacheObj* obj1, ShellCacheObj* obj2);
	
	//vector for all current contacts
	std::vector<TriContact> Contacts;

	std::ofstream cfile;
	std::string filename;
		
private:
	int Nx, Ny, Nz; //Number of voxels in each direction
	double rc;
	double thickness_sq;
	double thickness;
	
	
	int n_elements;
	int n_nodes;
	int boxshift; //For coordinate shift if the voxel cube is not starting at (0,0,0)
		
	std::vector<VoxelCell> Voxels; 
	std::vector<SphereNode> obj_voxel_map;
	std::vector<ShellCacheObj> *ShellCachePtr;
};
#endif
