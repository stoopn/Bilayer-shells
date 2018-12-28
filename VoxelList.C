/*
 *  VoxelList.h
 *  newshell_prj
 *
 *  Created by Norbert Stoop on 10.05.09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include "VoxelList.h"
// #include "inversion4x4Matrix.h"
#include <cmath>
#include <set>


void VoxelList::clear()
{
	for (int v=0; v<Voxels.size(); v++) {
		Voxels[v].clear_objects();
	}
}

VoxelList::VoxelList()
	{
	filename="Contacts_w.dat";
	cfile.open(filename.c_str());
	}

void VoxelList::add(ShellCacheObj* obj)
{
	int sidx;
	int nx=0;
	int ny=0;
	int nz=0;
	int idx[2], idy[2], idz[2];
	obj->VoxelIndices.clear();

	//Check if minima/maxima of bbox in each axis direction are different -> tells us if bbox overlaps multiple voxels
	idx[0] = obj->Bbox.mmin[0]/rc;
	idx[1] =  obj->Bbox.mmax[0]/rc;
	if (idx[0]<idx[1])
		nx=idx[1]-idx[0];
	idz[0] =  obj->Bbox.mmin[2]/rc;
	idz[1] =  obj->Bbox.mmax[2]/rc;
	if (idz[0]<idz[1])
		nz=idz[1]-idz[0];
	idy[0] =  obj->Bbox.mmin[1]/rc;
	idy[1] =  obj->Bbox.mmax[1]/rc;
	if (idy[0]<idy[1])
		ny=idy[1]-idy[0];

//	std::cout<<"Element is in "<<(nx+1)*(ny+1)*(nz+1)<<" voxels\n";
	for (int xd=0; xd<=nx; xd++)
		for (int yd=0; yd<=ny; yd++)
			for (int zd=0; zd<=nz; zd++)
			{
				sidx = (idx[0]+xd)*Nz*Ny + (idy[0]+yd)*Nz + idz[0]+zd + boxshift;
			//	std::cout<<"Adding element "<<obj->elem->id()<<" to voxel "<<sidx<<std::endl;
				Voxels[sidx].add(obj);
				obj->VoxelIndices.push_back(sidx);
				//std::cout<<"Element added to voxel. Now adding it to the obj_voxel_map:\n";
				//obj_voxel_map[obj->elem->id()].push_back(std::pair<ShellCacheObj*,int>(obj,sidx));
				//std::cout<<"... done\n";
			}
}

void VoxelList::update(SphereNode* obj)
{
	int sidx;
	int nx=0;
	int ny=0;
	int nz=0;
//	std::cout<<"Updating node "<<obj->node->id()<<std::endl;
	int idx[2], idy[2], idz[2];
	
	obj->VoxelIndices.clear();
	
	idx[0] = obj->Bbox.mmin[0]/rc;
	idx[1] =  obj->Bbox.mmax[0]/rc;
	if (idx[0]<idx[1])
		nx=idx[1]-idx[0];
	idz[0] =  obj->Bbox.mmin[2]/rc;
	idz[1] =  obj->Bbox.mmax[2]/rc;
	if (idz[0]<idz[1])
		nz=idz[1]-idz[0];
	idy[0] =  obj->Bbox.mmin[1]/rc;
	idy[1] =  obj->Bbox.mmax[1]/rc;
	if (idy[0]<idy[1])
		ny=idy[1]-idy[0];
	
//	std::cout<<"Node is in "<<(nx+1)*(ny+1)*(nz+1)<<" voxels\n";
	for (int xd=0; xd<=nx; xd++)
		for (int yd=0; yd<=ny; yd++)
			for (int zd=0; zd<=nz; zd++)
			{
				sidx = (idx[0]+xd)*Nz*Ny + (idy[0]+yd)*Nz + idz[0]+zd + boxshift;
			//	std::cout<<"Adding to ShellNode "<<obj->node->id()<<"'s VoxelIndices the new voxel "<<sidx<<std::endl;
				obj->VoxelIndices.push_back(sidx);
			}	
}


void VoxelList::clear_contacts()
{
	Contacts.clear();	
}

void VoxelList::initialize(RealVectorValue P1, RealVectorValue P2, int Nv, double thickness, int Nelements, int Nnodes, std::vector<ShellCacheObj> *scl)
{
	//Set a pointer to the systems ShellCache
	this->ShellCachePtr = scl;

	Nx = Nv;
	Ny = Nv;
	Nz = Nv;
	double Lx = P2(0)-P1(0);
	double Ly = P2(1)-P1(1);
	double Lz = P2(2)-P1(2);
	rc=Lz/Nv;
	
	boxshift = P1(0)/rc*Nz*Ny + P1(1)/rc*Nz + P1(2)/rc;
	
	Voxels.clear();
	Voxels.resize(Nx*Ny*Nz);
	
	
	if(obj_voxel_map.size()!=Nnodes){
		obj_voxel_map.resize(Nnodes);
		std::cout<<"obj_voxel_map.resize"<<std::endl;
	}
	
	thickness_sq=thickness*thickness;
	this->thickness = thickness;
	n_elements = Nelements;
	n_nodes = Nnodes;
	
	//std::cout<<"VoxelList initialized with thickness_sq="<<thickness_sq<<"\n";
}





void VoxelList::update_tri_tri_contacts() {
	//surf_surf_test.startmeasurementline();
	cfile<<"CONTACT_START\n";
	
	Contacts.clear();
	std::set<int> elem_handled_contacts;
	
	//Iterate over all non-ghost ShellCacheObjects
	for (int i=0; i<ShellCachePtr->size(); i++)
	{
		ShellCacheObj* obj1 = &(*ShellCachePtr)[i];
		//no ghosts
		if (obj1->is_ghost) continue;
		
		//If we are in more than one voxel, keep track of the checks done using the elem_handled_contacts hash.
		int nvoxels = obj1->VoxelIndices.size();
		if (nvoxels > 1) {
			elem_handled_contacts.clear();
		}
		//Now, iterate over all the voxels this element is in:
		for (int j=0; j<nvoxels; j++)
		{

			int voxelid = obj1->VoxelIndices[j];
			
			//Now, find all other non-ghost elements in the same voxel vid:
			for (int obj2id = 0; obj2id < Voxels[voxelid].n_objects(); obj2id++)
			{		
				ShellCacheObj* obj2 = Voxels[voxelid].get_object(obj2id);
				if (obj2->is_ghost)
					continue;
				
				//Continue if obj2's ID is not larger than obj1's (symmetric testing)
				if (obj1->elem->id() >= obj2->elem->id() )
					continue;
				
				//If we overlap in multiple voxels, check if we already handled this contact. If not, add it to the list of handled ones.
				if (nvoxels > 1) {
					if (elem_handled_contacts.count(obj2->elem->id())>0) 
						continue;
					else 
						elem_handled_contacts.insert(obj2->elem->id());
				}
				
				
				//exlude neighbour elements by a check of the normals
				if (ridge_check(obj1,obj2)){
					//std::cout<<"ridge"<<std::endl;
					break;
				}

				//perform a box test
				if (obj1->Bbox.test_overlap(obj2->Bbox)) {
					
					edge_edge_test(obj1, obj2);
					tri_tri_test(obj1, obj2);
					
				}
			}	
		}
	}
  cfile<<"CONTACT_END\n";
}






/*bool VoxelList::ridge_check(const ShellCacheObj* tri1, const ShellCacheObj* tri2) {
	bool in_area = false;
	if((*tri1).neighbour_elements.find(tri2->elem->id())!=(*tri1).neighbour_elements.end()) return true;	// if tri2 is a neighbour of tri1 contact evaluation is immediately suppressed


	if((*tri1).next_neighbour_elements.find(tri2->elem->id())!=(*tri1).next_neighbour_elements.end()){	// if tri2 is a next_neighbour of tri1 contact evaluation is done
		return (tri1->normal*tri2->normal >= -0.0);							// if their normals point against each others, i.e. they're part of a ridge
	}
	
	return false;												// if neither is the case, contact evaluation is done
}
*/

//returns a boolean, true = exlude from contact calc, false = do contact calc
bool VoxelList::ridge_check(const ShellCacheObj* tri1, const ShellCacheObj* tri2) {
	bool in_area = false;
	for (unsigned int i1=0; i1<3; ++i1) {
		Node &n1 = *(tri1->elem->get_node(i1));
		for (unsigned int i2=0; i2<3; ++i2) {
			Node &n2 = *(tri2->elem->get_node(i2));
			
			double d = (n1-n2).size_sq();

			//check if elements are neighbours
			if(d==0.0) return true;
			
			//check if elments distance is smaller than shell_thickness
			if(d < 8.0*thickness_sq) in_area = true;
//			if(tri1->elem->id()==0 && tri2->elem->id()==41) std::cout<<" dist_sq: "<<d<<" 30*thickness_sq: "<<30*thickness_sq<<"boolean d < 30.0*thickness_sq:"<<(d < 30.0*thickness_sq)<<std::endl;
		}
	}

	if(in_area){
		return (tri1->normal*tri2->normal >= -0.0);
	}
	
	return false;
}





//performs edge edge tests on the triangles obj1,obj2
bool VoxelList::edge_edge_test(ShellCacheObj* obj1, ShellCacheObj* obj2){
	RealVectorValue CPT;
	double s1,s2,distsq;
	bool has_contact = false;
	
	//edge edge checks
	for(int i=0;i<3;++i){
		for(int j=0;j<3;++j){
			line_line_cpp(obj1->current_position[i], obj1->current_position[(i+1)%3], obj2->current_position[j], obj2->current_position[(j+1)%3],CPT,s1,s2);
				
			//std::cout<<"seg1:"<<i<<","<<(i+1)%3<<", seg2:"<<j<<","<<(j+1)%3<<std::endl;
				
			distsq=CPT.size_sq();
				
			if(distsq<=thickness_sq){
				TriContact ct;
		
				//std::cout<<"Contact: i,j=("<<i<<","<<j<<"), obj1="<<obj1->elem->id()<<", obj2="<<obj2->elem->id()<<std::endl;
				if(i==0){
					ct.v1 = s1;
					ct.w1 = 0.0;
				}
				else if(i==1){
					ct.v1 = 1.0 - s1;
					ct.w1 = s1;
				}
				else if(i==2){
					ct.v1 = 0.0;
					ct.w1 = 1.0 - s1;
				}
					
				if(j==0){
					ct.v2 = s2;
					ct.w2 = 0.0;
				}
				else if(j==1){
					ct.v2 = 1.0 - s2;
					ct.w2 = s2;
				}
				else if(j==2){
					ct.v2 = 0.0;
					ct.w2 = 1.0 - s2;
				}
					
				ct.obj1 = obj1;
				ct.obj2 = obj2;
				ct.dist = CPT;
				ct.dist2 = distsq;
			
				Contacts.push_back(ct);
					
				has_contact = true;
			}
				
		}
	}
	
	return has_contact;
}



//performs point tri tests on the triangles obj1, obj2
bool VoxelList::tri_tri_test(ShellCacheObj* obj1, ShellCacheObj* obj2){
	RealVectorValue CPT;
	double v,w,distsq;
	bool has_contact = false;
	
	
	const double vpt[]={0.0,1.0,0.0};
	const double wpt[]={0.0,0.0,1.0};
	
	for(int i=0;i<3;++i){
		cpt_point_triangle(obj1->current_position[i], obj2->current_position[0], obj2->current_position[1], obj2->current_position[2],CPT, v, w);
		CPT = obj1->current_position[i] - CPT;
		distsq=CPT.size_sq();
		if(distsq<=thickness_sq){
			TriContact ct;
			
			ct.v1 = vpt[i];
			ct.w1 = wpt[i];
			ct.v2 = v;
			ct.w2 = w;
			
			ct.obj1 = obj1;
			ct.obj2 = obj2;
			ct.dist = CPT;
			ct.dist2 = distsq;
			
			Contacts.push_back(ct);
			
			has_contact = true;
		}
		
		cpt_point_triangle(obj2->current_position[i], obj1->current_position[0], obj1->current_position[1], obj1->current_position[2],CPT, v, w);
		CPT =  CPT - obj2->current_position[i];
		distsq=CPT.size_sq();
		if(distsq<=thickness_sq){
			TriContact ct;
//		std::cout<<"existing contact between: "<<"Obj1_id_it: "<<obj1->elem->id()<<" Obj2_id_it: "<<obj2->elem->id()<<" as :";
//		std::cout<<"dist: "<<CPT.size()<<" v1: "<<v<<" w1: "<<w<<" v2: "<<vpt[i]<<" w2: "<<wpt[i]<<std::endl;
			ct.v1 = v;
			ct.w1 = w;
			ct.v2 = vpt[i];
			ct.w2 = wpt[i];
			
			ct.obj1 = obj1;
			ct.obj2 = obj2;
			ct.dist = CPT;
			ct.dist2 = distsq;
			
			Contacts.push_back(ct);
			
			has_contact = true;
		}
	}
	
	return has_contact;
}




void VoxelList::print_info() 
{
	std::cout<<"VoxelList information: \n";
	std::cout<<"Number of voxels: "<<Voxels.size()<<std::endl;
	std::cout<<"Number of objects in voxels: "<<obj_voxel_map.size()<<std::endl;
	
	std::cout<<"Voxel information for each voxel:\n";
	for (int i=0; i<Voxels.size(); i++)
	{
		std::cout<<"Voxel "<<i<<" has elements\n ";
		for (int j=0; j<Voxels[i].n_objects(); j++)
		{
			ShellCacheObj* obj=Voxels[i].get_object(j);
			std::cout <<"Element ID="<<obj->elem->id()<<"with nodes "<<std::endl;
			std::cout<<"    "<<*(obj->elem->get_node(0))<<std::endl;
			std::cout<<"    "<<*(obj->elem->get_node(1))<<std::endl;
			std::cout<<"    "<<*(obj->elem->get_node(2))<<std::endl;
			std::cout<<"Bounding box ";
			obj->Bbox.print_info();
			std::cout<<std::endl;
		}
	}
}

 
