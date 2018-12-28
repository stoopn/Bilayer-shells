/*
 *  ExplicitShellSystem.C
 *  
 *
 *  Created by Norbert Stoop on 06.04.09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <fstream>
#include <cmath>
#include "dof_map.h"
#include "elem.h"
#include "fe_base.h"
#include "face_tri3_sd.h"
#include "fe_interface.h"
#include "equation_systems.h"
#include "shellelement.h"
#include "ExplicitShellsystem.h"
#include "libmesh_logging.h"
#include "mesh.h"
#include "numeric_vector.h"
#include "quadrature.h"
#include "sparse_matrix.h"

//** PARAMETER NAMES

//newmark
const std::string ExplicitShellSystem::P_NEWMARK_DT = "dt";
const std::string ExplicitShellSystem::P_NEWMARK_BETA = "newmark_beta";
const std::string ExplicitShellSystem::P_NEWMARK_GAMMA = "newmark_gamma";
const std::string ExplicitShellSystem::P_NEWMARK_DAMPING = "newmark_damping";

//boundary conditions
const std::string ExplicitShellSystem::P_BC_TYPE = "bc_type";
const std::string ExplicitShellSystem::P_GHOST_BC_TYPE = "ghost_bc_type";
const std::string ExplicitShellSystem::P_GHOST_BC_STIFFNESS = "ghost_bc_stiffness";

const std::string ExplicitShellSystem::P_GROWTH_FACTOR_KAPPABAR = "growth_factor_kappabar";
const std::string ExplicitShellSystem::P_GROWTH_FACTOR_L0 = "growth_factor_L0";

//spherical container
const std::string ExplicitShellSystem::P_CONTAINER_RADIUS = "container_radius";

//contact
// const std::string ExplicitShellSystem::P_RESOLVE_CONTACTS = "resolve_contacts";
// const std::string ExplicitShellSystem::P_CONTACT_STEPS = "contact_steps";
// const std::string ExplicitShellSystem::P_CONTACT_R_CRIT = "contact_r_crit";
// const std::string ExplicitShellSystem::P_LIMIT_SURFACE_EVALUATION = "limit_surface_evaluation";

//energies
const std::string ExplicitShellSystem::P_E_KIN = "e_kin";
const std::string ExplicitShellSystem::P_E_TENS = "e_tens";
const std::string ExplicitShellSystem::P_E_BEND = "e_bend";
const std::string ExplicitShellSystem::P_E_CONTACT = "e_contact";


ExplicitShellSystem::ExplicitShellSystem (EquationSystems& es, const std::string& name, const unsigned int number) 
	: Parent(es, name, number),
	_mass_assembled(false),
    uniform(false)	
{
	
	//add vectors and initialize them to zero
	this->add_vector("velocity",false);
	this->add_vector("acceleration",false);
	this->add_vector("lumped_mass",false);
	
	this->add_vector("newmark_u",false);
	this->add_vector("newmark_v",false);
	this->add_vector("newmark_a",false);
	
	this->add_vector("newmark_err",false);
		
	//** Initialize
    _use_cache = true;
	//newmark
	_newmark_step = 0;
	_newmark_adapt_counter = 0;
	
	//energies
	_Etens_total = 0.;
	_Ebend_total = 0.;
	_Eext_total = 0.;
	_Econtact_total = 0.;
	
	//time
	time = 0.0;
	deltat = es.parameters.get<double>(P_NEWMARK_DT);
}


ExplicitShellSystem::~ExplicitShellSystem() {
	
}


void ExplicitShellSystem::init_data() {
	Parent::init_data();
    //std::cout<<"here0\n";

    //std::cout<<"here011\n";

    //std::cout<<"here012\n";

	const MeshBase& mesh = this->get_mesh();
    //std::cout<<"here01\n";

    //std::cout<<"here02\n";

	const double &container_radius = param_->get<double>(P_CONTAINER_RADIUS);
    //std::cout<<"here03\n";

	//get voxelgrid size
	double radius = 1.5*container_radius;
	RealVectorValue P1(radius,-radius,-radius);
	RealVectorValue P2(-radius,radius,radius);
    //std::cout<<"here1\n";
	voxellist.initialize(P1,P2, 20, material.thickness, mesh.n_elem(), mesh.n_nodes(), &ShellCache);

    shellelement->set_growth(1.,1.); // This means no growth
    
    //std::cout<<"here2\n";
	//Iterate over the real system mesh once to initialize the obj_voxel_map and ShellCacheObj with their mesh counterparts:
	
	MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
	
	for ( ; el != end_el; ++el)
	{
		Elem* elem = *el;
		const Tri3SD* sdelem = static_cast<const Tri3SD*>(elem);
		if (sdelem->is_ghost()) {
			for (int nn=0; nn<elem->n_nodes(); ++nn)
			{
				Node* node = elem->get_node(nn);
				SphereNode* sp= voxellist.get_node(node->id());
				sp->node = node;
				ShellCache[elem->id()].is_ghost=true;
				ShellCache[elem->id()].elem = elem;
			}
		} else {
            
			for (int nn=0; nn<elem->n_nodes(); ++nn)
			{
				Node* node = elem->get_node(nn);
				SphereNode* sp= voxellist.get_node(node->id());
				sp->node = node;
				sp->is_ghost=false;
				ShellCache[elem->id()].is_ghost=false;
			}
        }
    }
    find_boundary_nodes();
}

void ExplicitShellSystem::find_boundary_nodes() {
// Finds the ids of nodes that are at the boundary of the mesh
    const MeshBase& mesh = this->get_mesh();

    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    
    this->boundary_nodes.resize(0);
    //loop over all ghost elements
    for ( ; el != end_el; ++el)
    {
        Elem* elem = *el;
        const Tri3SD* sdelem = static_cast<const Tri3SD*>(elem);
        
        if (!sdelem->is_ghost())
            continue;
        
        //find the side on which a real element sits:
        for (unsigned int s=0; s<elem->n_sides(); s++)
        {
            if (elem->neighbor(s) == NULL)
                continue;
            
            Tri3SD* nbelem = static_cast<Tri3SD*>(elem->neighbor(s));
            
            if (nbelem->is_ghost())
                continue;
            
            int id, dofidx;
            
            Node* n1 = sdelem->get_node(s);
            Node* n2 = sdelem->get_node((s+1)%3);
            this->boundary_nodes.push_back(n1->id());
            this->boundary_nodes.push_back(n2->id());
            
        }
    }
    
   std::sort(this->boundary_nodes.begin(), this->boundary_nodes.end());
    std::vector<int>::iterator last = std::unique(this->boundary_nodes.begin(), this->boundary_nodes.end());
    this->boundary_nodes.erase(last, this->boundary_nodes.end());
}



void ExplicitShellSystem::set_initial_volume() {
    const MeshBase& mesh = this->get_mesh();
   // std::cout<<"here3\n";
    this->volume0=0.;
    
    shellelement->set_growth(0.0, 1.0); // kappabar, L0
   // shellelement->set_uniform(uniform);
    //Determine initial volume (volume0):
    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    for ( ; el != end_el; ++el)
    {
        Elem* elem = *el;
        const Tri3SD* sdelem = static_cast<const Tri3SD*>(elem);
        if (sdelem->is_ghost()) {
            continue;
        }
        
        shellelement->reinit(elem);
        //Get volume...
        this->volume0 += shellelement->get_volume(true); //True means we need to setup the metric etc.
        
    }
    volume=volume0;
    std::cout<<"Volume at start: "<<this->volume0<<std::endl;
   
}

void ExplicitShellSystem::init_voxellist() {
	const MeshBase& mesh = this->get_mesh();

	//get voxelbox size
	double radius = 1.5*param_->get<double>(P_CONTAINER_RADIUS);
	const double &thickness = param_->get<double>(P_THICKNESS);
	
	
	RealVectorValue P1(radius,-radius,-radius);
	RealVectorValue P2(-radius,radius,radius);
	
	//triangle edge size estimation
	double n_obj = ShellCache.size()-1;
	double avg_edge_size = 0.0;
	unsigned int n_edge = 0;
	unsigned int n_elem = std::sqrt(mesh.n_elem());

	for(unsigned int i=0;i<n_elem;++i){
		unsigned int j = drand48()*n_obj;
		ShellCacheObj &scobj = ShellCache[j];
		
		avg_edge_size += (scobj.current_position[0] - scobj.current_position[1]).size();
		avg_edge_size += (scobj.current_position[1] - scobj.current_position[2]).size();
		avg_edge_size += (scobj.current_position[0] - scobj.current_position[2]).size();
		
		n_edge+=3;
	}
	
	avg_edge_size = 2.0*(avg_edge_size+thickness)/(double)n_edge;
	unsigned int n_voxels = 1.5*radius/avg_edge_size;
	
	voxellist.initialize(P1,P2, n_voxels, material.thickness, mesh.n_elem(), mesh.n_nodes(), &ShellCache);
	
	//First, iterate over all nodes. Instead of iterating over mesh nodes, we directly iterate over voxellist ShellNodes, which is equivalent.
	for (int nn=0; nn<voxellist.get_n_nodes(); nn++) {
		voxellist.update(voxellist.get_node(nn));		
	}
	
	//Now, iterate through all local elements and leave out the ghosts
	//int voxels_per_element=0;
	//int n_elem_ = 0;
	
	MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
	for ( ; el != end_el; ++el)
	{
		Elem* elem = *el;
		const Tri3SD* sdelem = static_cast<const Tri3SD*>(elem);
		
		if (sdelem->is_ghost())
			continue;
		
		voxellist.add(&ShellCache[elem->id()]);
		//voxels_per_element += ShellCache[elem->id()].VoxelIndices.size();
		//++n_elem_;
	}
	//std::cout<<"voxels per element = "<<(double)voxels_per_element/(double)n_elem_<<std::endl;
}


void ExplicitShellSystem::update_cached_position() {
	const MeshBase& mesh = this->get_mesh();
	
	for (int nn=0; nn<voxellist.get_n_nodes(); nn++)
	{
		//Update the current nodal position and the node_to_voxel mapping for this node
		update_current_nodal_position(voxellist.get_node(nn));
	}
	
	//Now, iterate through all local elements and leave out the ghosts
	MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
	
	for ( ; el != end_el; ++el)
	{
		Elem* elem = *el;
		const Tri3SD* sdelem = static_cast<const Tri3SD*>(elem);
				
		update_current_element_position(&ShellCache[elem->id()]);
		//	std::cout<<"Rebuilding voxel list for this element\n";
	}
}


void ExplicitShellSystem::mark_excluded_contacts(double radius) {
    
    const MeshBase& mesh = this->get_mesh();
    
    
    //Now, iterate through all local elements and leave out the ghosts
    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    RealVectorValue cv;
    cv.zero();
    /*
     if (center_node_idx<0) {
     std::cout<<"WARNING: No valid center node idx, assuming (0,0,0)\n";
     cv.zero();
     }  else {
     const Node* cnode = mesh.node_ptr(center_node_idx);
     cv = *cnode;
     }
     */
    std::cout<<"Position of center: "<<cv<<std::endl;
    std::cout<<"out1\n";
    
    for ( ; el != end_el; ++el)
    {
        Elem* elem = *el;
        const Tri3SD* sdelem = static_cast<const Tri3SD*>(elem);
        
        if (sdelem->is_ghost())
            continue;
        
        int totals=0;
        for (int nn=0; nn<elem->n_nodes(); ++nn)
        {
            Node* node = elem->get_node(nn);
            RealVectorValue dv = (*node)-cv;
            std::cout<<dv<<std::endl;
            if (dv.size_sq()<radius*radius) {
                SphereNode* sp= voxellist.get_node(node->id());
                sp->exclude=true;
                totals++;
            }
        }
        if (totals==3)
            ShellCache[elem->id()].exclude=true;
        
    }
    
}


//**********************
//** CACHE
//**
void ExplicitShellSystem::update_current_nodal_position(SphereNode* obj) {
  
	const Node& node = *(obj->node);
	NumericVector<Number>& vel = this->get_vector("velocity");
	
	obj->current_position(0) = node(0) + (*solution)(node.dof_number(0,u_var,0));
	obj->current_position(1) = node(1) + (*solution)(node.dof_number(0,v_var,0));
	obj->current_position(2) = node(2) + (*solution)(node.dof_number(0,w_var,0));
	obj->current_velocity(0) = vel(node.dof_number(0,u_var,0));
	obj->current_velocity(1) = vel(node.dof_number(0,v_var,0));
	obj->current_velocity(2) = vel(node.dof_number(0,w_var,0));

	obj->Bbox.update_from_sphere(obj->current_position, 0.5*material.thickness);
	obj->normalct=0;
	obj->normal.zero();
}

void ExplicitShellSystem::update_current_element_position(ShellCacheObj* obj) {
	const Elem* elem = obj->elem;
	SphereNode* sp;
	RealVectorValue v1 = obj->current_position[1]-obj->current_position[0];
	RealVectorValue v2 = obj->current_position[2]-obj->current_position[1];
	obj->normal=v1.cross(v2);
	//	std::cout<<"Elemnorm: "<<obj->normal<<std::endl;
	for (int nn=0; nn<elem->n_nodes(); nn++)
	{
	  sp = voxellist.get_node(elem->get_node(nn)->id());
	  obj->current_position[nn] = sp->current_position;
	  obj->current_velocity[nn] = sp->current_velocity;
	  //Now, we should also update the average facet normal at each triangle node (for normal contact check)
	  sp->normal+=obj->normal;
	  sp->normalct++;
	}
	obj->Bbox.update_from_triangle(obj->current_position[0], obj->current_position[1], obj->current_position[2], material.thickness);
}



//**********************
//** CONTACT HANDLING
//**

void ExplicitShellSystem::collision_detection() {
    //get parameters
    
    double &dt			= param_->set<double>(P_NEWMARK_DT);
    //const bool &resolve_contacts 	= param_->get<bool>(P_RESOLVE_CONTACTS);
    //const bool &limit_surf_eval 	= param_->get<bool>(P_LIMIT_SURFACE_EVALUATION);
    
    if (measure_now)
        _Econtact_total = 0.;
    
    init_voxellist();
    
    if (rebuild_voxel_list) {
        rebuild_voxel_list = false;
        init_voxellist();
    }
    voxellist.update_tri_tri_contacts();
    
    for (std::vector<TriContact>::iterator it=voxellist.Contacts.begin(); it!=voxellist.Contacts.end(); ++it) {
        handle_contact(*it);
        
    }
}


//handle point-surface contact
void ExplicitShellSystem::handle_contact(TriContact& ct) {
    
    Node* fnode;
    double mass;
    double E_ct;
    
    const Elem& el1 = *(ct.obj1->elem);
    const Elem& el2 = *(ct.obj2->elem);
    
    //distance
    double d = sqrt(ct.dist2);
    
    //force
    RealVectorValue force = material.contact_stiffness*(material.thickness-d)*(material.thickness-d)/(ct.dist2*d)*ct.dist;
    //energy
    if (measure_now){
        E_ct= 0.5*material.contact_stiffness*(material.thickness-d)*(material.thickness-d)/(ct.dist2);
        _Econtact_total += E_ct;
    }
    
    double u1 = 1.0 - ct.v1 - ct.w1;
    fnode = el1.get_node(0);
    mass =  -u1/(this->get_vector("lumped_mass")(fnode->dof_number(0,u_var,0)));
    if(measure_now) msys->add_nodal_value(fnode, 2, -mass*E_ct);
    rhs->add(fnode->dof_number(0,u_var,0), mass*force(0));
    rhs->add(fnode->dof_number(0,v_var,0), mass*force(1));
    rhs->add(fnode->dof_number(0,w_var,0), mass*force(2));
    
    fnode = el1.get_node(1);
    mass =  -ct.v1/(this->get_vector("lumped_mass")(fnode->dof_number(0,u_var,0)));
    if(measure_now) msys->add_nodal_value(fnode, 2, -mass*E_ct);
    rhs->add(fnode->dof_number(0,u_var,0), mass*force(0));
    rhs->add(fnode->dof_number(0,v_var,0), mass*force(1));
    rhs->add(fnode->dof_number(0,w_var,0), mass*force(2));
    
    fnode = el1.get_node(2);
    mass =  -ct.w1/(this->get_vector("lumped_mass")(fnode->dof_number(0,u_var,0)));
    if(measure_now) msys->add_nodal_value(fnode, 2, -mass*E_ct);
    rhs->add(fnode->dof_number(0,u_var,0), mass*force(0));
    rhs->add(fnode->dof_number(0,v_var,0), mass*force(1));
    rhs->add(fnode->dof_number(0,w_var,0), mass*force(2));
    
    
    double u2 = 1.0 - ct.v2 - ct.w2;
    fnode = el2.get_node(0);
    mass =  u2/(this->get_vector("lumped_mass")(fnode->dof_number(0,u_var,0)));
    if(measure_now) msys->add_nodal_value(fnode, 2, mass*E_ct);
    rhs->add(fnode->dof_number(0,u_var,0), mass*force(0));
    rhs->add(fnode->dof_number(0,v_var,0), mass*force(1));
    rhs->add(fnode->dof_number(0,w_var,0), mass*force(2));
    
    fnode = el2.get_node(1);
    mass =  ct.v2/(this->get_vector("lumped_mass")(fnode->dof_number(0,u_var,0)));
    if(measure_now) msys->add_nodal_value(fnode, 2, mass*E_ct);
    rhs->add(fnode->dof_number(0,u_var,0), mass*force(0));
    rhs->add(fnode->dof_number(0,v_var,0), mass*force(1));
    rhs->add(fnode->dof_number(0,w_var,0), mass*force(2));
    
    fnode = el2.get_node(2);
    mass =  ct.w2/(this->get_vector("lumped_mass")(fnode->dof_number(0,u_var,0)));
    if(measure_now) msys->add_nodal_value(fnode, 2, mass*E_ct);
    rhs->add(fnode->dof_number(0,u_var,0), mass*force(0));
    rhs->add(fnode->dof_number(0,v_var,0), mass*force(1));
    rhs->add(fnode->dof_number(0,w_var,0), mass*force(2));
    
}



// Build the residual and jacobian contributions on each element
void ExplicitShellSystem::compute_residual() {
	START_LOG("compute_residual()", "ExplicitShellSystem");
		
	const MeshBase& mesh = this->get_mesh();
    double nvolume=0.;
    
	if (measure_now) {
		_Ebend_total=0.;
		_Etens_total=0.;
		_Eext_total=0.;
		shellelement->measure_energy = true;
		msys->reset_measurements();
	} else
		shellelement->measure_energy = false;

	shellelement->set_growth(param_->get<double>(P_GROWTH_FACTOR_KAPPABAR), param_->get<double>(P_GROWTH_FACTOR_L0));
    shellelement->set_uniform(true);
    shellelement->set_volume_constrained(true);

	
	//loop over all non-ghost elements
	MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
	for ( ; el != end_el; ++el) {
		Elem* elem = *el;
		const Tri3SD* sdelem = static_cast<const Tri3SD*>(elem);
		
		if (sdelem->is_ghost())		
			continue;

		//Get the residual
		shellelement->reinit(elem);
		shellelement->assembly(true, false); 
		
		//Adds all additional loads as specified within the cache object (for instance, contact forces)
		shellelement->add_external_loads(get_cache_obj(elem)); 

		this->rhs->add_vector (shellelement->residual, shellelement->dof_indices);
		if (measure_now) {
		 	_Ebend_total += shellelement->Ebend;
			_Etens_total += shellelement->Etens;
			_Eext_total += shellelement->Eext;
            
			msys->add_element_value(elem, 0, shellelement->Etens/shellelement->ElRefArea);
			msys->add_element_value(elem, 1, shellelement->Ebend/shellelement->ElRefArea);
            msys->add_element_value(elem, 3, shellelement->kGauss/shellelement->ElDefArea);
            msys->add_element_value(elem, 4, shellelement->kMean/shellelement->ElDefArea);
            msys->add_element_value(elem, 5, shellelement->ElRefArea);
            msys->add_element_value(elem, 6, shellelement->ElDefArea);
		}
        
        nvolume += shellelement->get_volume(false); // false means the volume was already calculated
	}
    this->volume = nvolume;
  //  std::cout<<"Current volume: "<<this->volume<<", target volume: "<<this->volume0<<std::endl;
	STOP_LOG("compute_residual()", "ExplicitShellSystem");
}

// compute the mass of the elements
void ExplicitShellSystem::compute_lumped_mass() {
	libmesh_assert(_mass_assembled == false);

	const MeshBase& mesh = this->get_mesh();
	//this->update();

	//set mass vector to zero
	NumericVector<Number>& mass = this->get_vector("lumped_mass");
	mass.zero();
	
	MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
	
	//FIXME: Lumped mass is calculated from zero growth state
  //  std::cout<<"here1\n";
	shellelement->set_growth(1.0,1.0);
    //std::cout<<"here2\n";
	//loop over all non-ghost elements
	for ( ; el != end_el; ++el) {
		Elem* elem = *el;
	
		const Tri3SD* sdelem = static_cast<const Tri3SD*>(elem);
		
		if (sdelem->is_ghost())		
			continue;
      //  std::cout<<"here3\n";
		shellelement->reinit(elem);	
//	std::cout<<"here4\n";
		//Get the element mass
		shellelement->element_mass(); 
		mass.add_vector(shellelement->mass, shellelement->dof_indices);
	}

	mass.close();

	//invert all entries in mass vector
	for (int ii=0; ii<mass.size(); ++ii) {
		mass.set(ii, 1./mass(ii));
	}
	_mass_assembled = true;
}


//**********************
//** NEWMARK INTEGRATOR
//**
//u_t+1 = u_t + v_t*dt + ((1-2*beta)*a_t + 2*beta*a_t+1)*dt^2/2
//v_t+1 = v_t + ( (1-gamma)*a_t + gamma*a_t+1 )*dt

//newmark predict
void ExplicitShellSystem::newmark_predict() {
	//get parameters
/*	const double &dt = param_->get<double>(P_NEWMARK_DT);
	const double &beta = param_->get<double>(P_NEWMARK_BETA);
	const double &gamma = param_->get<double>(P_NEWMARK_GAMMA);
	
	NumericVector<Number>&  vel    = this->get_vector("velocity");
	NumericVector<Number>&  acc    = this->get_vector("acceleration");
	
	NumericVector<Number>&  newmark_u    = this->get_vector("newmark_u");
	NumericVector<Number>&  newmark_v    = this->get_vector("newmark_v");
	NumericVector<Number>&  newmark_a = this->get_vector("newmark_a");
	
	
	//update material
	material.thickness = param_->get<double>(P_THICKNESS);
	material.young = param_->get<double>(P_YOUNG_MODULUS);
	material.poisson = param_->get<double>(P_POISSON_RATIO);
	material.density = param_->get<double>(P_DENSITY);
	material.contact_stiffness = 1;
	material.init();
*/
    solution->close();
/*
 //   vel.close();
    
	//buffer u,v
//	newmark_u = *solution;
//	newmark_v = vel;
//	newmark_a = acc;

	//u_t+1 = u_t + v_t*dt + ..
	solution->add(dt, acc);  //Overdamped dynamics
	// .. + (1-2*beta)*a_t*dt^2/2 + ..
	//solution->add((0.5-beta)*dt*dt, acc);
	
	//v_t+1 = (1-gamma)*a_t*dt + ..
	//vel.add((1.0 -gamma)*dt, acc);
	*/
	rhs->zero();

	this->update();
/*
	_newmark_reject_timestep = false;*/
}

//newmark correct 
void ExplicitShellSystem::newmark_correct() {
	//get parameters
	const double &dt = param_->get<double>(P_NEWMARK_DT);
	const double &beta = param_->get<double>(P_NEWMARK_BETA);
	const double &gamma = param_->get<double>(P_NEWMARK_GAMMA);
	const double &damping = param_->get<double>(P_NEWMARK_DAMPING);
	//double &newmark_error	= param_->set<double>(P_NEWMARK_ERROR);
	
	/*
	NumericVector<Number>&  vel			= this->get_vector("velocity");
	NumericVector<Number>&  acc			= this->get_vector("acceleration");
	NumericVector<Number>&  mass		= this->get_vector("lumped_mass");
	
	NumericVector<Number>&  newmark_u = this->get_vector("newmark_u");
	NumericVector<Number>&  newmark_v = this->get_vector("newmark_v");
	NumericVector<Number>&  newmark_a = this->get_vector("newmark_a");
	NumericVector<Number>&  newmark_err = this->get_vector("newmark_err");
	*/
    rhs->close();
   // acc.close();
    solution->close();
   // vel.close();
    
	rhs->scale(-1.);  //Residual = -force !
	//acc.pointwise_mult(*rhs,mass);  //mass here is the inverse mass!!
	//acc.add(-damping, vel);
	
	// u_t+1 = .. + 2*beta*a_t+1*dt^2/2
	//if(beta!=0.0)
    solution->add(dt, *rhs);
	
	// .. + gamma*a_t+1*dt
	//vel.add(gamma*dt, acc);
	
}

//**********************
//** GHOST BC
//**
//apply constraints for the ghost nodes.
void ExplicitShellSystem::enforce_ghost_bc() {
	
	//get parameters
	const int &bc_type = param_->get<int>(P_BC_TYPE);
	const int &ghost_bc_type = param_->get<int>(P_GHOST_BC_TYPE);
	const double &ghost_bc_stiffness = param_->get<double>(P_GHOST_BC_STIFFNESS);
	
	
	//If we have free propagation of ghosts and free BCs, we do not need to apply any BCs here. 
	// For clamped BCs, it does not matter...
	
	if (bc_type==0 && ghost_bc_type == 0)
	  return;
	
	
	//boundary_conditions==0 means free boundary conditions, ie. mirror the node solution.
	const MeshBase& mesh = this->get_mesh();
	NumericVector<Number>&  vel    = this->get_vector("velocity");

	MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
	

	//loop over all ghost elements
	for ( ; el != end_el; ++el)
	{
		Elem* elem = *el;
		const Tri3SD* sdelem = static_cast<const Tri3SD*>(elem);
		
		if (!sdelem->is_ghost())		
			continue;
		
		//find the side on which a real element sits:
		for (unsigned int s=0; s<elem->n_sides(); s++)
		{
			if (elem->neighbor(s) == NULL)
				continue;
			
			Tri3SD* nbelem = static_cast<Tri3SD*>(elem->neighbor(s));
						
			if (nbelem->is_ghost())
				continue;
			
			int id, dofidx;
			RealVectorValue mirrored, mirrored_vel;
			ShellCacheObj* obj = &ShellCache[nbelem->id()];
			
			//** free bc
			if (bc_type==0)
			{
				
				// nbelem is now a real element. Find out on which side of sdelem it sits and 
				// calcualte the mirrored positions and velocities
				//
				//            u_gh =  u_1 + u_2 - u_3
				//             *
				//            / \
				//          /     \
				//     u_1 *-------* u_2
				//          \     /
				//            \ /
				//             *
				//            u_3
				
				
				if (sdelem->get_node(s) == nbelem->get_node(0)) {
					mirrored = (obj->current_position[0]-*(nbelem->get_node(0))) 
					+ (obj->current_position[2]-*(nbelem->get_node(2))) 
					- (obj->current_position[1]-*(nbelem->get_node(1)));
					
					mirrored_vel = obj->current_velocity[0] 		     
				    +  obj->current_velocity[2]                         
				    -  obj->current_velocity[1];
					
				} else if (sdelem->get_node(s) == nbelem->get_node(1)) {
					mirrored = (obj->current_position[0]-*(nbelem->get_node(0))) 
					+ (obj->current_position[1]-*(nbelem->get_node(1))) 
					- (obj->current_position[2]-*(nbelem->get_node(2)));
					
					mirrored_vel = obj->current_velocity[0]
					+  obj->current_velocity[1]
					-  obj->current_velocity[2];
					
				} else if (sdelem->get_node(s) == nbelem->get_node(2)) {
					mirrored = (obj->current_position[1]-*(nbelem->get_node(1))) 
					+ (obj->current_position[2]-*(nbelem->get_node(2))) 
					- (obj->current_position[0]-*(nbelem->get_node(0)));
					
					mirrored_vel = obj->current_velocity[1]
					+  obj->current_velocity[2]
					-  obj->current_velocity[0];
				}
			
				//** hard mirrored bc
				if(ghost_bc_type==1){
                    id = sdelem->get_node(MeshTools::Subdiv::prev[s])->dof_number(0,u_var,0);
					solution->set(id, mirrored(0));
					vel.set(id,0.);	rhs->set(id,0.);
					
					id = sdelem->get_node(MeshTools::Subdiv::prev[s])->dof_number(0,v_var,0);
					solution->set(id, mirrored(1));
					vel.set(id,0.); rhs->set(id,0.);
					
					id = sdelem->get_node(MeshTools::Subdiv::prev[s])->dof_number(0,w_var,0);
					solution->set(id, mirrored(2));
					vel.set(id,0.); rhs->set(id,0.);
					
				}
				
				//** soft mirrored bc
				else {
					id = sdelem->get_node(MeshTools::Subdiv::prev[s])->dof_number(0,u_var,0);
					rhs->add(id, -ghost_bc_stiffness*(mirrored(0)-(*solution)(id)));
					
					id = sdelem->get_node(MeshTools::Subdiv::prev[s])->dof_number(0,v_var,0);
					rhs->add(id, -ghost_bc_stiffness*(mirrored(1)-(*solution)(id)));
					
					id = sdelem->get_node(MeshTools::Subdiv::prev[s])->dof_number(0,w_var,0);
					rhs->add(id, -ghost_bc_stiffness*(mirrored(2)-(*solution)(id)));
				}
        
				
			//** clamped bc
			} else if (bc_type==1) {
				//Set all nodes of the ghost and of the neighbor to zero
				for (int nni=0; nni<3; nni++)
				{
					//set u,v,w of the nodes of sdelem to zero
					id = sdelem->get_node(nni)->dof_number(0,u_var,0);
					solution->set(id, 0.0); vel.set(id,0.0);
					id = sdelem->get_node(nni)->dof_number(0,v_var,0);
					solution->set(id, 0.0); vel.set(id,0.0);
					id = sdelem->get_node(nni)->dof_number(0,w_var,0);
					solution->set(id, 0.0); vel.set(id,0.0);
					
					//set u,v,w of the nodes of ndelem to 0.0
					id = nbelem->get_node(nni)->dof_number(0,u_var,0);
					solution->set(id, 0.0); vel.set(id,0.0);
					id = nbelem->get_node(nni)->dof_number(0,v_var,0);
					solution->set(id, 0.0); vel.set(id,0.0);
					id = nbelem->get_node(nni)->dof_number(0,w_var,0);
					solution->set(id, 0.0); vel.set(id,0.0);
				}
				
			}

		}
	}
	
	this->update();
}

//add uniform load

void ExplicitShellSystem::add_uniform_load(double load) {
	const MeshBase& mesh = this->get_mesh();
	
	MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
	
	//loop over all non-ghost elements
	for ( ; el != end_el; ++el) {
		Elem* elem = *el;
		const Tri3SD* sdelem = static_cast<const Tri3SD*>(elem);
		
		if (sdelem->is_ghost())		
			continue;
		
		ShellCacheObj* obj = &ShellCache[elem->id()];
		obj->load(2)=load;	
	}
} 



//calc Ekin, T and writes data to file.
void ExplicitShellSystem::write_measurements(std::string file_name) {
	//const int &collision_rate = param_->get<int>(P_ANDERSEN_COLLISION_RATE);
	
	//open file
	std::ofstream mfile;
	mfile.open(file_name.c_str(), std::ios::app|std::ios::out);
    if (!mfile.is_open()) {
		std::cout<<"At time "<<this->time<<": Measurement file could not be opened! Continuing...\n";
		return;	
    }
	
	NumericVector<Number>&  vel    = this->get_vector("velocity");
	NumericVector<Number>&  mass   = this->get_vector("lumped_mass");
	double Ekin_total = 0.0;

    vel.close();
    mass.close();

	int active=0;
	//loop over all non-ghost nodes
	for (int nn=0; nn<voxellist.get_n_nodes(); nn++) {
		SphereNode* sn=voxellist.get_node(nn);
		if (sn->is_ghost)
			continue;
		
		Node* node = sn->node;
		
		// Ekin = 0.5*m*|v|^2
		double ms = 1./mass(node->dof_number(0,u_var,0));
        double vu=vel(node->dof_number(0,u_var,0));
        double vv=vel(node->dof_number(0,v_var,0));
        double vw=vel(node->dof_number(0,w_var,0));
		Ekin_total += ms*(vu*vu+vv*vv+vw*vw);  
		
		active++;
	}
	Ekin_total*=0.5;
	
	param_->set<double>(P_E_KIN) = Ekin_total;
	param_->set<double>(P_E_TENS) = _Etens_total;
	param_->set<double>(P_E_BEND) = _Ebend_total;
	param_->set<double>(P_E_CONTACT) = _Econtact_total;
    double Etot =Ekin_total+this->_Etens_total+this->_Ebend_total+ this->_Eext_total;
	
    // Calculate the asymmetry of the boundary nodes:
    const MeshBase& mesh = this->get_mesh();
    std::vector<Point> pos;
    pos.resize(0);
    for (std::vector<int>::iterator it = this->boundary_nodes.begin(); it != this->boundary_nodes.end(); ++it)
    {
        const Node* node = mesh.node_ptr(*it);
        double u=(*solution)(node->dof_number(0,u_var,0));
        double v=(*solution)(node->dof_number(0,v_var,0));
        double w=(*solution)(node->dof_number(0,w_var,0));
        Point disp(u,v,w);
        pos.push_back(*node + disp);
    }
    // Calculate center of mass:
    Point COM(0,0,0);
    for (std::vector<Point>::iterator it = pos.begin(); it != pos.end(); ++it)
    {
        COM.add(*it);
    }
    COM/=pos.size();
    // Now obtain the 2d inertia tensor (x and y parts):
    double Ixx,Iyy,Ixy;
    Ixx=Iyy=Ixy=0;
    for (std::vector<Point>::iterator it = pos.begin(); it != pos.end(); ++it)
    {
        Point diff = *it - COM;
        double x=diff(0);
        double y=diff(1);
        //double z=diff(2);
        Ixx += y*y;
        Iyy += x*x;
        Ixy += -x*y;
    }
    Ixx/=pos.size();
    Iyy/=pos.size();
    Ixy/=pos.size();

    double l1,l2; //Eigenvalues
    l1 = ((Ixx+Iyy) + sqrt((Ixx-Iyy)*(Ixx-Iyy) + 4*Ixy*Ixy))/2.0;
    l2 = ((Ixx+Iyy) - sqrt((Ixx-Iyy)*(Ixx-Iyy) + 4*Ixy*Ixy))/2.0;
    
    double kappabar = param_->get<double>(P_GROWTH_FACTOR_KAPPABAR);
    double L0 = param_->get<double>(P_GROWTH_FACTOR_L0);
	//write file
	mfile<< std::scientific<<this->time<<"\t"<<Etot<<"\t"<<Ekin_total<<"\t"<<this->_Etens_total<<"\t"<<this->_Ebend_total<<"\t"
    <<this->_Eext_total<<"\t"
    << kappabar<<"\t"
    << L0<<"\t"
    << 1./sqrt(l1) <<"\t"
    << 1./sqrt(l2) <<"\t"
    << volume0 << "\t"
    << volume << "\t"
    << 0.5*mu_vol*(volume-volume0)*(volume-volume0)
    <<"\n";
	mfile.close();

    std::cout<<"Current kinetic energy: "<<Ekin_total<<"; total elastic energy: "<<Etot<<"\n";
	measure_now = false;
}

