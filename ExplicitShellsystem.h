/*
 *  shellsystem.h
 *  
 *
 *  Created by Norbert Stoop on 06.04.09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __explicitshellsystem_h__
#define __explicitshellsystem_h__

// Local Includes
#include "dense_matrix.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include "dense_vector.h"
#include "shellsystem.h"
#include "numeric_vector.h"
#include "shellelement.h"
#include "sparse_matrix.h"
#include <fstream>
#include "VoxelList.h"
#include "shell_util.h"
#include "mesh_subdiv_support.h"
// Forward Declarations


class ExplicitShellSystem : public ShellSystem
{

public:

	//*************
	//** CONSTANTS
	//**
	
	//newmark
	static const std::string P_NEWMARK_DT;
	static const std::string P_NEWMARK_MINIMAL_DT;
	static const std::string P_NEWMARK_BETA;
	static const std::string P_NEWMARK_GAMMA;
	static const std::string P_NEWMARK_DAMPING;
	
	//static const std::string P_NEWMARK_ADAPT_START;
	//static const std::string P_NEWMARK_ADAPT_INTERVAL;
	//static const std::string P_NEWMARK_ERROR;
	//static const std::string P_NEWMARK_ERROR_MIN;
	//static const std::string P_NEWMARK_ERROR_MAX;
	
	//boundary conditions
	static const std::string P_BC_TYPE;
	static const std::string P_GHOST_BC_TYPE;
	static const std::string P_GHOST_BC_STIFFNESS;
    
    static const std::string P_CONTAINER_RADIUS;
    static const std::string P_GROWTH_FACTOR_KAPPABAR;
    static const std::string P_GROWTH_FACTOR_L0;

    //contact
    static const std::string P_RESOLVE_CONTACTS;
    static const std::string P_CONTACT_STEPS;
    static const std::string P_CONTACT_R_CRIT;
    static const std::string P_LIMIT_SURFACE_EVALUATION;

	//energies
	static const std::string P_E_KIN;
	static const std::string P_E_TENS;
	static const std::string P_E_BEND;
	static const std::string P_E_CONTACT;
	
	
	
	//*************
	//** METHODS
	//**
	ExplicitShellSystem (EquationSystems& es, const std::string& name, const unsigned int number);
	virtual ~ExplicitShellSystem ();
    
    void find_boundary_nodes();
    
	typedef ExplicitShellSystem sys_type;
	typedef ShellSystem Parent;
	
	//initialization
	void init_data();
	
	//newmark
	void newmark_predict();
	void newmark_correct();
	bool newmark_adapt_timestep();
	void compute_residual();
	
	//force calculation
	void add_uniform_load(double load);
	void add_shear_load(double load);
	void enforce_ghost_bc();
	void compute_lumped_mass();  //Calculates the mass residual
	void assemble() { compute_residual(); }
    void set_initial_volume();

	
	//contact handling
	void collision_detection();
	void handle_contact(TriContact& ct);
	//void handle_surf_contact(SurfContact& ct);
    void mark_excluded_contacts(double radius);
	
	
	
	//andersen thermostat
	//void andersen_thermostat();
	
	
	//caching
	void update_cached_position();
	
	
	
	
	//write data to file
	void write_measurements(std::string file_name);
	
	
	
	//*************
	//** FIELDS
	//**
	
	bool uniform;
	
	bool _mass_assembled;
	double _Econtact_total;
	
	//int center_node_idx;
	bool rebuild_voxel_list;
    std::vector<int> boundary_nodes;
    
protected:
	VoxelList voxellist;
	
	
	void update_current_element_position(ShellCacheObj* selem); //Helper function for the time beeing...
	void update_current_nodal_position(SphereNode* obj);
	void handle_contacts(TriContact& ct);	
	void init_voxellist();
	
	//** newmark
	unsigned int _newmark_adapt_counter;
	unsigned int _newmark_step;
	bool _newmark_reject_timestep;
	

	
};





#endif
