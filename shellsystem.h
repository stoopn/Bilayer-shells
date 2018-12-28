/*
 *  shellsystem.h
 *  
 *
 *  Created by Norbert Stoop on 06.04.09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef __shellsystem_h__
#define __shellsystem_h__

// Local Includes
#include "libmesh_common.h"
#include "dense_matrix.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include "dense_vector.h"
#include "implicit_system.h"
#include "transient_system.h"
#include "numeric_vector.h"
#include "shellelement.h"
#include "sparse_matrix.h"
#include <fstream>
#include "mesh_tools.h"
#include "enum_xdr_mode.h"
#include "xdr_cxx.h"
#include "parallel_mesh.h"
#include "serial_mesh.h"
#include "equation_systems.h"
#include "shellmeasurementsystem.h"
#include "fe_base.h"

// Forward Declarations
//class FEBase;
//template <typename T> class NumericVector;

class ShellSystem : public ImplicitSystem
{
public:
	
	//*************
	//** CONSTANTS
	//**
	
	//integration
	static const std::string P_QUADRATURE_EXTRA_ORDER;
	
	static const std::string P_POISSON_RATIO;
	static const std::string P_YOUNG_MODULUS;
	static const std::string P_DENSITY;
	static const std::string P_THICKNESS;
	
	static const std::string P_CONTACT_STIFFNESS;
	
	
	
	
	
	ShellSystem (EquationSystems& es, const std::string& name, const unsigned int number);
	virtual ~ShellSystem ();

	typedef ShellSystem sys_type;
	typedef ImplicitSystem Parent;
  
	/**
	* Clear all the data structures associated with
	* the system. 
	*/
	void clear ();

	/**
	* Reinitializes the member data fields associated with
	* the system, so that, e.g., \p assemble() may be used.
	*/
	void reinit ();
	
	
	void init_data ();
		
	AutoPtr<ShellElement> shellelement;
	void init_shell_element();
	void update_time(double ntime);

	ShellCacheObj* get_cache_obj(const Elem *elem) const;
	bool has_cache_obj(const Elem *elem) const;
	void add_cache_obj(const Elem *elem, ShellCacheObj* obj);
	
	
	ShellMaterial material;
	unsigned int u_var, v_var, w_var;
	double gfactor;
	bool print_data;
	bool clamped;
    double growth_prefact;

	/**
	* For time-dependent problems, this is the time t for which the current
	* nonlinear_solution is defined.
	*/
	Real time;
	/**
	* For time-dependent problems, this is the amount delta t to advance the
	* solution in time.
	*/
	Real deltat;

	/**
	 * Flag to tell the system if it should measure now
	 */
	bool measure_now;
	/*void write_measurements(std::string file_name) {std::cout<<"write_measurement: Not implemented in base class\n"; };*/
	void set_measurement_system(ShellMeasurementSystem* sys);
	
	
	Parameters& parameters(){
		return *param_;
	}
	Parameters& parameters(Parameters& param){
		param_ = &param;
		return *param_;
	}
	
    double vol_multiplier;
    double mu_vol;
    double volume;
    double volume0;
	
protected:
	int _extra_order;
	std::vector<ShellCacheObj> ShellCache;
	bool _use_cache;
	int number_of_non_ghosts;
		
	
	/**
	 * Total potential energy of tension and bending
	 */
	double _Etens_total;
	double _Ebend_total;
	double _Eext_total;
	ShellMeasurementSystem* msys;
	
    EquationSystems* eqsys;  //Pointer to the equation system
	
	//pointer to equationsystems->parameters
	Parameters* param_;
};


#endif
