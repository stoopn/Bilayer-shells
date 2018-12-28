/*
 *  shellsystem.h
 *  
 *
 *  Created by Norbert Stoop on 06.04.09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
  

#ifndef __shellmeasurementsystem_h__
#define __shellmeasurementsystem_h__

// Local Includes
#include "libmesh_common.h"
#include "explicit_system.h"
#include "numeric_vector.h"
#include <fstream>
#include "mesh_tools.h"
#include "serial_mesh.h"
#include "equation_systems.h"
#include "dof_map.h"
#include "dof_map_subdiv.h"
#include <dense_vector.h>
#include <dense_subvector.h>


// Forward Declarations

//template <typename T> class NumericVector;

class ShellMeasurementSystem : public ExplicitSystem
{
public:
	
	ShellMeasurementSystem (EquationSystems& es, const std::string& name, const unsigned int number);
	virtual ~ShellMeasurementSystem ();

	typedef ShellMeasurementSystem sys_type;
	typedef ExplicitSystem Parent;
  
	DenseVector<Real> avg_weights;	
	void init_data ();
	void clear();	
	unsigned int e_tens_var, e_bend_var, e_contact_var, mass_var, kMean_var, kGauss_var, El_DefArea_var, El_RefArea_var;
	void build_average_weights();
	void add_element_value(const Elem *el, int var, double value);
	void add_nodal_value(const Node *node, int var, double value);

	void prepare_measurement();
	void reset_measurements();

	int n_dofs;
	int n_vars;
	int n_var_dofs;
};


#endif
