/*
 *  shellmeasurementsystem.C
 *  
 *
 *  Created by Norbert Stoop on 06.04.09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <fstream>
#include "dof_map.h"
#include "elem.h"
#include "fe_base.h"
#include "equation_systems.h"
#include "face_tri3_sd.h"
#include "libmesh_logging.h"
#include "mesh.h"
#include "numeric_vector.h"
#include "shellmeasurementsystem.h"

ShellMeasurementSystem::ShellMeasurementSystem (EquationSystems& es, const std::string& name, const unsigned int number)
  : Parent(es, name, number)	 
{
	AutoPtr<DofMap> subdivdofmap( new DofMapSubdiv(this->number()));
	std::cout<<"Replacing ShellMeasurementSystem DOF map with Subdiv DOF map\n";
	this->replace_dof_map(subdivdofmap);

}


ShellMeasurementSystem::~ShellMeasurementSystem() { this->clear(); }

void ShellMeasurementSystem::init_data()
{
	std::cout<<"Measurement system being initialized\n";

	const MeshBase& mesh = this->get_mesh();
	
	this->add_variable ("e_tens", FOURTH,SUBDIV);
	e_tens_var = 0;
	this->add_variable ("e_bend", FOURTH,SUBDIV);
	e_bend_var = 1;	
	this->add_variable ("e_contact", FOURTH,SUBDIV);
	e_contact_var = 2;
    this->add_variable ("kGauss", FOURTH,SUBDIV);
    kGauss_var = 3;
    this->add_variable ("kMean", FOURTH,SUBDIV);
    kMean_var = 4;
    this->add_variable ("ElRefArea", FOURTH,SUBDIV);
    El_RefArea_var = 5;
    this->add_variable ("ElDefArea", FOURTH,SUBDIV);
    El_DefArea_var = 6;
    
	Parent::init_data();

	n_dofs = this->get_dof_map().n_dofs();
	n_vars = this->get_dof_map().n_variables();
	n_var_dofs = n_dofs/n_vars;
	build_average_weights();

}

void ShellMeasurementSystem::clear()
{
  std::cout<<"Clearing measurement system\n";
  Parent::clear();
  std::cout<<"Measurement system cleared\n";
}


void ShellMeasurementSystem::add_element_value(const Elem *el, int var, double value)
{
	int dofidx0 =  el->get_node(0)->dof_number(this->number(),var,0);
	int dofidx1 =  el->get_node(1)->dof_number(this->number(),var,0);
	int dofidx2 =  el->get_node(2)->dof_number(this->number(),var,0);
   	solution->add(dofidx0,value);
	solution->add(dofidx1,value);
	solution->add(dofidx2,value);
}

void ShellMeasurementSystem::add_nodal_value(const Node *node, int var, double value)
{
	int dofidx =  node->dof_number(1,var,0);
	solution->add(dofidx,value);
}


void ShellMeasurementSystem::build_average_weights()
{
	
	avg_weights.resize(n_var_dofs);
	
	const MeshBase& mesh = this->get_mesh();
	
	MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
	const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
		
	for ( ; el != end_el; ++el)
	{
		Elem* elem = *el;
		const Tri3SD* sdelem = static_cast<const Tri3SD*>(elem);
		
		if (sdelem->is_ghost())		
			continue;
		
		int dofidx0 =  elem->get_node(0)->dof_number(this->number(),0,0);
		int dofidx1 =  elem->get_node(1)->dof_number(this->number(),0,0);
		int dofidx2 =  elem->get_node(2)->dof_number(this->number(),0,0);

		avg_weights(dofidx0) += 1;
		avg_weights(dofidx1) += 1;
		avg_weights(dofidx2) += 1;
	}		

}


void ShellMeasurementSystem::prepare_measurement()
{
	//First, average the accumulated values with the correct weighting:
	
	for (int j=0; j<n_vars; j++)
	for (int i = 0; i<n_var_dofs; i++)
	{	
		int idx=j*n_var_dofs+i;
		if (avg_weights(i)>0)
			solution->set(idx, (*solution)(idx)/avg_weights(i));
		
	}
    solution->close();
	this->update();  
	solution->close();
   // std::cout<<"measurement system solution closed"<<std::endl;
    //solution->print();
}

void ShellMeasurementSystem::reset_measurements()
{
    solution->close();
	solution->zero();
	this->update();  
	solution->close();
}
