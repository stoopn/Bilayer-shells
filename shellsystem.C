/*
 *  shellsystem.C
 *  
 *
 *  Created by Norbert Stoop on 06.04.09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "dof_map.h"
#include "elem.h"
#include "fe_base.h"
#include "face_tri3_sd.h"
#include "fe_interface.h"
#include "equation_systems.h"
#include "shellelement.h"
#include "shellsystem.h"
#include "libmesh_logging.h"
#include "mesh.h"
#include "numeric_vector.h"
#include "quadrature.h"
#include "sparse_matrix.h"
#include "getpot.h"
#include "parallel.h"
#include "dof_map.h"
#include "dof_map_subdiv.h"


//integration
const std::string ShellSystem::P_QUADRATURE_EXTRA_ORDER = "quadrature_extra_order";
const std::string ShellSystem::P_POISSON_RATIO = "poisson_ratio";
const std::string ShellSystem::P_YOUNG_MODULUS = "young_modulus";
const std::string ShellSystem::P_THICKNESS = "thickness";
const std::string ShellSystem::P_DENSITY = "density";



ShellSystem::ShellSystem (EquationSystems& es, const std::string& name, const unsigned int number)
  : Parent(es, name, number)
{
	param_ = &es.parameters;
	
	//integration
	_extra_order = param_->get<int>(P_QUADRATURE_EXTRA_ORDER);
    eqsys = &es; //Set the pointer
	_use_cache = true;
	AutoPtr<DofMap> subdivdofmap( new DofMapSubdiv(this->number()));
	std::cout<<"Replacing ShellSystem DOF map with Subdiv DOF map\n";
	this->replace_dof_map(subdivdofmap);
	measure_now = false;
 }


//Destructor
ShellSystem::~ShellSystem () {
	param_ = NULL;
	this->clear();
}


void ShellSystem::clear()
{
	Parent::clear();
	for (int cc=0; cc<ShellCache.size(); ++cc)
	{
		if (ShellCache[cc].element_fe != NULL)
			delete ShellCache[cc].element_fe;
		ShellCache[cc].element_fe = NULL;
		ShellCache[cc].initialized = false;
	}
}

void ShellSystem::set_measurement_system(ShellMeasurementSystem* sys)
{
	this->msys = sys;	
}

void ShellSystem::reinit ()
{
  Parent::reinit();
}


void ShellSystem::init_data ()
{
  const unsigned int dim = this->get_mesh().mesh_dimension();

  this->add_variable ("u", FOURTH,SUBDIV);
  u_var = 0;
  this->add_variable ("v", FOURTH,SUBDIV);
  v_var = 1;
  this->add_variable ("w", FOURTH,SUBDIV);
  w_var = 2;

	
  material.poisson	    = param_->get<double>(P_POISSON_RATIO);
  material.thickness    = param_->get<double>(P_THICKNESS); 
  material.young		= param_->get<double>(P_YOUNG_MODULUS); 
  material.density		= param_->get<double>(P_DENSITY); 
  material.contact_stiffness	= 1;
  
  material.init();

	
  Parent::init_data();
  //Shellelem is an AutoPtr and should destroy itself when it is time... Pray it is so...
  shellelement = AutoPtr<ShellElement>(new ShellElement(*this, _extra_order, _use_cache));
  
  if (_use_cache) 
  {
	  const MeshBase& mesh = this->get_mesh();
	  //Count all non-ghosted elements for the cache setup
	  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
	  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
	  number_of_non_ghosts = 0;
	  int numberii=0;
	  for ( ; el != end_el; ++el)
	  {
		  Elem* elem = *el;
		  const Tri3SD* sdelem = static_cast<const Tri3SD*>(elem);
		  ++number_of_non_ghosts;
		  if (!sdelem->is_ghost())
			  ++numberii;
	  }
	  std::cout<<"Number of non-ghosts: "<<numberii<<std::endl;
	  ShellCache.resize(mesh.n_elem());
	  for (int cc=0; cc<ShellCache.size(); ++cc)
	  {
		  ShellCache[cc].element_fe = NULL;
		  ShellCache[cc].initialized = false;
		  ShellCache[cc].elem = mesh.elem(cc);
	  }
  }
    std::cout<<"ShellSystem initialized\n";
 }

void ShellSystem::update_time(double ntime)
{
  this->time = ntime;
}


ShellCacheObj* ShellSystem::get_cache_obj(const Elem *elem) const
{
	libmesh_assert(elem != NULL);
	return const_cast<ShellCacheObj*>(&(ShellCache[elem->id()]));
}

bool ShellSystem::has_cache_obj(const Elem *elem) const
{
	return ShellCache[elem->id()].initialized;
}

void ShellSystem::add_cache_obj(const Elem* elem, ShellCacheObj* obj)
{
	ShellCache[elem->id()] = *obj;
}

