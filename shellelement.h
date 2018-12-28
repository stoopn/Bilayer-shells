/*
 *  utils.h
 *  
 *
 *  Created by Norbert Stoop on 04.04.09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef __shellelement_h__
#define __shellelement_h__ 
 
// Basic include file needed for the mesh functionality.
#include "dense_matrix.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include "dense_vector.h"
#include "vector_value.h"
#include "tensor_value.h"
#include "elem.h"
#include "fe_base.h"
#include "BoundingBox.h"
#include "numeric_vector.h"

// Forward Declarations
//class FEBase;
class ShellSystem;

//template <typename T> class NumericVector;

class ShellCacheObj
{
public:
  ShellCacheObj() {initialized = false; element_fe=NULL; load.zero(); exclude=false;}
  ~ShellCacheObj() { if (initialized) { delete element_fe; element_fe = NULL; initialized=false;  }}
	const Elem* elem;
    FEBase* element_fe;
	std::vector<unsigned int> dof_indices;
	std::vector<unsigned int> dof_indices_u;
	std::vector<unsigned int> dof_indices_v;
	std::vector<unsigned int> dof_indices_w;

	RealVectorValue current_position[3];
	RealVectorValue current_velocity[3];
	bool initialized;
	BoundingBox Bbox;
	RealVectorValue load;
	RealVectorValue normal;
	bool is_ghost;
	bool exclude;
	std::vector<int> VoxelIndices;
	std::vector<Node*> node_patch;
	std::set<int> neighbour_elements;
	std::set<int> next_neighbour_elements;
};


class SurfaceMetric 
{
public:
    void print()
	{
		std::cout<<"Acov[3]: \n";
		acov[0].print(std::cout);
		acov[1].print(std::cout);
		acov[2].print(std::cout);
		std::cout<<"Acov2[3]: \n";
		acov2[0].print(std::cout);
		acov2[1].print(std::cout);
		acov2[2].print(std::cout);
		std::cout<<"amcov, amcon: \n";
		amcov.print(std::cout);
		amcon.print(std::cout);
	}
	SurfaceMetric() 
	{
		amcon.resize(2,2);
		amcov.resize(2,2);
	}
	RealVectorValue acov[3];
    RealVectorValue acon[3];
	RealVectorValue acov2[3];
	DenseMatrix<Real> amcon;
	DenseMatrix<Real> amcov;
    double Concentration;  //Value of the concentration field at the specific quadrature point
	double det;
    Point xyz;

};

class ShellMaterial 
{
public:
	void init() { 
		prefact1 = young*thickness/(1.-poisson*poisson); 
		prefact2 = young*thickness*thickness*thickness/(12.*(1.-poisson*poisson));
	}
	double young;
	double poisson;
	double thickness;
	double prefact1;
	double prefact2;
	double density;
	double contact_stiffness;
}; 


class ShellElement 
{
 public:
	ShellElement(const ShellSystem &sys, int extra_order=1, bool use_cache=false);
	~ShellElement();
	
	void reinit(const Elem *e);
	
	void elem_reinit(); //FIXME:: Not needed
	void elem_fe_reinit(); //FIXME:: Just a wrapper for fe->reinit(elem);
	void set_growth(double kappabar, double L0)
    {
        _gfactor_kappabar = kappabar; _gfactor_L0 = L0;
        if (kappabar!=0. || L0 !=1.)
        {
            is_growing=true;
        }
        else
        {
            is_growing=false;
        }
    }
	
	void add_external_loads(ShellCacheObj* sobj);
    double get_volume(bool rebuild);

	bool request_jacobian;
	DenseVector<Number> mass;
	DenseVector<Number> residual;
	DenseMatrix<Number> jacobian;
	
	DenseSubVector<Number> * subresidual_u;
	DenseSubVector<Number> * subresidual_v;
	DenseSubVector<Number> * subresidual_w;

								
	std::vector<unsigned int> dof_indices;
	std::vector<unsigned int> dof_indices_u;
	std::vector<unsigned int> dof_indices_v;
	std::vector<unsigned int> dof_indices_w;
	
	SurfaceMetric ar;
	SurfaceMetric ac;
	
	void setup_surface_metrics(int qp);
	void calculate_growth_transformation(Point p, double C);
	void assembly(bool request_jacobian, bool request_residual);
	void setupHmat(RealTensorValue& H);
	void enforce_bc(bool request_jacobian, bool request_residual);
	void add_load();
	void element_mass();
	bool measure_energy;
	void set_uniform(bool uniform) {_uniform=uniform; };
    void set_volume_constrained(bool is_constrained) {volume_constrained=is_constrained; };

	double Etens;
	double Ebend;
	double Eext;
    double kGauss;
    double kMean;
    double ElDefArea;
    double ElRefArea;
    
    unsigned int c_var;
    double volume;
    bool is_growing;
    bool volume_constrained;
protected:
	void cross_mat(RealTensorValue& M, RealVectorValue& v);
	RealTensorValue cross_mat(RealVectorValue& v);
	void dyad_prod(RealTensorValue& w, RealVectorValue& u, RealVectorValue& v);
	RealTensorValue dyad_prod(RealVectorValue& u, RealVectorValue& v);
	
	const ShellSystem* sys;
	FEBase* element_fe;
	QBase *element_qrule;
	const Elem *elem;
	double _gfactor_kappabar;
    double _gfactor_L0;
	
	unsigned int n_dofs;
	unsigned int n_var_dofs;
	
	double H(int b, int a, int g, int d);

	DenseMatrix<Real> GI, GIT;	
	RealTensorValue GV;
    unsigned int u_var;
	unsigned int v_var;
	unsigned int w_var;
	
	int gaga;
	
	ShellMaterial material;
	bool _uniform;
	
	double _use_cache;
	int dim;
	FEType fe_type;
	//ShellCacheObj* sobj;
};

#endif
