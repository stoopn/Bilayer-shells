/*
 *  utils.h
 *  
 *
 *  Created by Norbert Stoop on 04.04.09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

// Basic include file needed for the mesh functionality.
#include "shellelement.h"
#include "shellsystem.h"
#include "dense_submatrix.h"
#include "dense_subvector.h"
#include "fe.h"
#include "elem.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "boundary_info.h"
#include "libmesh_logging.h"

ShellElement::ShellElement(const ShellSystem &sys, int extra_order, bool use_cache) : 
sys(&sys),
_use_cache(use_cache)
{
	material = sys.material;
	is_growing = false;
    volume_constrained = false;
	_uniform=false;

	u_var=sys.u_var;
	v_var=sys.v_var;
	w_var=sys.w_var;
	material.init();
	
	subresidual_u = new DenseSubVector<Number>(residual);
	subresidual_v = new DenseSubVector<Number>(residual);
	subresidual_w = new DenseSubVector<Number>(residual);

	GI.resize(2,2);
	
	fe_type = sys.variable_type(0);
	this->dim=2;
	
	element_qrule = fe_type.default_quadrature_rule(dim,extra_order).release();
	if (!_use_cache) 
	{
		element_fe = FEBase::build(dim, fe_type).release();
		element_fe->attach_quadrature_rule(element_qrule);
	}
	Etens=0.;
	Ebend=0.;
    kGauss=0;
    kMean=0;
    ElDefArea = 0;
    ElRefArea = 0;
	measure_energy = false;
}


ShellElement::~ShellElement ()
{
	
	delete subresidual_u;
	delete subresidual_v;
	delete subresidual_w;
		

	if (!_use_cache) {
		if (element_fe != NULL)
			delete element_fe;
		element_fe = NULL;
	}
	
	delete element_qrule;
    element_qrule = NULL;
	
}

double ShellElement::get_volume(bool rebuild)
{
    if (rebuild)
    {
        const std::vector< Real > weights = element_qrule->get_weights();
        
        volume=0;
        for (unsigned int qp=0; qp<element_qrule->n_points(); ++qp)
        {
            setup_surface_metrics(qp);
            volume+=weights[qp]*ac.det*ac.acov[2]*ac.xyz*1./3.;
        }
    }
    return volume;
    
}

void ShellElement::add_external_loads(ShellCacheObj* sobj)
{
	const std::vector<Real>& JxW = element_fe->get_JxW();
	const std::vector<std::vector<Real> >& phi = element_fe->get_phi();
	DenseSubVector<Number> &Ru = *subresidual_u;
	DenseSubVector<Number> &Rv = *subresidual_v;
	DenseSubVector<Number> &Rw = *subresidual_w;
	
	if (measure_energy) 
		Eext = 0.;
	
//std::cout<<"Adding the following loads: "<<sobj->load<<std::endl;

	for (unsigned int qp=0; qp<element_qrule->n_points(); ++qp)
	{
		
		for (unsigned int i=0; i<n_var_dofs; i++) 
		{
			Ru(i) += JxW[qp]*phi[i][qp]*sobj->load(0);
			Rv(i) += JxW[qp]*phi[i][qp]*sobj->load(1);
			Rw(i) += JxW[qp]*phi[i][qp]*sobj->load(2);
		}
		if (measure_energy) 
		{
			double work=0.;
			for (unsigned int l=0; l<n_var_dofs; l++) 
			{
				double w = sys->current_solution(dof_indices_w[l]);
				work+=sobj->load(2)*w*phi[l][qp]*JxW[qp];
				
			}
			Eext += work;
		}
	}
	
}


void ShellElement::reinit(const Elem* elem)
{
  START_LOG("reinit", "ShellElement");
	//this->elem=elem;
    //std::cout<<"here31: ID="<<elem->id()<<"\n";

	if (_use_cache && sys->has_cache_obj(elem) )
	{
		ShellCacheObj* cache_obj = sys->get_cache_obj(elem);
		dof_indices = cache_obj->dof_indices;
		dof_indices_u = cache_obj->dof_indices_u;
		dof_indices_v = cache_obj->dof_indices_v;
		dof_indices_w = cache_obj->dof_indices_w;
		element_fe = cache_obj->element_fe;
      //  std::cout<<"here32\n";
	
	} else {
       // std::cout<<"Reinitializing element "<<elem->id()<<"\n";

		sys->get_dof_map().dof_indices(elem, dof_indices);
		sys->get_dof_map().dof_indices(elem, dof_indices_u, u_var);
		sys->get_dof_map().dof_indices(elem, dof_indices_v, v_var);
		sys->get_dof_map().dof_indices(elem, dof_indices_w, w_var);

		// Reinitializing the FE objects is definit ely necessary
      //  std::cout<<"here34\n";

		if (_use_cache) {
			ShellCacheObj* cache_obj = sys->get_cache_obj(elem);
			cache_obj->element_fe = FEBase::build(dim, fe_type).release();
			cache_obj->element_fe->attach_quadrature_rule(element_qrule);
			cache_obj->elem = elem;
        //    std::cout<<"here35\n";

			cache_obj->element_fe->reinit(elem);
			cache_obj->dof_indices = dof_indices;
			cache_obj->dof_indices_u = dof_indices_u;
			cache_obj->dof_indices_v = dof_indices_v;
			cache_obj->dof_indices_w = dof_indices_w;
			cache_obj->initialized = true;
			element_fe = cache_obj->element_fe;
          //  std::cout<<"here36\n";
	
		} else {
          //  std::cout<<"here37\n";
			element_fe->reinit(elem);
		}
	}
	n_var_dofs = dof_indices_u.size(); 
	n_dofs = dof_indices.size();
	
	Etens = 0.;
    kGauss= 0;
    kMean = 0;
	Ebend = 0.;
	Eext =  0.;
    ElDefArea = 0;
    ElRefArea = 0;
	STOP_LOG("reinit", "ShellElement");
  //  std::cout<<"here37\n";

}

double ShellElement::H(int b, int a, int g, int d) 
{
	return material.poisson*ar.amcon(b,a)*ar.amcon(g,d) + 0.5*(1.-material.poisson)*(ar.amcon(b,g)*ar.amcon(a,d)+ar.amcon(b,d)*ar.amcon(a,g));
}


//CHECKED AND CORRECT
void ShellElement::setupHmat(RealTensorValue& H)
{
	H(0,0)          = ar.amcon(0,0)*ar.amcon(0,0);
	H(0,1) = H(1,0) = material.poisson*ar.amcon(0,0)*ar.amcon(1,1)+(1.-material.poisson)*ar.amcon(0,1)*ar.amcon(0,1);
	H(0,2) = H(2,0) = ar.amcon(0,0)*ar.amcon(0,1);
	H(1,1)          = ar.amcon(1,1)*ar.amcon(1,1);
	H(1,2) = H(2,1) = ar.amcon(1,1)*ar.amcon(0,1);
	H(2,2)          = 0.5*( (1.-material.poisson)*ar.amcon(0,0)*ar.amcon(1,1) + (1.+material.poisson)*ar.amcon(0,1)*ar.amcon(0,1) );
}

//CHECKED AND CORRECT
void ShellElement::setup_surface_metrics(int qp)
{
	START_LOG("setup_surface_metrics", "ShellElement");

  //  const std::vector<std::vector<RealGradient> >& dphi = element_fe->get_dphi();
    const std::vector<std::vector<Real> >& phi = element_fe->get_phi();
    const std::vector< Point >& xyz = element_fe->get_xyz();

	const std::vector<std::vector<RealGradient> >& dphi = element_fe->get_dphi();
	const std::vector<std::vector<RealTensor> >& d2phi = element_fe->get_d2phi();
	const std::vector< RealGradient >& dxyzdxi = element_fe->get_dxyzdxi();
	const std::vector< RealGradient >& dxyzdeta = element_fe->get_dxyzdeta();
	const std::vector< RealGradient >& d2xyzdxi2 = element_fe->get_d2xyzdxi2();
	const std::vector< RealGradient >& d2xyzdeta2 = element_fe->get_d2xyzdeta2();
	const std::vector< RealGradient >& d2xyzdxideta = element_fe->get_d2xyzdxideta();
	RealGradient disp, ddispdxi, ddispdeta, d2dispdxi2, d2dispdeta2, d2dispdxideta;
	Real det, invdet;
		
	ar.acov[0] = dxyzdxi[qp];
	ar.acov[1] = dxyzdeta[qp];
	ar.acov[2] = ar.acov[0].cross(ar.acov[1]);
	ar.det = ar.acov[2].size();
	ar.acov[2]/=ar.det;
	
	ar.acov2[0] = d2xyzdxi2[qp];
	ar.acov2[2] = d2xyzdxideta[qp];
	ar.acov2[1] = d2xyzdeta2[qp];
	
	// Covariant metric tensor:       
    ar.amcov(0,0) = ar.acov[0]*ar.acov[0];
    ar.amcov(1,0) = ar.amcov(0,1) = ar.acov[0]*ar.acov[1];
    ar.amcov(1,1) = ar.acov[1]*ar.acov[1];
	//Contravariant metric tensor:
	det = (ar.amcov(0,0)*ar.amcov(1,1) - ar.amcov(1,0)*ar.amcov(1,0));      
	invdet = 1./det;
    ar.amcon(0,0) = invdet*ar.amcov(1,1);
    ar.amcon(1,0) = ar.amcon(0,1) = -invdet*ar.amcov(0,1);
    ar.amcon(1,1) = invdet*ar.amcov(0,0);
	
    ar.xyz = xyz[qp];

	//Ok, now deal with the current configuration's metric. We first need to find the current 
	//solution approximant at the Gauss point gp

	double u,v,w;
	RealVectorValue sv;
	for (unsigned int l=0; l<n_var_dofs; l++)
	{
		//std::cout<<"dof_indices_w["<<l<<"]: "<<dof_indices_w[l]<<", with sol: "<<current_soln(dof_indices_w[l])<<std::endl;
		sv(0) = sys->current_solution(dof_indices_u[l]);
		sv(1) = sys->current_solution(dof_indices_v[l]);
		sv(2) = sys->current_solution(dof_indices_w[l]);
		ddispdxi.add_scaled(sv, dphi[l][qp](0));
		ddispdeta.add_scaled(sv, dphi[l][qp](1));
		d2dispdxi2.add_scaled(sv, d2phi[l][qp](0,0));
		d2dispdeta2.add_scaled(sv, d2phi[l][qp](1,1));
		d2dispdxideta.add_scaled(sv, d2phi[l][qp](0,1));
        disp.add_scaled(sv,phi[l][qp]); //osman
	}
	ac.acov[0] = ar.acov[0] + ddispdxi;
	ac.acov[1] = ar.acov[1] + ddispdeta;
	ac.acov[2] = ac.acov[0].cross(ac.acov[1]);
	ac.det = ac.acov[2].size();
	ac.acov[2]/=ac.det;
	
	ac.acov2[0] = ar.acov2[0] + d2dispdxi2;
	ac.acov2[2] = ar.acov2[2] + d2dispdxideta;
	ac.acov2[1] = ar.acov2[1] + d2dispdeta2;
	
	// Covariant metric tensor:       
	ac.amcov(0,0) = ac.acov[0]*ac.acov[0];
	ac.amcov(1,0) = ac.amcov(0,1) = ac.acov[0]*ac.acov[1];
	ac.amcov(1,1) = ac.acov[1]*ac.acov[1];
    
    ac.xyz = xyz[qp]+disp;

	//Contravariant metric tensor:
    det = (ac.amcov(0,0)*ac.amcov(1,1) - ac.amcov(1,0)*ac.amcov(1,0));
    invdet = 1./det;
    ac.amcon(0,0) = invdet*ac.amcov(1,1);
    ac.amcon(1,0) = ac.amcon(0,1) = -invdet*ac.amcov(0,1);
    ac.amcon(1,1) = invdet*ac.amcov(0,0);
    
    // Contraviariant tangent vectors
    ac.acon[0] = ac.amcon(0,0)*ac.acov[0]+ac.amcon(0,1)*ac.acov[1];
    ac.acon[1] = ac.amcon(1,0)*ac.acov[0]+ac.amcon(1,1)*ac.acov[1];
    
	STOP_LOG("setup_surface_metrics", "ShellElement");
}


void ShellElement::calculate_growth_transformation(Point p, double C)
{
/*
    double x=p(0);
	double y=p(1);
	double time = 1.;
	double k1t=0.0;
	//k1t = _gfactor;
	double k2t=_gfactor;
	if (_uniform)
	  k1t=k2t;

	double r2=x*x+y*y;
*/
	//	double sinp2=y*y/r2;
	//double cosp2=x*x/r2;
	/*double sinp2=1./(1.+(x/y)*(x/y));
	double cosp2=1./(1.+(y/x)*(y/x));
	double cospsinp=x*y/r2;
	double cos2p = cosp2 - sinp2;
	double sin2p = 2.*cospsinp;
	
	double a1x=ar.acov[0](0);
	double a1y=ar.acov[0](1);
	double a2x=ar.acov[1](0);
	double a2y=ar.acov[1](1);
	double tt1=(1.+k1t)*(1.+k2t);
	double tt2=a1x*a2y-a1y*a2x;
	double tti= 1./(tt1*tt2);
	double tt3=(k1t-k2t)*cospsinp;*/
	
	/*	GI(0,0) = tti*(  a1x*(  a2x*tt3 + a2y*((1.+k2t)*cosp2 + (1.+k1t)*sinp2) )
					+ a1y*( -a2y*tt3 - a2x*((1.+k1t)*cosp2 + (1.+k2t)*sinp2) )
					);
	GI(0,1) = tti*(  a2x*(  a2x*tt3 + a2y*((1.+k2t)*cosp2 + (1.+k1t)*sinp2) )
					+ a2y*( -a2y*tt3 - a2x*((1.+k1t)*cosp2 + (1.+k2t)*sinp2) )
					);
	GI(1,0) = tti*(  a1x*( -a1x*tt3 - a1y*((1.+k2t)*cosp2 + (1.+k1t)*sinp2) )
					+ a1y*(  a1y*tt3 + a1x*((1.+k1t)*cosp2 + (1.+k2t)*sinp2) )
					);
	GI(1,1) = tti*(  a2x*( -a1x*tt3 - a1y*((1.+k2t)*cosp2 + (1.+k1t)*sinp2) )
					+ a2y*(  a1y*tt3 + a1x*((1.+k1t)*cosp2 + (1.+k2t)*sinp2) )
					);
	*/

	//Note this is the inverse of G!
	/*
     
     GI(0,0) = (1.+growth_prefact*C);
	GI(0,1) = 0;

	GI(1,0) = 0;
	GI(1,1) = (1.+growth_prefact*C);

*/
    
	//Debugging: Grow in patch coordinate system:
	/*	GI(0,0) = 1./(1.+k1t);
	GI(0,1) = 0.;
	GI(1,0) = 0.;
	GI(1,1) = 1./(1.+k2t);
	*/
	//GV is the matrix in Voigt Notation for the transformation G^-T X G^-1 for some matrix X.
	//This can be written as GV*XV, where XV is X as a matrix operator in Voigt notation and GV as below. We are allowed to
	//do that becase the above transformation still results in a symmetric matrix.	
	/*
	GV(0,0) = GI(0,0)*GI(0,0);
	GV(0,1) = GI(1,0)*GI(1,0);
	GV(0,2) = 2.*GI(0,0)*GI(1,0);
	
	GV(1,0) = GI(0,1)*GI(0,1);
	GV(1,1) = GI(1,1)*GI(1,1);
	GV(1,2) = 2.*GI(0,1)*GI(1,1);
	
	GV(2,0) = GI(0,0)*GI(0,1);
	GV(2,1) = GI(1,0)*GI(1,1);
	GV(2,2) = GI(0,1)*GI(1,0) + GI(0,0)*GI(1,1);
*/
//	std::cout<<"Growth GV:\n";
//	GV.print(std::cout);
}


//FIXME:: Currently only for clamped systems
void ShellElement::enforce_bc(bool get_residual, bool get_jacobian) 
{
	/*
	const Real penalty = 50000.;
	double vvalue=0.0;
	const std::vector<Real>& JxW = element_fe->get_JxW();
	const std::vector<std::vector<Real> >& phi = element_fe->get_phi();
	const MeshBase& mesh = sys->get_mesh();
	
	std::cout<<"Enforcing BCs: Don't do it like this (shellelement.C)";
	exit(0);
	
	DenseSubVector<Number> &u = *subsolution_u;
	DenseSubVector<Number> &v = *subsolution_v;
	DenseSubVector<Number> &w = *subsolution_w;
	DenseSubVector<Number> &Ru = *subresidual_u;
	DenseSubVector<Number> &Rv = *subresidual_v;
	DenseSubVector<Number> &Rw = *subresidual_w;

	DenseSubMatrix<Number> &Kuu = *subjacobian_uu;
	DenseSubMatrix<Number> &Kvv = *subjacobian_vv;
	DenseSubMatrix<Number> &Kww = *subjacobian_ww;
	int edge_found=0;
	const Tri3SD* sdelem = static_cast<const Tri3SD*>(elem);
	int	val=sdelem->get_ordered_node(0)->valence();
	*/
			//Debugging: 
/*	if (get_jacobian)
	{
		for (unsigned int i=0; i<n_var_dofs; i++) 
		{
		Kww(i,i) += penalty;
		}
	}
*/
	//IMPORTANT: I think for nonlinear FE, enforcing BCs by penalties needs addition of the penalty to the diagonal
	// *plus* penalty*(current_value - desired_value) in the residual. The linearizes system J*x=R then finds x=current_value-desired_value
	// The next Newton iterate solution s then has s -> s - x, with s=current_solution, ie the new s will be s = current_value - current_value + desired_value
	// In other words, the linear constraint is fulfilled in each linear step exactly, as it should be. 
/*	for (unsigned int s=0; s<elem->n_sides(); s++)
	{
		//short int bcid = mesh.boundary_info->boundary_id(elem->get_node(s));
		Tri3SD* nbelem = static_cast<Tri3SD*>(elem->neighbor(s));
		
		if (!nbelem->is_ghost())
			continue;
				
		if (get_residual)
		{
			//Debugging
			
*/		
			/*	//The first node of the triangle
				Ru(0) += penalty*(u(0)-vvalue); 
				Rv(0) += penalty*(v(0)-vvalue); 
				Rw(0) += penalty*(w(0)-vvalue); 
				//The second node of the triangle
				Ru(1) += penalty*(u(1)-vvalue); 
				Rv(1) += penalty*(v(1)-vvalue); 
				Rw(1) += penalty*(w(1)-vvalue); 
				//The third node of the triangle
				Ru(val) += penalty*(u(val)-vvalue); 
				Rv(val) += penalty*(v(val)-vvalue); 
				Rw(val) += penalty*(w(val)-vvalue); 
			*/	// Find out on which side the ghost sits:
	/*		
				if (sdelem->get_node(s) == sdelem->get_ordered_node(0)) 
				{
					//The first node of the triangle
					Ru(0) += penalty*(u(0)-vvalue); 
					Rv(0) += penalty*(v(0)-vvalue); 
					Rw(0) += penalty*(w(0)-vvalue); 
					//The second node of the triangle
					Ru(1) += penalty*(u(1)-vvalue); 
					Rv(1) += penalty*(v(1)-vvalue); 
					Rw(1) += penalty*(w(1)-vvalue); 
			
				//The node on the ghost has to be equal to the val node (the interior one)
					Ru(2) += penalty*(u(2)-(-u(val))); 
					Rv(2) += penalty*(v(2)-(-v(val))); 
					Rw(2) += penalty*(w(2)-(-w(val))); 
				//	Ru(val) += -penalty*(u(2)-(-u(val))); 
				//	Rv(val) += -penalty*(v(2)-(-v(val))); 
				//	Rw(val) += -penalty*(w(2)-(-w(val))); 
	*/				
			/*
					Ru(2) += penalty*(u(2)-vvalue); 
					Rv(2) += penalty*(v(2)-vvalue); 
					Rw(2) += penalty*(w(2)-vvalue);
			*/
	/*			} else if (sdelem->get_node(s) == sdelem->get_ordered_node(1)) 
				{
					Ru(1) += penalty*(u(1)-vvalue); 
					Rv(1) += penalty*(v(1)-vvalue); 
					Rw(1) += penalty*(w(1)-vvalue);
					Ru(val) += penalty*(u(val)-vvalue); 
					Rv(val) += penalty*(v(val)-vvalue); 
					Rw(val) += penalty*(w(val)-vvalue);
				
					Ru(val+1) += penalty*(u(val+1)-vvalue); 
					Rv(val+1) += penalty*(v(val+1)-vvalue); 
					Rw(val+1) += penalty*(w(val+1)-vvalue); 
				//	Ru(0) += -penalty*(u(val+1)-(-u(0))); 
				//	Rv(0) += -penalty*(v(val+1)-(-v(0))); 
				//	Rw(0) += -penalty*(w(val+1)-(-w(0))); 
	*/			
				/*
					Ru(val+1) += penalty*(u(val+1)-vvalue); 
					Rv(val+1) += penalty*(v(val+1)-vvalue); 
					Rw(val+1) += penalty*(w(val+1)-vvalue); 
				*/
	/*			} else 
				{
					Ru(val) += penalty*(u(val)-vvalue); 
					Rv(val) += penalty*(v(val)-vvalue); 
					Rw(val) += penalty*(w(val)-vvalue);
					Ru(0) += penalty*(u(0)-vvalue); 
					Rv(0) += penalty*(v(0)-vvalue); 
					Rw(0) += penalty*(w(0)-vvalue); 
					
					//Ghost sits on the left of the triangle
					Ru(val-1) += penalty*(u(val-1)-vvalue); 
					Rv(val-1) += penalty*(v(val-1)-vvalue); 
					Rw(val-1) += penalty*(w(val-1)-vvalue); 
				//	Ru(1) += -penalty*(u(val-1)-(-u(1))); 
				//	Rv(1) += -penalty*(v(val-1)-(-v(1))); 
				//	Rw(1) += -penalty*(w(val-1)-(-w(1))); 
	*/				
				/*
					//Ghost sits on the left of the triangle
					Ru(val-1) += penalty*(u(val-1)-vvalue); 
					Rv(val-1) += penalty*(v(val-1)-vvalue); 
					Rw(val-1) += penalty*(w(val-1)-vvalue); 
				*/
	/*			}
		}
		if (get_jacobian)
		{
	*/		
			/*	//The first node of the triangle
		 		Kuu(0,0) += penalty;			
				Kvv(0,0) += penalty;
				Kww(0,0) += penalty;
				//The second node of the triangle
				Kuu(1,1) += penalty;
				Kvv(1,1) += penalty;
				Kww(1,1) += penalty;
				//The third node of the triangle
				Kuu(val,val) += penalty;
				Kvv(val,val) += penalty;
				Kww(val,val) += penalty;
			*/	// Find out on which side the ghost sits:	
	/*		
			if (sdelem->get_node(s) == sdelem->get_ordered_node(0)) 
				{
					Kuu(0,0) += penalty;			
					Kvv(0,0) += penalty;
					Kww(0,0) += penalty;
					Kuu(1,1) += penalty;
					Kvv(1,1) += penalty;
					Kww(1,1) += penalty;
			
					Kuu(2,2) += penalty;
					Kvv(2,2) += penalty;
					Kww(2,2) += penalty;
				//	Kuu(val,val) += penalty;
				//	Kvv(val,val) += penalty;
				//	Kww(val,val) += penalty;
					
				//	Kuu(2,val) += -penalty;
				//	Kvv(2,val) += -penalty;
				//	Kww(2,val) += -penalty;
				//	Kuu(val,2) += -penalty;
				//	Kvv(val,2) += -penalty;
				//	Kww(val,2) += -penalty;
			
				} else if (sdelem->get_node(s) == sdelem->get_ordered_node(1)) 
				{
					Kuu(1,1) += penalty;
					Kvv(1,1) += penalty;
					Kww(1,1) += penalty;
					Kuu(val,val) += penalty;
					Kvv(val,val) += penalty;
					Kww(val,val) += penalty;
			
					Kuu(val+1,val+1) += penalty;
					Kvv(val+1,val+1) += penalty;
					Kww(val+1,val+1) += penalty;
				//	Kuu(0,0) += penalty;
				//	Kvv(0,0) += penalty;
				//	Kww(0,0) += penalty;
					
				//	Kuu(val+1,0) += -penalty;
				//	Kvv(val+1,0) += -penalty;
				//	Kww(val+1,0) += -penalty;
				//	Kuu(0,val+1) += -penalty;
				//	Kvv(0,val+1) += -penalty;
				//	Kww(0,val+1) += -penalty;					
			
				} else 
				{
					Kuu(0,0) += penalty;			
					Kvv(0,0) += penalty;
					Kww(0,0) += penalty;
					Kuu(val,val) += penalty;
					Kvv(val,val) += penalty;
					Kww(val,val) += penalty;
					
					Kuu(val-1,val-1) += penalty;
					Kvv(val-1,val-1) += penalty;
					Kww(val-1,val-1) += penalty;
				//	Kuu(1,1) += penalty;
				//	Kvv(1,1) += penalty;
				//	Kww(1,1) += penalty;
					
				//	Kuu(val-1,1) += -penalty;
				//	Kvv(val-1,1) += -penalty;
				//	Kww(val-1,1) += -penalty;
				//	Kuu(1,val-1) += -penalty;
				//	Kvv(1,val-1) += -penalty;
				//	Kww(1,val-1) += -penalty;					
			
				}
			
			
		}
		
	}

*/	
}


	
//Converts a cross product of v with some other vector into a matrix multiplication
//with the matrix M and the vector. 
void ShellElement::cross_mat(RealTensorValue& M, RealVectorValue& v)
{
	M(0,0) =    0.;  M(0,1) = -v(2);  M(0,2) =  v(1);
	M(1,0) =  v(2);  M(1,1) =    0.;  M(1,2) = -v(0);
	M(2,0) = -v(1);  M(2,1) =  v(0);  M(2,2) =    0.;
}
RealTensorValue ShellElement::cross_mat(RealVectorValue& v)
{
	RealTensorValue M;
	M(0,0) =    0.;  M(0,1) = -v(2);  M(0,2) =  v(1);
	M(1,0) =  v(2);  M(1,1) =    0.;  M(1,2) = -v(0);
	M(2,0) = -v(1);  M(2,1) =  v(0);  M(2,2) =    0.;
	return M;
}


//Constructs the dyadic product w = u x v where x is the tensor product
void ShellElement::dyad_prod(RealTensorValue& w, RealVectorValue& u, RealVectorValue& v)
{
	w(0,0) = u(0)*v(0); w(0,1) = u(0)*v(1); w(0,2) = u(0)*v(2); 
	w(1,0) = u(1)*v(0); w(1,1) = u(1)*v(1); w(1,2) = u(1)*v(2); 
	w(2,0) = u(2)*v(0); w(2,1) = u(2)*v(1); w(2,2) = u(2)*v(2); 
}
RealTensorValue ShellElement::dyad_prod(RealVectorValue& u, RealVectorValue& v)
{
	RealTensorValue w;
	w(0,0) = u(0)*v(0); w(0,1) = u(0)*v(1); w(0,2) = u(0)*v(2); 
	w(1,0) = u(1)*v(0); w(1,1) = u(1)*v(1); w(1,2) = u(1)*v(2); 
	w(2,0) = u(2)*v(0); w(2,1) = u(2)*v(1); w(2,2) = u(2)*v(2); 
	return w;
}



void ShellElement::assembly(bool get_residual, bool get_jacobian)
{
	START_LOG("assembly", "ShellElement");
	
	const std::vector< Point>& xyz = element_fe->get_xyz();
	const std::vector<Real>& JxW = element_fe->get_JxW();
	const std::vector<std::vector<Real> >& phi = element_fe->get_phi();
	const std::vector<std::vector<RealGradient> >& dphi = element_fe->get_dphi();
	const std::vector<std::vector<RealTensor> >& d2phi = element_fe->get_d2phi();
    const std::vector< Real > weights = element_qrule->get_weights();

	if (get_residual) {
		residual.resize(n_dofs);
		subresidual_u->reposition(u_var*n_var_dofs, n_var_dofs);
		subresidual_v->reposition(v_var*n_var_dofs, n_var_dofs);
		subresidual_w->reposition(w_var*n_var_dofs, n_var_dofs);
	}

	if (measure_energy) {
		Etens = 0.;
		Ebend = 0.;
        kGauss = 0;
        kMean = 0;
        ElDefArea = 0;
        ElRefArea = 0;
	}
	
	DenseSubVector<Number> &Ru = *subresidual_u;
	DenseSubVector<Number> &Rv = *subresidual_v;
	DenseSubVector<Number> &Rw = *subresidual_w;

	RealVectorValue a1Ca2, a2Ca3, a3Ca1, vb11A, vb11B, vb22A, vb22B, vb12A, vb12B;
	
	RealTensorValue HMAT; //H in Voigt-Notation
	RealVectorValue NF; //nf in Voigt-Notation
	RealVectorValue MF; //mf in Voigt-Notation
	RealVectorValue EPS; //eps in Voigt-Notation (membrane strains)
	RealVectorValue K; //rho in Voigt-Notation (bending strains)
	RealTensorValue BmT; //Membrane strain differential operator
	RealTensorValue BbT; //Bending strain differential operator
    RealVectorValue EPSD, KD, ABAR, BBAR;
	RealVectorValue ResM, ResB;
	
    volume=0;
    for (unsigned int qp=0; qp<element_qrule->n_points(); ++qp)
	{
		setup_surface_metrics(qp);
        //Determine new element volume contribution and volume constraint factors:
        volume+=JxW[qp]*ac.acov[2]*ac.xyz*1./3.;
        RealVectorValue NN0, NN1, FEVect;
        if (this->volume_constrained)
        {
            double v_multiplier =  - sys->mu_vol*(sys->volume - sys->volume0);
            //std::cout<<"Multiplier: "<<v_multiplier<<std::endl;
            //multiplier term
            NN0 =  - v_multiplier/3.*( (ac.xyz*ac.acov[2])*ac.acon[0] - (ac.xyz*ac.acon[0])*ac.acov[2] );
            NN1 =  - v_multiplier/3.*( (ac.xyz*ac.acov[2])*ac.acon[1] - (ac.xyz*ac.acon[1])*ac.acov[2] );
            FEVect = v_multiplier/3.*ac.acov[2];
        }
        
		setupHmat(HMAT);
	
        double kappabar=_gfactor_kappabar;
        double L0=_gfactor_L0;  // \Lambda_0  // We always set m=1
        ABAR(0) = L0*L0*ar.amcov(0,0);
        ABAR(1) = L0*L0*ar.amcov(1,1);
        ABAR(2) = L0*L0*ar.amcov(0,1);
        
        BBAR(0) = L0*ar.acov2[0]*ar.acov[2] - L0*L0*kappabar*ar.amcov(0,0);
        BBAR(1) = L0*ar.acov2[1]*ar.acov[2] - L0*L0*kappabar*ar.amcov(1,1);
        BBAR(2) = L0*ar.acov2[2]*ar.acov[2] - L0*L0*kappabar*ar.amcov(1,0);
        
		
        K(0) = BBAR(0) - ac.acov2[0]*ac.acov[2];
        K(1) = BBAR(1) - ac.acov2[1]*ac.acov[2];
        K(2) = BBAR(2) - ac.acov2[2]*ac.acov[2];
        K(2) *= 2.; //Voigt notation
        
        EPS(0) = ac.amcov(0,0)-ABAR(0);
        EPS(1) = ac.amcov(1,1)-ABAR(1);
        EPS(2) = ac.amcov(0,1)-ABAR(2);
        EPS(2) *= 2.; //Voigt notation
        EPS *= 0.5; //stretching strains have prefact 1/2
        
        KD(0) = L0*L0*material.prefact2*K(0);
        KD(1) = L0*L0*material.prefact2*K(1);
        KD(2) = L0*L0*material.prefact2*K(2);

		EPSD(0) = material.prefact1*EPS(0);
		EPSD(1) = material.prefact1*EPS(1);
		EPSD(2) = material.prefact1*EPS(2);
		
		NF = (HMAT*EPSD);
		MF = (HMAT*KD);

		
		//NF and MF are in Voigt notation (nf_11, nf_22, nf_12), ie. without factor 2 in third component. This is due to the form of HMAT.

		if (measure_energy)
		{
		            
            const Real w = ac.det * weights[qp];
            const Real inv_det = 1. / (ac.amcov(0,0) * ac.amcov(1,1) - ac.amcov(0,1)*ac.amcov(1,0));
            kGauss +=  w*inv_det * (ac.acov2[0]*ac.acov[2] * ac.acov2[1]*ac.acov[2] - ac.acov2[2]*ac.acov[2]*ac.acov2[2]*ac.acov[2]);
            kMean +=  w*inv_det * 0.5 * (ac.acov2[0]*ac.acov[2] * ac.amcov(1,1) + ac.acov2[1]*ac.acov[2] * ac.amcov(0,0) - 2*ac.acov2[2]*ac.acov[2] * ac.amcov(0,1));
            
            ElDefArea += w;
            ElRefArea += JxW[qp];

			Etens += JxW[qp]*0.5*NF*EPS;
			Ebend += JxW[qp]*0.5*MF*K;
		}
	
		//We can precalculate the vectorial part of the bending discrete operator. We need it in any case for
		//the residual and the material part of the  bending jacobian
		
			a2Ca3 = ac.acov[1].cross(ac.acov[2]);
			a3Ca1 = ac.acov[2].cross(ac.acov[0]);
		
			vb11A = (ac.acov2[0].cross(ac.acov[1]) + (ac.acov2[0]*ac.acov[2])*a2Ca3)/ac.det;
			vb11B = (ac.acov[0].cross(ac.acov2[0]) + (ac.acov2[0]*ac.acov[2])*a3Ca1)/ac.det;
	
			vb22A = (ac.acov2[1].cross(ac.acov[1]) + (ac.acov2[1]*ac.acov[2])*a2Ca3)/ac.det;
			vb22B = (ac.acov[0].cross(ac.acov2[1]) + (ac.acov2[1]*ac.acov[2])*a3Ca1)/ac.det;

			vb12A = (ac.acov2[2].cross(ac.acov[1]) + (ac.acov2[2]*ac.acov[2])*a2Ca3)/ac.det;
			vb12B = (ac.acov[0].cross(ac.acov2[2]) + (ac.acov2[2]*ac.acov[2])*a3Ca1)/ac.det;

		//Loop over all basis functions (or the patch of the subdivs)					
		for (unsigned int i=0; i<n_var_dofs; i++) 
		{
			//Okay, now lets use Simo&Fox' elegant formulation of the discrete operators. We nned them for the residual
			//and for the material part of the membrane jacobian

			RealVectorValue row;
			row = ac.acov[0]*dphi[i][qp](0);
			BmT(0,0) = row(0); BmT(1,0) = row(1); BmT(2,0) = row(2);
			row = ac.acov[1]*dphi[i][qp](1);
			BmT(0,1) = row(0); BmT(1,1) = row(1); BmT(2,1) = row(2);
			row = 0.5*(ac.acov[0]*dphi[i][qp](1) + ac.acov[1]*dphi[i][qp](0));
			BmT(0,2) = row(0); BmT(1,2) = row(1); BmT(2,2) = row(2);
			
			BmT(0,2)*=2.; BmT(1,2)*=2.; BmT(2,2)*=2.;
			
			//The same for the bending:
			row = - ac.acov[2]*d2phi[i][qp](0,0) + (dphi[i][qp](0)*vb11A + dphi[i][qp](1)*vb11B);
			BbT(0,0) = row(0); BbT(1,0) = row(1); BbT(2,0) = row(2);
			row = - ac.acov[2]*d2phi[i][qp](1,1) + (dphi[i][qp](0)*vb22A + dphi[i][qp](1)*vb22B);
			BbT(0,1) = row(0); BbT(1,1) = row(1); BbT(2,1) = row(2);
			//Check again if this is really correct (ie. the factor of 2!)
			row = - ac.acov[2]*d2phi[i][qp](0,1) + (dphi[i][qp](0)*vb12A + dphi[i][qp](1)*vb12B);
		//	BbT(0,2) = 2.*row(0); BbT(1,2) = 2.*row(1); BbT(2,2) = 2.*row(2);
			BbT(0,2) = row(0); BbT(1,2) = row(1); BbT(2,2) = row(2);
			
			BbT(0,2)*=2.; BbT(1,2)*=2.; BbT(2,2)*=2.;

			
			if (get_residual) 
			{
                ResM=BmT*NF;
				ResB=BbT*MF;
				
				//Total residual:
				//CKECKED AND CORRECT
                if (this->volume_constrained)
                {
                    RealVectorValue ResN = NN0*dphi[i][qp](0) + NN1*dphi[i][qp](1) - FEVect*phi[i][qp];
                    Ru(i) += JxW[qp]*(ResM(0)+ResB(0))+ResN(0)*ac.det*weights[qp];
                    Rv(i) += JxW[qp]*(ResM(1)+ResB(1))+ResN(1)*ac.det*weights[qp];
                    Rw(i) += JxW[qp]*(ResM(2)+ResB(2))+ResN(2)*ac.det*weights[qp];
                } else {
                    Ru(i) += JxW[qp]*(ResM(0)+ResB(0));
                    Rv(i) += JxW[qp]*(ResM(1)+ResB(1));
                    Rw(i) += JxW[qp]*(ResM(2)+ResB(2));
                }
            }
		}
	}
	STOP_LOG("assembly", "ShellElement");
}




void ShellElement::add_load() 
{
	const std::vector<Real>& JxW = element_fe->get_JxW();
	const std::vector<std::vector<Real> >& phi = element_fe->get_phi();
	DenseSubVector<Number> &Rw = *subresidual_w;

	double Qload=1.0;
	for (unsigned int qp=0; qp<element_qrule->n_points(); ++qp)
	{

		for (unsigned int i=0; i<n_var_dofs; i++) 
		{
			Rw(i) += JxW[qp]*phi[i][qp]*Qload;
		}
	}
}

void ShellElement::element_mass() 
{
	
	const std::vector< Point>& xyz = element_fe->get_xyz();
	const std::vector<Real>& JxW = element_fe->get_JxW();
	const std::vector<std::vector<Real> >& phi = element_fe->get_phi();
	
	mass.resize(n_dofs);
	
	DenseSubVector<Number> Mu(mass);
	DenseSubVector<Number> Mv(mass);
	DenseSubVector<Number> Mw(mass);
	
	Mu.reposition(u_var*n_var_dofs, n_var_dofs);
	Mv.reposition(v_var*n_var_dofs, n_var_dofs);
	Mw.reposition(w_var*n_var_dofs, n_var_dofs);
	
	for (unsigned int qp=0; qp<element_qrule->n_points(); ++qp)
	{
		
		for (unsigned int i=0; i<n_var_dofs; i++) 
		{
			double massg=0.;
			for (unsigned int j=0; j<n_var_dofs; j++) 
				massg += phi[i][qp]*phi[j][qp];	
				
			Mu(i) += JxW[qp]*material.density*material.thickness*massg;
			Mv(i) += JxW[qp]*material.density*material.thickness*massg;
			Mw(i) += JxW[qp]*material.density*material.thickness*massg;
		}
		
	}
}
