
#include <iostream>
#include <algorithm>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh.h"
#include "mesh.h"
#include "mesh_generation.h"
#include "mesh_tools.h"
#include "mesh_subdiv_support.h"
#include "mesh_modification.h"
#include "mesh_refinement.h"
#include "equation_systems.h"
#include "fe.h"
#include "dof_map.h"
#include "dof_map_subdiv.h"
#include "utility.h"
#include "getpot.h"

// For systems of equations the \p DenseSubMatrix
// and \p DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "ExplicitShellsystem.h"
#include "vtk_io.h"
#include <fstream>
#include "perf_log.h"
#include "shell_util.h"
#include <signal.h>

//Mathematica definitions for code export
#include "mdefs.h"

#include "petsc_macro.h"
EXTERN_C_FOR_PETSC_BEGIN
# include <petscts.h>
EXTERN_C_FOR_PETSC_END

//apply random initial conditions
void random_initial_conditions(EquationSystems& es, const std::string& system_name) {
    libmesh_assert(system_name == "Shell");
    
    ShellSystem & system = es.get_system<ShellSystem>("Shell");
    MeshBase& mesh = system.get_mesh();
    
    for (MeshBase::node_iterator it=mesh.nodes_begin() ; it!=mesh.nodes_end(); ++it) {
        Node* node=*it;
        
        system.solution->set(node->dof_number(0,0,0),(drand48()-0.5)*es.parameters.get<double>("ic_rnd_u"));
        system.solution->set(node->dof_number(0,1,0),(drand48()-0.5)*es.parameters.get<double>("ic_rnd_v"));
        system.solution->set(node->dof_number(0,2,0),(drand48()-0.5)*es.parameters.get<double>("ic_rnd_w"));
    }
    
    std::cout<<"System perturbed with random initial conditions!"<<std::endl;
}

EquationSystems* EQS_PTR=NULL;

void save_SIGNAL(int sig)
{
  
  OStringStream file_name;
  if (!EQS_PTR)
    {
      std::cout<<"Aborting.."<<std::endl;
      exit(0);
      return;
    }
  std::cout<<"Aborting... Checkpointing data\n";
  const ExplicitShellSystem& system = EQS_PTR->get_system<ExplicitShellSystem> ("Shell");

  file_name.str("");
  file_name<<EQS_PTR->parameters.get<std::string>("Data dir");
  file_name<<"system_dump.xdr";
//  system.write(file_name.str(), libMeshEnums::WRITE);
  std::cout<<"Dump written. Exiting normally\n";
  exit(0);
  return;
}

void project_to_1_sphere(MeshBase& mesh)
{
    MeshBase::const_node_iterator           nd = mesh.nodes_begin();
    const MeshBase::const_node_iterator end_nd = mesh.nodes_end();
    for ( ; nd != end_nd; ++nd)
    {
        Node* node = *nd;
        *node /=  node->size();
    }
}


// The main program.
int main (int argc, char** argv)
{

	LibMeshInit init (argc, argv);
    GetPot infile("sh.in");
    GetPot command_line(argc, argv);
    bool read_solution=false;
    double starttime=0.0;
    double time;
    bool ghosted, resolve_contacts;
    int plot_each, measure_each;

    if (command_line.search(1,"-read_solution"))
    {
        read_solution=true;
    }
    
    std::string input_mesh = infile("Input_Mesh", "input.off");
    double dt = infile("dt", 0.1);
    int Nsteps = infile("nsteps", 1000);
    std::srand(infile("seed", 0));
    
    const unsigned int dim = 2;
    Mesh mesh(dim);
    mesh.read(input_mesh);

    MeshRefinement mesh_refinement (mesh);
    mesh_refinement.uniformly_refine(infile("nrefine",0));
    MeshTools::Modification::flatten(mesh);
    MeshTools::Modification::all_tri(mesh);
    
    double voxelrad;
    if (infile("project_to_sphere",0)) {
        project_to_1_sphere(mesh);
    voxelrad=1;
    }
    
    const double scalef = infile("scalef", 100);
    
    ghosted = infile("ghosted",false);
    std::cout<<"Scaling mesh by factor scalef="<<scalef<<std::endl;
    MeshTools::Modification::scale(mesh, scalef, scalef, scalef);
    MeshTools::Subdiv::prepare_subdiv_mesh(mesh,ghosted);
    mesh.print_info();
    std::cout<<"D1\n"<<std::endl;

	EquationSystems equation_systems(mesh);
    std::cout<<"D2\n"<<std::endl;

	EQS_PTR = &equation_systems;

    equation_systems.parameters.set<Real>("container_radius") = voxelrad*scalef;
	equation_systems.parameters.set<Real>("dt")              = infile("dt", 1.0e-5);
	equation_systems.parameters.set<Real>("newmark_damping") = infile("Newmark_Damping", 0.1);
	equation_systems.parameters.set<Real>("newmark_gamma")   = infile("Newmark_Gamma", 0.5);
    equation_systems.parameters.set<Real>("newmark_beta")   = infile("Newmark_Beta", 0.25);
   
    //Volume penalty and quench:
    double penalty_mu_vol              = infile("mu_vol", 1.0);
    double reduced_volume_rate         = infile("reduced_volume_rate", 0.0);
    double max_reduced_volume = infile("max_reduced_volume", 1.0);
    double min_reduced_volume = infile("min_reduced_volume", 1.0);
    double reduced_volume_quench_start_time = infile("reduced_volume_quench_start_time", 0.0);
    
    //growth parameters
    double growth_rate_kappabar   = infile("growth_rate_kappabar", 0.0);
    double growth_rate_L0   = infile("growth_rate_L0", 0.0);
    double max_kappabar   = infile("max_kappabar", 10.0);
    double max_L0   = infile("max_L0", 10.0);
    double min_L0   = infile("min_L0", -10.0);
    double min_kappabar   = infile("min_kappabar", -10.0);
    double kappabar= infile("kappabar", 0.0);
    double L0= infile("L0", 1);
    bool fixed_growth=false;
    double growth_factor_kappabar=kappabar;
    double growth_factor_L0=L0;
    equation_systems.parameters.set<double>("growth_factor_kappabar") = growth_factor_kappabar;
    equation_systems.parameters.set<double>("growth_factor_L0") = L0; //start without growth;
    
    //random initial conditions
    equation_systems.parameters.set<double>("ic_rnd_u")						=	infile("ic_rnd_u", 0.001);
    equation_systems.parameters.set<double>("ic_rnd_v")						=	infile("ic_rnd_v", 0.001);
    equation_systems.parameters.set<double>("ic_rnd_w")						=	infile("ic_rnd_w", 0.001);
    
    std::string data_dir = infile("datadir","./");
    equation_systems.parameters.set<std::string>("datadir")  = data_dir;
    equation_systems.parameters.set<int>("quadrature_extra_order") = infile("Quadrature_Extra_Order", 0);
    equation_systems.parameters.set<Real>("poisson_ratio")              = infile("Poisson_Ratio", 0.5);
    equation_systems.parameters.set<Real>("young_modulus")              = infile("Young_Modulus", 1.0);
    equation_systems.parameters.set<Real>("thickness")              = infile("Thickness", 1.0);
    equation_systems.parameters.set<Real>("density")              = infile("Density", 1.0);

    equation_systems.parameters.set<bool>("ghosted")						=	ghosted;
    equation_systems.parameters.set<int>("bc_type")						=	0;			//0=free, 1=clamped
    equation_systems.parameters.set<int>("ghost_bc_type")						=	0;			//0=prop, 1=mirrored, 2=softmirrored
    equation_systems.parameters.set<double>("ghost_bc_stiffness")					=	0.0;
    equation_systems.parameters.set<bool>("resolve_contacts")					=	infile("Resolve_Contacts", false);
    
   
    
    equation_systems.parameters.set<std::string>("measurement_filename")				=	data_dir+"measurement.dat";


    plot_each = infile("ploteach",1);
    measure_each = infile("Measure_Each", 50);

    std::cout<<"D3\n"<<std::endl;

	// Create an equation systems object and declare its variables.
	ExplicitShellSystem & system = equation_systems.add_system<ExplicitShellSystem> ("Shell");
	ShellMeasurementSystem & meas_system = equation_systems.add_system<ShellMeasurementSystem>("ShellMeasurements");
    system.attach_init_function(random_initial_conditions);
	std::cout<<"Initializing equation system\n";
	equation_systems.init ();
    
    
	equation_systems.print_info();
    equation_systems.parameters.print();
	PerfLog perf_log("newshell");
	system.set_measurement_system(&meas_system);
    std::cout<<"Computing shell system mass  \n";
	system.compute_lumped_mass();
	std::cout<<"ShellSystem mass computed\n";
    system.set_initial_volume();
    double initial_volume = system.volume0;
    std::cout<<"Initial volume computed\n";
    system.mu_vol = penalty_mu_vol;
    
    double R0=scalef; //radius of initial hemisphere
    double naturalrad=R0; //start from initial hemisphere

    
    
 	for (int step=0; step<Nsteps; step++)
	{
	  //system.rebuild_voxel_list=!(step%30);

        if (!(step%50)) {
			std::cout<<"Solving for timestep "<<step<< "  (time="<<system.time<<")"<<std::endl;
            std::cout<<"Growth factors (kappabar,L0): "<<growth_factor_kappabar<<", "<<growth_factor_L0<<std::endl;
           
            std::cout<<"(Target red. vol., target vol., curr. vol, constr. energy): ("<<system.volume0/initial_volume<<", "
            <<system.volume0<<", "<<system.volume<<", "<< 0.5*system.mu_vol*(system.volume0-system.volume)*(system.volume0-system.volume)
            <<")"<<std::endl;
            
        }
		if (!(step%measure_each) || !(step%plot_each)) {
			system.measure_now = true;
		}
       // std::cout<<"Main 1\n"<<std::endl;
		system.newmark_predict();
//        std::cout<<"Main 2\n"<<std::endl;

      //  std::cout<<"predicted\n";
        system.enforce_ghost_bc(); //Apply clamped (1) or free(0) boundary conditions
      //  std::cout<<"bc enforced\n";
  //      std::cout<<"Main 3\n"<<std::endl;

        system.compute_residual();
    //    std::cout<<"Main 4\n"<<std::endl;

        //system.collision_detection();
		system.newmark_correct();
      //  std::cout<<"Main 5\n"<<std::endl;
        if (growth_rate_kappabar!=0 || growth_rate_L0!=0) {
            growth_factor_kappabar += growth_rate_kappabar*dt;
            growth_factor_L0 += growth_rate_L0*dt;
            if (growth_factor_kappabar>max_kappabar)
                growth_factor_kappabar=max_kappabar;
            if (growth_factor_kappabar<min_kappabar)
                growth_factor_kappabar=min_kappabar;
            if (growth_factor_L0>max_L0)
                growth_factor_L0=max_L0;
            if (growth_factor_L0<min_L0)
                growth_factor_L0=min_L0;
        }
            equation_systems.parameters.set<double>("growth_factor_kappabar") = growth_factor_kappabar;
            equation_systems.parameters.set<double>("growth_factor_L0") = growth_factor_L0;
        
        if (system.time>=reduced_volume_quench_start_time && reduced_volume_rate != 0.0) {
            system.volume0 += initial_volume*reduced_volume_rate*dt;
            if (system.volume0 > initial_volume*max_reduced_volume)
                system.volume0 = initial_volume*max_reduced_volume;
            if (system.volume0 < initial_volume*min_reduced_volume)
                system.volume0 = initial_volume*min_reduced_volume;
        }
        
		if (!(step%plot_each)) {
			OStringStream file_name;
			file_name<< data_dir;
			file_name << "result-";
			OSSRealzeroright(file_name,8,0, step);
			file_name<<".pvtu";
           // std::cout<<"preparing measurement sys\n";
			meas_system.prepare_measurement();
            //std::cout<<"prepared measurement sys\n";
			VTKIO (mesh).write_equation_systems (file_name.str(), equation_systems);
            
            //std::cout<<"EQS written\n";
            
        }
		
		if (system.measure_now ) {
            //std::cout<<"writing measurement\n";

			system.write_measurements(equation_systems.parameters.get<std::string>("measurement_filename"));
            //std::cout<<"wrote measurement\n";

			system.measure_now = false;
 		}

		system.update_time(system.time+system.deltat);
		
	}
	
	std::cout<<"System solved\n";
	OStringStream file_name;
	file_name<< data_dir;
	file_name<<"final-solution.pvtu";
	VTKIO (mesh).write_equation_systems(file_name.str(), equation_systems);
	//TODO: Dump solution to data_dir!
	file_name.str("");
	file_name<<data_dir;
	file_name<<"system_dump.xdr";
	//system.write(file_name.str(), libMeshEnums::WRITE);
	return 0;
}

