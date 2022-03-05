#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<math.h>
#include<iostream>
#include<vector>
#include<map>
#include<string>

#include<sys/stat.h>

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/typetraits.hh>
#include<dune/common/timer.hh>
#include <dune/common/float_cmp.hh> 	

#include<dune/grid/io/file/vtk.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/yaspgrid.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>

#include<dune/pdelab/finiteelementmap/qkdg.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/gridfunctionspace/subspace.hh>
#include<dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/p0.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/instationaryfilenamehelper.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/instationary/onestep.hh>
#include<dune/pdelab/function/callableadapter.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/instationary/onestep.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/finiteelement/l2orthonormal.hh>
#include<dune/pdelab/finiteelement/qkdglagrange.hh>
#include<dune/pdelab/solver/newton.hh>
#include<dune/pdelab/function/discretegridviewfunction.hh>
#include <dune/pdelab/localoperator/convectiondiffusiondg.hh>
#include <dune/pdelab/finiteelementmap/opbfem.hh>
#include <dune/pdelab/stationary/linearproblem.hh>

#include"timecapsule.hh"
#include"schemes.hh"
#include"navier-stokes-lop.hh"
#include"driver_coupled.hh"

#define STRUCTURED
#define SIMPLEX
#define DEGREE 3

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
    {
    if(helper.rank()==0)
      std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
    }

    // open ini file and parse it in
    Dune::ParameterTree ptree;
    const std::string config_filename("navier-stokes-rayleigh-benard.ini");
    try {
      Dune::ParameterTreeParser::readINITree(config_filename, ptree);
      Dune::ParameterTreeParser ptreeparser;
      ptreeparser.readOptions(argc,argv,ptree);
    }
    catch(...) {
      std::cerr << "The configuration file \"navier-stokes-rayleigh-benard.ini\" could not be read. "
        "Exiting..." << std::endl;
      exit(1);
    } 

    // construct grid
    std::string filename = ptree.get<std::string>("grid.meshfile");
    const int refinement = ptree.get<int>("grid.refinement");
#if HAVE_UG // see if we have UG
    typedef Dune::UGGrid<2> Grid;
#else 
    std::cout << "Example requires UG grid!" << std::endl;
#endif
#if HAVE_UG
#ifdef STRUCTURED
    Dune::StructuredGridFactory<Grid> factory;
    Dune::FieldVector<double,2> lowerLeft(0.0);
    auto upperRight = ptree.get<Dune::FieldVector<double,2> >("grid.extend");
    auto cells = ptree.get<std::array<unsigned int,2> >("grid.cells");
#ifdef CUBE
    auto gridp = factory.createCubeGrid(lowerLeft,upperRight,cells);
#else
    auto gridp = factory.createSimplexGrid(lowerLeft,upperRight,cells);
#endif
#else
    Dune::GridFactory<Grid> factory;
    Dune::GmshReader<Grid>::read(factory,filename,true,true);
    std::shared_ptr<Grid> gridp(factory.createGrid());
#endif
    
    // refine grid
    Dune::Timer timer;
    gridp->globalRefine(refinement);
    std::cout << "Time for mesh refinement " << timer.elapsed()
              << " seconds" << std::endl;

    // types & constants
    typedef Grid::LeafGridView GV;
    typedef double RF;
    const int dim = GV::dimension;

    // do computation on a grid view
    GV gv=gridp->leafGridView();

    // now define the specific problem by supplying
    // a) boundary condition type for a given point at the boundary
    // b) information about matching(!) constraints for individual solution components
    // c) boundary values and initial condition function
    
    // make a scalar boundary condition type function for the whole system
    const RF domainX = 1.0; // extend of domain in X
    const RF domainY = 1.0; // extend of domain in Y
    const RF eps=1e-7;       // accuracy for boundary detection
    auto bctypelambda = [&](const auto& x){ // Dirichlet for x-component of velocity
      return Dune::PDELab::NavierStokesBoundaryCondition::noslip; // means Dirichlet
    };

    // define matching constraints on solution components
    // we assume that type of b.c. does not depend on time
    auto buxlambda = [&](const auto& x){ // Dirichlet for x-component of velocity
      return true;
    };
    auto bux = Dune::PDELab::makeBoundaryConditionFromCallable(gv,buxlambda);
    auto buylambda = [&](const auto& x){ // Dirichlet for y-component of velocity
      return true;
    };
    auto buy = Dune::PDELab::makeBoundaryConditionFromCallable(gv,buylambda);
    auto bu = Dune::PDELab::CompositeConstraintsParameters<decltype(bux),decltype(buy)>(bux,buy);
    auto bplambda = [](const auto& x){return false;};
    auto bp = Dune::PDELab::makeBoundaryConditionFromCallable(gv,bplambda);
    auto bconstraints = Dune::PDELab::CompositeConstraintsParameters<decltype(bu),decltype(bp)>(bu,bp);

    // time capsule that allows to have time-dependent grid function
    TimeCapsule tc(0.0);
    
    // make combined Dirichlet and initial value function
    auto gulambda = [&](const auto& x){
      Dune::FieldVector<RF,dim> u(0.0);
      return u;
    };
    auto gu = Dune::PDELab::makeInstationaryGridFunctionFromCallable(gv,gulambda,tc);
    auto rho_0 = ptree.get<RF>("problem.rho_0"); // rho = rho_0 - alpha*temperature
    auto gplambda = [&](const auto& x){
      RF s = rho_0*(0.5*domainY-x[1]);
      return s;
    };
    auto gp = Dune::PDELab::makeInstationaryGridFunctionFromCallable(gv,gplambda,tc);
    auto g = Dune::PDELab::CompositeGridFunction<decltype(gu),decltype(gp)>(gu,gp);


    // Parameters for the transport problem

    // first we have the boundary condition type
    auto Tbctypelambda = [&](const auto& x){
      if (x[0]<eps || x[0]>domainX-eps) return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann;
      return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet;
    };

    // ... and combined Dirichlet and initial value function for temperature distribution
    TimeCapsule tc2(0.0); // another time capsule for the temperature
    auto width = 0.05*domainX;
    auto Tglambda = [&](const auto& x){
      RF T = 0.0;
      RF ramp;
      if (x[0]>=8.0/21.0 && x[0]<=9.0/21.0)
        ramp = 0.025*domainY * (1.0 + std::sin( x[0]*21.0*M_PI ));
      else
        ramp = 0.025*domainY;
      if (x[1]>=ramp) return T;
      T = (1.0-x[1]/ramp)*(1.0-x[1]/ramp);
      return T;
    };
    auto Tg = Dune::PDELab::makeInstationaryGridFunctionFromCallable(gv,Tglambda,tc2);

    // call the general driver
#ifdef CUBE 
    auto scheme = TaylorHood_21_Quadrilateral(gv);
#else
#if (DEGREE==2)
    auto scheme = TaylorHood_21_Triangle(gv);
#endif
#if (DEGREE==3)
    auto scheme = TaylorHood_32_Triangle(gv);
#endif
#endif

    // call driver
    driver_coupled(gv,scheme,bctypelambda,bconstraints,g,Tbctypelambda,Tg,ptree);
    
#endif

   }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  return 1;
  }
}
