
#include "config.h"
#include "dune/common/parallel/mpihelper.hh"
#include "dune/common/parametertree.hh"
#include "dune/common/parametertreeparser.hh"
#include "dune/pdelab/backend/istl.hh"
#include "dune/grid/yaspgrid.hh"
#include "dune/pdelab/finiteelementmap/qkfem.hh"
#include "dune/pdelab/constraints/conforming.hh"
#include "dune/testtools/gridconstruction.hh"
#include "heatequation_poisson_operator.hh"
#include "dune/pdelab/function/callableadapter.hh"
#include "dune/pdelab/gridoperator/onestep.hh"
#include "dune/pdelab/gridoperator/gridoperator.hh"
#include "heatequation_mass_operator.hh"
#include "dune/pdelab/newton/newton.hh"
#include "dune/pdelab/gridfunctionspace/vtk.hh"
#include "dune/grid/io/file/vtk/vtksequencewriter.hh"
#include "dune/grid/io/file/vtk/subsamplingvtkwriter.hh"
#include "string"
#include "dune/codegen/vtkpredicate.hh"
#include <random>



int main(int argc, char** argv){  
  try
  {    

    if (argc != 2){
      std::cerr << "This program needs to be called with an ini file" << std::endl;
      return 1;
    }
    
    // Initialize basic stuff...    
    Dune::MPIHelper& mpihelper = Dune::MPIHelper::instance(argc, argv);
    Dune::ParameterTree initree;
    Dune::ParameterTreeParser::readINITree(argv[1], initree);
    using RangeType = double;
    double time = 0.0;
    
    // Setup grid (view)...    
    using Grid = Dune::YaspGrid<2, Dune::EquidistantCoordinates<RangeType, 2>>;
    using GV = Grid::LeafGridView;
    using DF = Grid::ctype;
    IniGridFactory<Grid> factory(initree);
    std::shared_ptr<Grid> grid = factory.getGrid();
    GV gv = grid->leafGridView();
    
    // Set up finite element maps...    
    using CG1_FEM = Dune::PDELab::QkLocalFiniteElementMap<GV, DF, RangeType, 1>;
    CG1_FEM cg1_fem(gv);
    
    // Set up grid function spaces...    
    using VectorBackendCG1 = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
    using DirichletConstraintsAssember = Dune::PDELab::ConformingDirichletConstraints;
    using CG1_dirichlet_GFS = Dune::PDELab::GridFunctionSpace<GV, CG1_FEM, DirichletConstraintsAssember, VectorBackendCG1>;
    CG1_dirichlet_GFS cg1_dirichlet_gfs_(gv, cg1_fem);
    cg1_dirichlet_gfs_.name("cg1_dirichlet_gfs_");
    
    // Set up constraints container...    
    auto bctype_0000_lambda = [&](const auto& x){ return 1.0; };
    auto bctype_0000 = Dune::PDELab::makeBoundaryConditionFromCallable(gv, bctype_0000_lambda);
    using CG1_dirichlet_GFS_CC = CG1_dirichlet_GFS::ConstraintsContainer<RangeType>::Type;
    CG1_dirichlet_GFS_CC cg1_dirichlet_gfs__cc;
    cg1_dirichlet_gfs__cc.clear();
    auto bctype_0001_lambda = [&](const auto& x){ return 1.0; };
    auto bctype_0001 = Dune::PDELab::makeBoundaryConditionFromCallable(gv, bctype_0001_lambda);
    Dune::PDELab::constraints(bctype_0001, cg1_dirichlet_gfs_, cg1_dirichlet_gfs__cc);
    
    // Set up grid grid operators...    
    using LOP_POISSON = PoissonOperator<CG1_dirichlet_GFS, CG1_dirichlet_GFS>;
    LOP_POISSON lop_poisson(cg1_dirichlet_gfs_, cg1_dirichlet_gfs_, initree);
    lop_poisson.setTime(0.0);
    using MatrixBackend = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    using GO_poisson = Dune::PDELab::GridOperator<CG1_dirichlet_GFS, CG1_dirichlet_GFS, LOP_POISSON, MatrixBackend, DF, RangeType, RangeType, CG1_dirichlet_GFS_CC, CG1_dirichlet_GFS_CC>;
    using LOP_MASS = MassOperator<CG1_dirichlet_GFS, CG1_dirichlet_GFS>;
    using GO_mass = Dune::PDELab::GridOperator<CG1_dirichlet_GFS, CG1_dirichlet_GFS, LOP_MASS, MatrixBackend, DF, RangeType, RangeType, CG1_dirichlet_GFS_CC, CG1_dirichlet_GFS_CC>;
    using IGO = Dune::PDELab::OneStepGridOperator<GO_poisson,GO_mass>;
    cg1_dirichlet_gfs_.update();
    int generic_dof_estimate =  4 * cg1_dirichlet_gfs_.maxLocalSize();
    int dofestimate = initree.get<int>("istl.number_of_nnz", generic_dof_estimate);
    MatrixBackend mb(dofestimate);
    GO_poisson go_poisson(cg1_dirichlet_gfs_, cg1_dirichlet_gfs__cc, cg1_dirichlet_gfs_, cg1_dirichlet_gfs__cc, lop_poisson, mb);
    std::cout << "gfs with " << cg1_dirichlet_gfs_.size() << " dofs generated  "<< std::endl;
    std::cout << "cc with " << cg1_dirichlet_gfs__cc.size() << " dofs generated  "<< std::endl;
    LOP_MASS lop_mass(cg1_dirichlet_gfs_, cg1_dirichlet_gfs_, initree);
    GO_mass go_mass(cg1_dirichlet_gfs_, cg1_dirichlet_gfs__cc, cg1_dirichlet_gfs_, cg1_dirichlet_gfs__cc, lop_mass, mb);
    std::cout << "gfs with " << cg1_dirichlet_gfs_.size() << " dofs generated  "<< std::endl;
    std::cout << "cc with " << cg1_dirichlet_gfs__cc.size() << " dofs generated  "<< std::endl;
    IGO igo(go_poisson, go_mass);
    
    // Set up solution vectors...    
    using V_POISSON = Dune::PDELab::Backend::Vector<CG1_dirichlet_GFS,DF>;
    V_POISSON x_poisson(cg1_dirichlet_gfs_);
    x_poisson = 0.0;
    auto lambda_0000 = [&](const auto& is, const auto& xl){ auto x=is.geometry().global(xl); return cos(6.283185307179586 * lop_poisson.getTime()) * cos(3.141592653589793 * x[0]) * cos(3.141592653589793 * x[0]) * cos(3.141592653589793 * x[1]) * cos(3.141592653589793 * x[1]);; };
    auto func_0000 = Dune::PDELab::makeInstationaryGridFunctionFromCallable(gv, lambda_0000, lop_poisson);
    Dune::PDELab::interpolate(func_0000, cg1_dirichlet_gfs_, x_poisson);
    
    // Set up (non)linear solvers...    
    using LinearSolver = Dune::PDELab::ISTLBackend_SEQ_SuperLU;
    LinearSolver ls(false);
    using SNP = Dune::PDELab::NewtonMethod<IGO, LinearSolver>;
    SNP snp(igo, ls);

    // Do visualization...    
    using VTKSW = Dune::VTKSequenceWriter<GV>;
    using VTKWriter = Dune::SubsamplingVTKWriter<GV>;
    Dune::RefinementIntervals subint(initree.get<int>("vtk.subsamplinglevel", 1));
    VTKWriter vtkwriter(gv, subint);
    std::string vtkfile = initree.get<std::string>("wrapper.vtkcompare.name", "output");
    VTKSW vtkSequenceWriter(std::make_shared<VTKWriter>(vtkwriter), vtkfile);
    CuttingPredicate predicate;
    Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter, cg1_dirichlet_gfs_, x_poisson, Dune::PDELab::vtk::defaultNameScheme(), predicate);
    vtkSequenceWriter.write(time, Dune::VTK::appendedraw);
    
    // Set up instationary stuff...    
    using OSM = Dune::PDELab::OneStepMethod<RangeType, IGO, SNP, V_POISSON, V_POISSON>;
    using TSM = Dune::PDELab::OneStepThetaParameter<RangeType>;
    TSM tsm(initree.get<double>("instat.theta",1.0));
    OSM osm(tsm,igo,snp);
    
    double T = initree.get<double>("instat.T", 1.0);
    double dt = initree.get<double>("instat.dt", 0.1);
    int step_number(0);int output_every_nth = initree.get<int>("instat.output_every_nth", 1);
    while (time<T-1e-8){
      // Assemble constraints for new time step
  lop_poisson.setTime(time+dt);
  Dune::PDELab::constraints(bctype_0000, cg1_dirichlet_gfs_, cg1_dirichlet_gfs__cc);


      // Do time step
      V_POISSON x_poissonnew(x_poisson);
      osm.apply(time, dt, x_poisson, func_0000, x_poissonnew);
    
      // Accept new time step
      x_poisson = x_poissonnew;
      time += dt;
    
      step_number += 1;
      if (step_number%output_every_nth == 0){
        // Output to VTK File
        vtkSequenceWriter.vtkWriter()->clear();
        Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter, cg1_dirichlet_gfs_, x_poisson,
                                             Dune::PDELab::vtk::defaultNameScheme(), predicate);
        vtkSequenceWriter.write(time, Dune::VTK::appendedraw);
      }
    }
    
    
    // Maybe print residuals and matrices to stdout...    
    if (initree.get<bool>("printresidual", false)) {
      using Dune::PDELab::Backend::native;
      V_POISSON r(x_poisson);
      // Setup random input
      std::size_t seed = 0;
      auto rng = std::mt19937_64(seed);
      auto dist = std::uniform_real_distribution<>(-1., 1.);
      for (auto& v : x_poisson)
        v = dist(rng);
      r=0.0;
      go_poisson.residual(x_poisson, r);
      Dune::printvector(std::cout, native(r), "residual vector", "row");
    }
    if (initree.get<bool>("printmatrix", false)) {
      using Dune::PDELab::Backend::native;
      V_POISSON r(x_poisson);
      // Setup random input
      std::size_t seed = 0;
      auto rng = std::mt19937_64(seed);
      auto dist = std::uniform_real_distribution<>(-1., 1.);
      for (auto& v : x_poisson)
        v = dist(rng);
      using M = typename GO_poisson::Traits::Jacobian;
      M m(go_poisson);
      go_poisson.jacobian(x_poisson,m);
      using Dune::PDELab::Backend::native;
      Dune::printmatrix(std::cout, native(m),"global stiffness matrix","row",9,1);
    }
    
    // Maybe calculate errors for test results...    
    bool testfail(false);
    
    // Return statement...    
    return testfail;
    
  }  
  catch (Dune::Exception& e)
  {    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }  
  catch (std::exception& e)
  {    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }  
}

