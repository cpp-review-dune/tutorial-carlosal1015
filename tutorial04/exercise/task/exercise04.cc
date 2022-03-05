// -*- tab-width: 4; indent-tabs-mode: nil -*-
// Beware of line number changes, they may corrupt docu!
/** \file
    \brief Solve Poisson equation with P1 conforming finite elements
*/
// always include the config file
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// C includes
#include<sys/stat.h>
// C++ includes
#include<cmath>
#include<iostream>
// dune-common includes
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/timer.hh>
// dune-grid includes
#include<dune/grid/onedgrid.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/grid/io/file/vtk.hh>
#include<dune/grid/io/file/gmshreader.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
// dune-istl included by pdelab
// dune-pdelab includes
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/functionutilities.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/function/callableadapter.hh>
#include<dune/pdelab/function/sqr.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/subspace.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/variablefactories.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/instationary/onestep.hh>
#include<dune/pdelab/newton/newton.hh>

// include all components making up this tutorial
#include"wavefem.hh"
#include"driver.hh"

//===============================================================
// Main program with grid setup
//===============================================================
int main(int argc, char** argv)
{
  try{
    // Maybe initialize Mpi
    Dune::MPIHelper&
      helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
      std::cout << "Parallel code run on "
                << helper.size() << " process(es)" << std::endl;

    // open ini file
    Dune::ParameterTree ptree;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree("exercise04.ini",ptree);
    ptreeparser.readOptions(argc,argv,ptree);

    // read ini file
    const int refinement = ptree.get<int>("grid.refinement");
    const int degree = ptree.get<int>("fem.degree");

#if GRIDDIM == 1
    {
      // read grid parameters from input file
      using DF = Dune::OneDGrid::ctype;
      auto a = ptree.get<DF>("grid.oned.a");
      auto b = ptree.get<DF>("grid.oned.b");
      auto N = ptree.get<unsigned int>("grid.oned.elements");

      // create equidistant intervals
      using Intervals = std::vector<DF>;
      Intervals intervals(N+1);
      for(unsigned int i=0; i<N+1; ++i)
        intervals[i] = a + DF(i)*(b-a)/DF(N);

      // Construct grid
      using Grid = Dune::OneDGrid;
      Grid grid(intervals);
      grid.globalRefine(refinement);

      // call generic function
      using GV = Dune::OneDGrid::LeafGridView;
      GV gv = grid.leafGridView();
      if (degree==1) {
        using FEM = Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,1>;
        FEM fem(gv);
        driver(gv,fem,ptree);
      }
      if (degree==2) {
        using FEM = Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,2>;
        FEM fem(gv);
        driver(gv,fem,ptree);
      }
      //add degree 3 for higher order exercise
    }
#endif

#if GRIDDIM == 2
    // UGGrid section
    {
      const int dim=2;

      using Grid = Dune::UGGrid<dim>;
      using DF = Grid::ctype;
      Dune::FieldVector<DF,dim> L;
      //upper right
      L[0] = ptree.get("grid.structured.LX",(double)1.0);
      L[1] = ptree.get("grid.structured.LY",(double)1.0);
      std::array<unsigned int,dim> N;
      N[0] = ptree.get("grid.structured.NX",(unsigned int)10);
      N[1] = ptree.get("grid.structured.NY",(unsigned int)10);
      //lower left
      Dune::FieldVector<double,dim> lowerleft(0.0);

      // build a structured simplex grid
      auto gridp = Dune::StructuredGridFactory<Grid>::createSimplexGrid(lowerleft, L, N);

      gridp->globalRefine(refinement);
      using GV = Grid::LeafGridView;
      GV gv=gridp->leafGridView();
      if (degree==1) {
        using FEM = Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,1>;
        FEM fem(gv);
        driver(gv,fem,ptree);
      }
      if (degree==2) {
        using FEM = Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,2>;
        FEM fem(gv);
        driver(gv,fem,ptree);
      }
      //add degree 3 for higher order exercise
    }
#endif

#if GRIDDIM == 3
    {
      const int dim=3;

      using Grid = Dune::UGGrid<dim>;
      using DF = Grid::ctype;
      Dune::FieldVector<DF,dim> L;
      //upper right
      L[0] = ptree.get("grid.structured.LX",(double)1.0);
      L[1] = ptree.get("grid.structured.LY",(double)1.0);
      L[2] = ptree.get("grid.structured.LZ",(double)1.0);
      std::array<unsigned int,dim> N;
      N[0] = ptree.get("grid.structured.NX",(unsigned int)1);
      N[1] = ptree.get("grid.structured.NY",(unsigned int)1);
      N[2] = ptree.get("grid.structured.NZ",(unsigned int)1);
      //lower left
      Dune::FieldVector<double,dim> lowerleft(0.0);

      // build a structured simplex grid
      auto gridp = Dune::StructuredGridFactory<Grid>::createSimplexGrid(lowerleft, L, N);

      gridp->globalRefine(refinement);
      using GV = Grid::LeafGridView;
      GV gv=gridp->leafGridView();
      if (degree==1) {
        using FEM = Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,1>;
        FEM fem(gv);
        driver(gv,fem,ptree);
      }
      if (degree==2) {
        using FEM = Dune::PDELab::PkLocalFiniteElementMap<GV,DF,double,2>;
        FEM fem(gv);
        driver(gv,fem,ptree);
      }
    }
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
