// -*- tab-width: 4; indent-tabs-mode: nil -*-
// Beware of line number changes, they may corrupt docu!
/** \file
    \brief Solve Poisson equation with P1 conforming finite elements
*/
// always include the config file
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
// C++ includes
#include<math.h>
#include<iostream>
// dune-common includes
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/timer.hh>
// dune-geometry includes
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
// dune-grid includes
#include<dune/grid/onedgrid.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#include<dune/alugrid/dgf.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
// dune-istl included by pdelab
// dune-pdelab includes
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/function/callableadapter.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/variablefactories.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/solver/newton.hh>

// include all components making up this tutorial
#include"nonlinearpoissonfv.hh"
#include"problem.hh"
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
    ptreeparser.readINITree("tutorial02.ini",ptree);
    ptreeparser.readOptions(argc,argv,ptree);

    // read ini file
    const int dim = ptree.get<int>("grid.dim");
    const int refinement = ptree.get<int>("grid.refinement");

    // in 1d use OneDGrid
    if (dim==1)
      {
        // read grid parameters from input file
        typedef Dune::OneDGrid::ctype DF;
        const DF a = ptree.get<DF>("grid.oned.a");
        const DF b = ptree.get<DF>("grid.oned.b");
        const unsigned int N = ptree.get<int>("grid.oned.elements");

        // create equidistant intervals
        typedef std::vector<DF> Intervals;
        Intervals intervals(N+1);
        for(unsigned int i=0; i<N+1; ++i)
          intervals[i] = a + DF(i)*(b-a)/DF(N);

        // Construct grid
        typedef Dune::OneDGrid Grid;
        Grid grid(intervals);
        grid.globalRefine(refinement);

        // call generic function
        typedef Dune::OneDGrid::LeafGridView GV;
        GV gv=grid.leafGridView();
        driver(gv,ptree);
      }

    // YaspGrid section
    if (dim==2)
      {
        const int dim=2;
        typedef Dune::YaspGrid<dim> Grid;
        typedef Grid::ctype DF;
        Dune::FieldVector<DF,dim> L;
        L[0] = ptree.get("grid.structured.LX",(double)1.0);
        L[1] = ptree.get("grid.structured.LY",(double)1.0);
        std::array<int,dim> N;
        N[0] = ptree.get("grid.structured.NX",(int)10);
        N[1] = ptree.get("grid.structured.NY",(int)10);
        std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(new Grid(L,N));
        gridp->globalRefine(refinement);
        typedef Grid::LeafGridView GV;
        GV gv=gridp->leafGridView();
        driver(gv,ptree);
      }
    if (dim==3)
      {
        const int dim=3;
        typedef Dune::YaspGrid<dim> Grid;
        typedef Grid::ctype DF;
        Dune::FieldVector<DF,dim> L;
        L[0] = ptree.get("grid.structured.LX",(double)1.0);
        L[1] = ptree.get("grid.structured.LY",(double)1.0);
        L[2] = ptree.get("grid.structured.LZ",(double)1.0);
        std::array<int,dim> N;
        N[0] = ptree.get("grid.structured.NX",(int)1);
        N[1] = ptree.get("grid.structured.NY",(int)1);
        N[2] = ptree.get("grid.structured.NZ",(int)1);
        std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(new Grid(L,N));
        gridp->globalRefine(refinement);
        typedef Grid::LeafGridView GV;
        GV gv=gridp->leafGridView();
        driver(gv,ptree);
      }
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
