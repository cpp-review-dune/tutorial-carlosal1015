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
#include<dune/grid/io/file/vtk/vtkwriter.hh>
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
#include<dune/pdelab/finiteelementmap/pkfem.hh>
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
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

// include all components making up this tutorial
#include"poissonp1.hh"
#include"driver.hh"

//#undef HAVE_DUNE_ALUGRID
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

    // Read parameters from ini file
    Dune::ParameterTree ptree;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree("tutorial00.ini",ptree);
    ptreeparser.readOptions(argc,argv,ptree);

    // make grid
    const int dim = ptree.get("grid.dim",(int)2);
    const int refinement = ptree.get<int>("grid.refinement");
    if (dim==1)
      {
        // read grid parameters from input file
        typedef Dune::OneDGrid::ctype ctype;
        const ctype a = ptree.get<ctype>("grid.oned.a");
        const ctype b = ptree.get<ctype>("grid.oned.b");
        const unsigned int N = ptree.get<int>("grid.oned.elements");

        // create equidistant intervals
        typedef std::vector<ctype> Intervals;
        Intervals intervals(N+1);
        for(unsigned int i=0; i<N+1; ++i)
          intervals[i] = a + ctype(i)*(b-a)/ctype(N);

        // Construct grid
        typedef Dune::OneDGrid Grid;
        Grid grid(intervals);
        grid.globalRefine(refinement);

        // call generic function
        driver(grid.leafGridView(),ptree);
      }
    if (dim==2)
      {
#if HAVE_DUNE_ALUGRID
  typedef Dune::ALUGrid<2,2,Dune::simplex,
                        Dune::nonconforming> Grid;
#elif HAVE_UG
  typedef Dune::UGGrid<2> Grid;
#else  // ! (HAVE_UG || HAVE_DUNE_ALUGRID)
  std::cout << "Example requires a simplex grid!" << std::endl;
#endif
#if (HAVE_UG || HAVE_DUNE_ALUGRID)
  std::string filename = ptree.get("grid.twod.filename",
                                   "unitsquare.msh");
  Dune::GridFactory<Grid> factory;
  Dune::GmshReader<Grid>::read(factory,filename,true,true);
  std::shared_ptr<Grid> gridp(factory.createGrid());
  Dune::Timer timer;
  gridp->globalRefine(refinement);
  std::cout << "Time for mesh refinement " << timer.elapsed()
            << " seconds" << std::endl;
  driver(gridp->leafGridView(),ptree);
#endif
      }
    if (dim==3)
      {
#if HAVE_DUNE_ALUGRID
        typedef Dune::ALUGrid<3,3,Dune::simplex,Dune::nonconforming> Grid;
#elif HAVE_UG
        typedef Dune::UGGrid<3> Grid;
#else  // ! (HAVE_UG || HAVE_DUNE_ALUGRID)
        std::cout << "This example requires a simplex grid!" << std::endl;
#endif
#if (HAVE_UG || HAVE_DUNE_ALUGRID)
        std::string filename = ptree.get("grid.threed.filename","unitcube.msh");
        Dune::GridFactory<Grid> factory;
        Dune::GmshReader<Grid>::read(factory,filename,true,true);
        std::shared_ptr<Grid> gridp(factory.createGrid());
        Dune::Timer timer;
        gridp->globalRefine(refinement);
        std::cout << "Time for mesh refinement " << timer.elapsed() << " seconds" << std::endl;
        driver(gridp->leafGridView(),ptree);
#endif
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
