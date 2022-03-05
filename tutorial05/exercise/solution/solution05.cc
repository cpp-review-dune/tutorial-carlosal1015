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
#include<dune/pdelab/finiteelementmap/pkfem.hh>
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
#include<dune/pdelab/adaptivity/adaptivity.hh>

// include files needed for the exercise
#include<dune/pdelab/common/functionutilities.hh>
#include<dune/pdelab/function/sqr.hh>

// include all components making up this tutorial
#include"../../../tutorial01/src/nonlinearpoissonfem.hh"
#include"nonlinearpoissonfemestimator.hh"
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
    ptreeparser.readINITree("tutorial05.ini",ptree);
    ptreeparser.readOptions(argc,argv,ptree);

    // read ini file
    const int dim = 2;
    std::string gridmanager = "ug";
    // const int dim = ptree.get<int>("grid.dim");
    // std::string gridmanager = ptree.get("grid.manager","ug");
    const int degree = ptree.get<int>("fem.degree");

    // use UG grid if available and selected
    if (dim==2 && gridmanager=="ug")
      {
#if HAVE_UG
        typedef Dune::UGGrid<2> Grid;
        std::string filename = ptree.get("grid.twod.filename",
                                         "ldomain.msh");
        Dune::GridFactory<Grid> factory;
        Dune::GmshReader<Grid>::read(factory,filename,true,true);
        std::shared_ptr<Grid> gridp(factory.createGrid());
        if (degree==1) driver<Grid,1>(*gridp,ptree);
        if (degree==2) driver<Grid,2>(*gridp,ptree);
        if (degree==3) driver<Grid,3>(*gridp,ptree);
        if (degree==4) driver<Grid,4>(*gridp,ptree);
#else
        std::cout << "You selected ug as grid manager "
          << "but ug was not found during installation"
          << std::endl;
#endif
      }

    if (!(gridmanager=="ug") && !(gridmanager=="alu"))
      std::cout << "Example requires an unstructured grid!" << std::endl;

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
