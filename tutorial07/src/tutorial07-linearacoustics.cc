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
#include<dune/grid/onedgrid.hh>
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
#include<dune/pdelab/gridfunctionspace/subspace.hh>
#include<dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
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

#include"hyperbolicdg.hh"

//===============================================================
// Include your hyperbolic model and problem to solve
//===============================================================
/// tex: include
//Linear Acoustics
#include"linearacoustics.hh"//model
#include"linearacousticsproblem.hh"
#include"numericalflux.hh"
/// tex: include

#include"driver.hh"

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

    // open ini file
    Dune::ParameterTree ptree;
    Dune::ParameterTreeParser ptreeparser;

    ptreeparser.readINITree("tutorial07-linearacoustics.ini",ptree);

    ptreeparser.readOptions(argc,argv,ptree);

    // read ini file
    const int dim = ptree.get<int>("grid.dim");
    const int refinement = ptree.get<int>("grid.refinement",(int)0);
    const int degree = ptree.get<int>("fem.degree");

    // parallel overlapping yaspgrid version
    if (ptree["grid.manager"] == "yasp")
      {
        if (dim == 1)
          {
            constexpr int dim = 1;
            // read grid parameters from input file
            using DF = Dune::OneDGrid::ctype;
            auto a = 0;
            std::array<double,dim> lower_left; for (int i=0; i<dim; i++) lower_left[i]=0.0;
            auto upper_right = ptree.get<Dune::FieldVector<double,dim> >("grid.L");
            auto cells = ptree.get<std::array<int,dim> >("grid.N");

            auto b = upper_right[0];
            auto N = cells[0];
            // OneD Grid
            // create equidistant intervals
            using Intervals = std::vector<DF>;
            Intervals intervals(N+1);
            for(int i=0; i<N+1; ++i)
              // store start of every interval and the end of grid
              intervals[i] = a + DF(i)*(b-a)/DF(N);

            // Construct grid
            using Grid = Dune::OneDGrid;
            Grid grid(intervals);
            grid.globalRefine(refinement);

            // call generic function
            using GV = Dune::OneDGrid::LeafGridView;
            GV gv = grid.leafGridView();

            /// tex: promodflux
            //create problem (setting)
            using PROBLEM = Problem<GV,GV::Grid::ctype>;
            PROBLEM problem;

            //create model on a given setting
            using MODEL = Model<PROBLEM >;
            MODEL model(problem);

            //create numerical flux
            using NUMFLUX = VariableFluxVectorSplitting<MODEL>;
            NUMFLUX numflux(model);
            /// tex: promodflux

            if (degree==0)
              {
                using FEM = Dune::PDELab::QkDGLocalFiniteElementMap<GV::Grid::ctype,double,0,dim,
                                                                    Dune::PDELab::QkDGBasisPolynomial::legendre>;
                FEM fem;
                driver<GV,FEM, NUMFLUX>(gv,fem,numflux,ptree);
              }
            if (degree==1)
              {
                using FEM = Dune::PDELab::QkDGLocalFiniteElementMap<GV::Grid::ctype,double,1,dim,
                                                                    Dune::PDELab::QkDGBasisPolynomial::legendre>;
                FEM fem;
                driver<GV,FEM, NUMFLUX>(gv,fem,numflux,ptree);
              }
            if (degree==2)
              {
                /// tex: fem
                using FEM = Dune::PDELab::QkDGLocalFiniteElementMap
                  <GV::Grid::ctype, double, 2, dim,
                   Dune::PDELab::QkDGBasisPolynomial::legendre>;
                FEM fem;
                driver<GV,FEM, NUMFLUX>(gv,fem,numflux,ptree);
                /// tex: fem
              }
            return 0;

          }


        if (dim == 2)
          {
            const int dim=2;

            std::bitset<dim> periodic(false);
            int overlap=1;
            std::array<double,dim> lower_left; for (int i=0; i<dim; i++) lower_left[i]=0.0;
            auto upper_right = ptree.get<Dune::FieldVector<double,dim> >("grid.L");
            auto cells = ptree.get<std::array<int,dim> >("grid.N");

            // make grid
            using GM = Dune::YaspGrid<dim>;
            GM grid(upper_right,cells,periodic,overlap,helper.getCommunicator());
            grid.refineOptions(false); // keep overlap in cells
            grid.globalRefine(refinement);
            // grid view
            using GV = GM::LeafGridView ;
            GV gv=grid.leafGridView();

            //create problem (setting)
            using PROBLEM = Problem<GV,GV::Grid::ctype>;
            PROBLEM problem;

            //create model on a given setting
            using MODEL = Model<PROBLEM>;
            MODEL model(problem);

            //create numerical flux
            using NUMFLUX = VariableFluxVectorSplitting<MODEL>;
            NUMFLUX numflux(model);

            if (degree==0)
              {
                using FEM = Dune::PDELab::QkDGLocalFiniteElementMap<GV::Grid::ctype,double,0,dim,
                                                                    Dune::PDELab::QkDGBasisPolynomial::legendre>;
                FEM fem;
                driver<GV,FEM, NUMFLUX>(gv,fem,numflux,ptree);
              }
            if (degree==1)
              {
                using FEM = Dune::PDELab::QkDGLocalFiniteElementMap<GV::Grid::ctype,double,1,dim,
                                                                    Dune::PDELab::QkDGBasisPolynomial::legendre>;
                FEM fem;
                driver<GV,FEM, NUMFLUX>(gv,fem,numflux,ptree);
              }
            if (degree==2)
              {
                using FEM = Dune::PDELab::QkDGLocalFiniteElementMap<GV::Grid::ctype,double,2,dim,
                                                                    Dune::PDELab::QkDGBasisPolynomial::legendre>;
                FEM fem;
                driver<GV,FEM, NUMFLUX>(gv,fem,numflux,ptree);
              }
            return 0;
          }
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
