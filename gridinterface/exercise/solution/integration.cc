#include "config.h"

#include <array>
#include <iostream>

#include <dune/common/filledarray.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/printgrid.hh>
#include <dune/pdelab/common/quadraturerules.hh>


// This is the function f to integrate
double f(Dune::FieldVector<double, 2> x)
{
  return std::exp(x[0])*std::exp(x[1]);
}

int main(int argc, char** argv)
{
    try
    {
        Dune::MPIHelper::instance(argc, argv);

        // typedef for grid of preferred type
        static const int dim = 2;
        typedef Dune::YaspGrid<dim> GridType;

        // define the extensions of the domain: a unit square with N by N cells
        Dune::FieldVector<double,dim> lowerleft(0.0);
        Dune::FieldVector<double,dim> upperright(1.0);
        auto N = Dune::filledArray<dim, unsigned int>(4);

        // build a structured grid
        auto grid = Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerleft, upperright, N);

        // Extract the grid view we would like to use
        // You could use a different grid view here -- the grid itself is not used anymore
        auto gridview = grid->leafGridView();
        double integral = 0.0;

        // iterate over all entities of the grid
        for (const auto& e : elements(gridview))
        {
            auto geo = e.geometry();

            auto rule = Dune :: PDELab :: quadratureRule (geo, 3);
            for ( const auto & qp : rule )
            {
                integral += qp.weight()*f(geo.global(qp.position()))*geo.integrationElement(qp.position());
            }
        }
        std::cout << "Computed the integral: " << integral << std::endl;
    }
    // catch exceptions
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
