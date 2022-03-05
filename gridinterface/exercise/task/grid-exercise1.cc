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

int main(int argc, char** argv)
{
    try
    {
        const auto& helper = Dune::MPIHelper::instance(argc, argv);

        // typedef for grid of preferred type
        static const int dim = 2;
        typedef Dune::YaspGrid<dim> GridType;

        // define the extensions of the domain: a unit square with N by N cells
        Dune::FieldVector<double,dim> lowerleft(0.0);
        Dune::FieldVector<double,dim> upperright(1.0);
        auto N  = Dune::filledArray<dim, unsigned int>(4);

        // build a structured grid
        auto grid = Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerleft, upperright, N);

        // print the grid! This will generate a file 'printgrid_0.png' that shows some
        // detailed information about the grid (see exercise sheet).
        Dune::printGrid(*grid, helper);

        // Extract the grid view we would like to use
        // You could use a different grid view here -- the grid itself is not used anymore
        auto gridview = grid->leafGridView();

        // iterate over all entities of the grid
        for (const auto& e : elements(gridview))
        {
            auto geo = e.geometry();
            auto gt = geo.type();

            std::cout << std::endl;
            std::cout << "Visiting a cell with geometry type: " << gt << std::endl;
            // Calculate the center of the cell by summing up all corner vectors
            // and dividing by the number of corners.  Print the results.
            // Also print the center as returned by center() method of the geometry.
            // *** ADD CODE HERE (1/4) ***

            std::cout << std::endl;
            // iterate over intersections of current entity
            for (const auto& is : intersections(gridview, e))
            {
                auto fgeo = is.geometry();
                auto fgt = fgeo.type();
                auto fgeo_self = is.geometryInInside();

                std::cout << " Intersection endpoints: ";
                // Output the corner and center coordinates of the intersection
                // (the methods corner and center return global coordinates)
                // *** ADD CODE HERE (2/4) ***

                if (is.neighbor())       // intersection with a neighboring cell
                {
                    // The coordinate of the center of an intersection with respect to the intersection itself is 0.5.
                    // You can as well extract it from the reference element:
                    auto localcenter = referenceElement(fgeo).position(0,0);

                    // Furthermore, the very same point can be expressed
                    // 1.) in global coordinates with respect to the domain origin,
                    // 2.) in coordinates with respect to the current element or
                    // 3.) in coordinates with resepct to the neighboring element.
                    // *** ADD CODE HERE (3/4) ***
                }
                else if (is.boundary())  // intersection with the boundary
                {
                    std::cout << " Boundary intersection center = 0.5";
                    // get the global coordinates of the center of the intersection
                    // and the coordinates with respect to the owning element
                    // *** ADD CODE HERE (4/4) ***
                    std::cout << std::endl;
                }
            }
        }
    }
    // catch exceptions
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
