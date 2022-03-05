#include "config.h"

#include <array>
#include <iostream>
#include <string>

#include <dune/common/filledarray.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/onedgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/utility/tensorgridfactory.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/printgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

int main(int argc, char** argv)
{
    try
    {
        // for parallel grids, MPI needs to be correctly initialized!
        Dune::MPIHelper::instance(argc, argv);

        // define a static variable for the grid dimension
        static const int dim = 2;

        // This is a list of the grid types we want to explore in this exercise.
        // Exactly one of those lines should not be a comment.
        //
        // Note: You can build an OneDGrid by adjusting dim above and
        // using the StructuredGridFactory.
        typedef Dune::YaspGrid<dim> GridType;
        //typedef Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<double, dim> > GridType;
        //typedef Dune::YaspGrid<dim, Dune::TensorProductCoordinates<double, dim> > GridType;
        //typedef Dune::OneDGrid GridType;
        //typedef Dune::UGGrid<dim> GridType;
        //typedef Dune::ALUGrid<dim, dim, Dune::simplex,Dune::nonconforming> GridType;

        // Build an equidistant grid:
        // Comment out if you do not want to build an equidistant grid.
        // Define the extensions of the domain: a unit square with N by N cells.
        Dune::FieldVector<double,dim> lowerleft(0.0);
        Dune::FieldVector<double,dim> upperright(1.0);
        auto N = Dune::filledArray<dim, unsigned int>(4);
        // build the equidistant grid
        auto grid = Dune::StructuredGridFactory<GridType>::createCubeGrid(lowerleft, upperright, N);

        // Build an unstructured grid:
        // Remove comments if you want to load a grid from a gmsh file.
//        std::string mshfile = "ldomain.msh";
//        auto grid = Dune::GmshReader<GridType>::read(mshfile);

        // Build a tensorproduct grid
        // Comment out, unless you build a tensorproduct grid
//        Dune::TensorGridFactory<GridType> factory;
//        for(int i=0; i<dim; ++i)
//            // This defines the coordinates in the ith direction to be a sequence
//            // that fills the interval [0,1] with 10 elements, such that the first
//            // mesh width is 0.5 and the mesh width follow a geometric series.
//            factory.geometricFillRange(i, 10, 1., 0.5);
//        auto grid = factory.createGrid();

        // refine the grid once globally
//        grid->globalRefine(1);

        // output some rudimentary information on the grid
        std::cout << "Successfully constructed a grid with" << std::endl;
        std::cout << "  " << grid->size(0) << " elements" << std::endl;
        std::cout << "  " << grid->size(dim) << " vertices" << std::endl;

        // Extract the grid view we would like to use
        auto gridview = grid->leafGridView();

        // get a vector with the indices of all elements iteration order
        std::vector<int> data(gridview.size(0));
        const auto& indexset = gridview.indexSet();
        for(const auto& e: elements(gridview))
            data[indexset.index(e)] = indexset.index(e);

        // write this information to a vtk writer
        std::string outputname = "elementdata";
        Dune::VTKWriter<GridType::LeafGridView> vtkwriter(gridview);
        vtkwriter.addCellData(data,"indexing");
        vtkwriter.write(outputname, Dune::VTK::appendedraw);
   }
    // catch exceptions
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
