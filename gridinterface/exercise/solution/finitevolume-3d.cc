// finitevolume.cc
// Implementation of the explicit finite volume method
// for the transport equation

#include "config.h"

#include <array>
#include <iostream>

#include <dune/common/filledarray.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include "fv.hh"

// Class encapsulating the equation data
// The velocity field is constant in dircetion (1, 1).  There is no
// inflow from the boundary.
template<typename ctype, int dim, typename Vars>
class SimpleTransport
{
public:
  // the velocity field
  Dune::FieldVector<double,dim> u(const Dune::FieldVector<ctype,dim>& global) const
  {
    Dune::FieldVector<double,dim> r(0.0);
    r = 1.0;
    return r;
  }

  // Dirichlet boundary condition on inflow boundaries
  double c_in(const Dune::FieldVector<ctype,dim>& global) const
  {
    return 0.0;
  }
};

// initialize unknowns with initial value
template<typename GV, typename V>
void initialize (const GV& grid_view, V& c)
{
  // create mapper
  Dune::MultipleCodimMultipleGeomTypeMapper<GV> mapper(grid_view, Dune::mcmgElementLayout());

  for (const auto& cell : elements(grid_view))
    {
      // get global coordinate of cell center
      auto global = cell.geometry().center();

      if (global[0]>=0.1 && global[0]<=0.2 &&
          global[1]>=0.1 && global[1]<=0.2 &&
          global[2]>=0.1 && global[2]<=0.2)
        c[mapper.index(cell)] = 1.0;
      else
        c[mapper.index(cell)] = 0.0;
    }
}

int main(int argc, char** argv)
{
  try{

    // initialize MPI
    Dune::MPIHelper::instance(argc, argv);

    // instantiate YaspGrid on the unit square with 100 by 100 cells:
    const int dim = 3;
    using Grid = Dune::YaspGrid<dim>;

    // the physical extent (here: 1.0 in each direction)
    auto L = Dune::FieldVector<double,dim>(1.0);

    // the number of cells per direction
    auto N = Dune::filledArray<dim, int>(32);

    // We don't want a periodic grid in any direction
    auto periodic = std::bitset<dim>(false);

    Grid grid(L, N, periodic, 0);

    using GridView = Grid::LeafGridView;
    GridView grid_view = grid.leafGridView();

    // create concentration vector
    using ScalarField = std::vector<double>;
    ScalarField c(grid_view.size(0));

    // initialize concentration
    initialize(grid_view, c);

    using Parameters = SimpleTransport<double,GridView::dimension,ScalarField>;
    Parameters params;;

    FiniteVolume<GridView,Parameters,ScalarField> fv(grid_view, params);

    // Create VTKSequenceWriter for graphical output
    Dune::VTKSequenceWriter<GridView> vtk_writer(
      std::make_shared<Dune::VTKWriter<GridView> >(grid_view),
      "concentration"
      );
    // add cell data to vtk writer
    vtk_writer.addCellData(c,"c");

    const double dt = 0.3*(L[0]/N[0]);
    double t = 0.0;
    for (int step = 0; step < 3 * N[0]; ++step)
      {
        // output solution:
        std::cout << "Solution at timestep " << step << ", time " << t << std::endl;

        // use the write() method to write the data into the file
        vtk_writer.write(t);

        // explicit update solution from t to t+dt
        fv.update_concentration(c,dt);

        // update time
        t += dt;
      }
    return 0;
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
