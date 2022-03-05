#include <dune/grid/somegrid.hh>

void iterate_the_grid()
{
    // we have a grid
    Dune::SomeGrid grid(parameters);

    // iterate over the interior Border Partition of level 2
    auto levelGV = grid.levelGridView(2);
    for (const auto &face : facets(levelGV,Dune::Partitions::InteriorBorder) {
        ...
    }

    // iterate over all partition of the leaf
    auto leafGV = grid.leafGridView();
    for (const auto &node : entities(leafGV,Dune::Codim<DIM>,Dune::Partitions::All) {
        ...
    }
}
