#include <dune/grid/somegrid.hh>

void iterate_the_grid()
{
    // we have a grid
    Dune::SomeGrid grid(parameters);
    // get iterator types
    typedef Dune::SomeGrid::
        LeafGridView::Codim<0>::Iterator LeafIterator;
    typedef Dune::SomeGrid::
        LevelGridView::Codim<0>::Iterator LevelIterator;

    // iterate over level 2
    for (LevelIterator it = grid.levelGridView(2).template begin<0>();
         it != grid.levelGridView(2).template end<0>(); ++it) {
        ...
    }
    // iterate over the leaf
    for (LeafIterator it = grid.leafGridView().template begin<0>();
         it != grid.leafGridView().template end(); ++it) {
        ...
    }
}
