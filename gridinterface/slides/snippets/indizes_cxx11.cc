#include <dune/grid/common/mcmgmapper.hh>
...

typedef Dune::SomeGrid::LeafGridView GridView;
...

/* create a mapper*/
// Layout description
Dune::MCMGLayout layout =
     [](Dune::GeometryType gt, int griddim) {
         return gt.dim() == griddim;
     };

// mapper for elements (codim=0) on leaf
using Mapper =
    Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
Mapper mapper(gridview,layout);

/* setup sparsity pattern */
// iterate over the leaf
for (const auto& entity : elements(gridview))
{
    int index = mapper.index(entity);
    // iterate over all intersections of this cell
    for (const auto& i : intersections(gridview,entity))
    {
        // neighbor intersection
        if (i.neighbor()) {
            int nindex = mapper.index(i.outside());
            matrix[index].insert(nindex);
        }
    }
}
