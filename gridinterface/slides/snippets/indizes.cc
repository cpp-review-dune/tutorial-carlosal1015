#include <dune/grid/common/mcmgmapper.hh>
...

typedef Dune::SomeGrid::LeafGridView GridView;
typedef GridView::template Codim<0>::Iterator Iterator;
typedef GridView::IntersectionIterator
    IntersectionIterator;
...

/* create a mapper*/
// Layout description (equivalent to Dune::MCMGElementLayout)
template<int dim>
struct CellData {
    bool contains (Dune::GeometryType gt) {
        return gt.dim() == dim;
    }
};

// mapper for elements (codim=0) on leaf
typedef
Dune::MultipleCodimMultipleGeomTypeMapper<GridView,CellData> Mapper;
Mapper mapper(gridview);

/* setup sparsity pattern */
// iterate over the leaf
for (Iterator it = gridview.template begin<0>();
     it != gridview.template end<0>(); ++it)
{
    int index = mapper.map(*it);

    IntersectionIterator iitend = gridview.iend(*it);
    for (IntersectionIterator iit = gridview.ibegin(*it);
         iit != iitend; ++iit) {
        // neighbor intersection
        if (iit->neighbor()) {
            int nindex = mapper.map(*(iit->outside()));
            matrix[index].insert(nindex);
        }
    }
}
