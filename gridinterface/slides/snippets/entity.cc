template<class GridView>
void do_something(const GridView & grid)
{
    // get types
    typedef typename GridView::template Codim<0>::Iterator
        Iterator;
    typedef typename GridView::Grid::template Codim<0>::Entity
        Entity;
    // iterate over the leaf
    for (Iterator it = grid.template begin<0>();
         it != grid.template end<0>(); ++it) {
    Entity & e = *it;
    }
}
