template <class Grid>
void size_of_the_grid(const Grid & grid)
{
    // get iterator types
    typedef typename Grid::LeafGridView::template Codim<0>::Iterator
        Iterator;
    typedef typename Grid::ctype ctype;

    ctype size = 0.0;

    // iterate over the leaf
    for (Iterator it = grid.leafGridView().template begin<0>();
         it != grid.leafGridView().template end<0>(); ++it) {
        // sum up the volumina
        size += it->geometry().volume();
    }
}
