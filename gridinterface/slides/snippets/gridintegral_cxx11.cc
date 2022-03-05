template <class Grid>
void size_of_the_grid(const Grid & grid)
{
    // get coordinate types
    typedef typename Grid::ctype ctype;
    ctype size = 0.0;

    // iterate over the leaf
    for (auto it = grid.leafGridView().template begin<0>();
         it != grid.leafGridView().template end<0>(); ++it) {
        // sum up the volumina
        size += it->geometry().volume();
    }
}
