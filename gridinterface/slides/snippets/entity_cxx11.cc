template<class GridView>
void do_something(const GridView & grid)
{
    // iterate over the leaf
    for (auto it = grid.template begin<0>();
         it != grid.template end<0>(); ++it) {
        // get (reference to) entity
        // type is GridView::Grid::Codim<0>::Entity
        auto e = *it;
    }
}
