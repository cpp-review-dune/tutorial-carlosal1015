template<class GridView>
void do_something(const GridView &grid)
{
    // iterate over the grid 
    for (auto entity : entities(gv,DUNE::Codim<0>))
    {
      ...
    }
}
