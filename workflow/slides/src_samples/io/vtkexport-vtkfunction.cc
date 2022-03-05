#include<dune/grid/io/file/vtk/vtkwriter.hh>

template<class GridView> class XY_VTKFunction
  : public Dune::VTKWriter<GridView>::VTKFunction
{
    dune_static_assert(GridView::dimension == 2,
                       "Illegal GridView dimension");

  public:
    typedef typename GridView::template Codim<0>::Entity Entity;
    typedef Dune::FieldVector<typename GridView::ctype,
            GridView::dimension> CoordType;

    virtual int ncomps() const { return 1; }
    virtual std::string name() const { return "xy_elementwise"; }
    virtual double evaluate(int comp, const Entity &e,
                            const CoordType &xi) const
    {
        auto coord=e.geometry().global(xi);
        return coord[0]*coord[1];
    }
};

template<class GridView>
void functiondata(const GridView& gridview) {
    Dune::VTKWriter<GridView> vtkwriter(gridview);
    vtkwriter.addVertexData(
        Dune::make_shared<XY_VTKFunction<GridView> >());
    vtkwriter.write("functiondata",Dune::VTK::ascii);
}
