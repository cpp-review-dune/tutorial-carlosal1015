#include<dune/grid/io/file/vtk/vtkwriter.hh>

template<class GridView, class Functor>
void vertexdata(const GridView& gridview, const Functor& f) {
  // allocate a vector for the data
  std::vector<double> solution(gridview.size(GridView::dimension));

  // iterate through all entities of codim dim
  for (const auto& v : vertices(gridview))
    {
      // evaluate functor and store value
      solution[gridview.indexSet().index(v)] = f(v.geometry().corner(0));
    }

  // generate a VTK file
  Dune::VTKWriter<GridView> vtkwriter(gridview);
  vtkwriter.addVertexData(solution,"data");
  vtkwriter.write("vertexdata",Dune::VTK::appendedraw);
}
