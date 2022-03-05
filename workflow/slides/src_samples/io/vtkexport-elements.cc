template<class GridView>
void elementdata(const GridView& gridview)
{
    // allocate a vector for the data
    std::vector<double> solution(gridview.size(0));

    // iterate through all entities of codim 0
    for (const auto& e : elements(gridview))
    {
        // get global coordinate of cell center
        auto global = e.geometry().center();
        // evaluate function and store value
        solution[gridview.indexSet().index(e)] = exp(global[0]*global[1]);
    }

    // generate a VTK file
    Dune::VTKWriter<GridView> vtkwriter(gridview);
    vtkwriter.addCellData(solution, "data");
    vtkwriter.write("elementdata", Dune::VTK::appendedraw);
}
