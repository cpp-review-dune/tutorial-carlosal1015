/********************************************************/
// Beware of line number changes, they may corrupt docu!
//! \brief Driver function to set up and solve the problem
/********************************************************/

template<class GV>
void driver (const GV& gv, Dune::ParameterTree& ptree)
{
  // dimension and important types
  const int dim = GV::dimension;
  typedef typename GV::Grid::ctype DF; // type for ccordinates
  typedef double RF;                   // type for computations

  // make user functions
  auto flambda = [](const auto& x){
    return Dune::FieldVector<RF,1>(-2.0*x.size());};
  auto f = Dune::PDELab::
    makeGridFunctionFromCallable(gv,flambda);
  auto glambda = [](const auto& x){
    RF s=0.0; for (std::size_t i=0; i<x.size(); i++) s+=x[i]*x[i];
    return s;};
  auto g = Dune::PDELab::
    makeGridFunctionFromCallable(gv,glambda);
  auto blambda = [](const auto& x){return true;};
  auto b = Dune::PDELab::
    makeBoundaryConditionFromCallable(gv,blambda);

  // Make grid function space
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,RF,1> FEM;
  FEM fem(gv);
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  gfs.name("P1");

  // Assemble constraints
  typedef typename GFS::template
    ConstraintsContainer<RF>::Type CC;
  CC cc;
  Dune::PDELab::constraints(b,gfs,cc); // assemble constraints
  std::cout << "constrained dofs=" << cc.size() << " of "
            << gfs.globalSize() << std::endl;

  // A coefficient vector
  using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
  Z z(gfs); // initial value

  // Make a grid function out of it
  typedef Dune::PDELab::DiscreteGridFunction<GFS,Z> ZDGF;
  ZDGF zdgf(gfs,z);

  // Fill the coefficient vector
  Dune::PDELab::interpolate(g,gfs,z);
  Dune::PDELab::set_nonconstrained_dofs(cc,0.0,z);

  // Make a local operator
  typedef PoissonP1<decltype(f),FEM> LOP;
  LOP lop(f,fem.find(*gv.template begin<0>()));

  // Make a global operator
  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  MBE mbe(1<<(dim+1)); // guess nonzeros per row
  typedef Dune::PDELab::GridOperator<
    GFS,GFS,  /* ansatz and test space */
    LOP,      /* local operator */
    MBE,      /* matrix backend */
    RF,RF,RF, /* domain, range, jacobian field type*/
    CC,CC     /* constraints for ansatz and test space */
    > GO;
  GO go(gfs,cc,gfs,cc,lop,mbe);

  // Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO> LS;
  LS ls(100,3);

  // Assemble and solve linear problem
  typedef Dune::PDELab::
    StationaryLinearProblemSolver<GO,LS,Z> SLP;
  SLP slp(go,ls,z,1e-10);
  slp.apply(); // here all the work is done!

  // Write VTK output file
  Z w(gfs); // Lagrange interpolation of exact solution
  Dune::PDELab::interpolate(g,gfs,w);
  ZDGF wdgf(gfs,w);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  typedef Dune::PDELab::VTKGridFunctionAdapter<ZDGF> VTKF;
  vtkwriter.addVertexData(std::shared_ptr<VTKF>(new
                                         VTKF(zdgf,"fesol")));
  vtkwriter.addVertexData(std::shared_ptr<VTKF>(new
                                         VTKF(wdgf,"exact")));
  vtkwriter.write(ptree.get("output.filename","output"),
                  Dune::VTK::appendedraw);
}
