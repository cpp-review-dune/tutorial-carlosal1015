/********************************************************/
// Beware of line number changes, they may corrupt docu!
//! \brief Driver function to set up and solve the problem
/********************************************************/

template<typename GV>
void driver (const GV& gv, Dune::ParameterTree& ptree)
{
  // dimension and important types
  const int dim = GV::dimension;
  typedef typename GV::Grid::ctype DF; // type for ccordinates
  typedef double RF;                   // type for computations

  // make user functions
  RF eta = ptree.get("problem.eta",(RF)1.0);
  Problem<RF> problem(eta);
  auto glambda = [&](const auto& e, const auto& x){return problem.g(e,x);};
  auto g = Dune::PDELab::makeGridFunctionFromCallable(gv,glambda);

  // Make grid function space
  typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> FEM;
  FEM fem(Dune::GeometryTypes::cube(dim));
  typedef Dune::PDELab::NoConstraints CON;
  typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  gfs.name("Q0");

  // A coefficient vector
  using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
  Z z(gfs); // initial value

  // Make a grid function out of it
  typedef Dune::PDELab::DiscreteGridFunction<GFS,Z> ZDGF;
  ZDGF zdgf(gfs,z);

  // Fill the coefficient vector
  Dune::PDELab::interpolate(g,gfs,z);

  // Make a local operator
  typedef NonlinearPoissonFV<Problem<RF> > LOP;
  LOP lop(problem);

  // Make a global operator
  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
  MBE mbe(2*dim+1); // guess nonzeros per row
  typedef Dune::PDELab::EmptyTransformation CC;
  typedef Dune::PDELab::GridOperator<
    GFS,GFS,  /* ansatz and test space */
    LOP,      /* local operator */
    MBE,      /* matrix backend */
    RF,RF,RF, /* domain, range, jacobian field type*/
    CC,CC     /* constraints for ansatz and test space */
    > GO;
  GO go(gfs,gfs,lop,mbe);

  // Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO> LS;
  LS ls(100,2);

  // solve nonlinear problem
  Dune::PDELab::NewtonMethod<GO,LS> newton(go,ls);
  newton.setParameters(ptree.sub("newton"));
  newton.apply(z);

  // Write VTK output file
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  typedef Dune::PDELab::VTKGridFunctionAdapter<ZDGF> VTKF;
  vtkwriter.addCellData(std::shared_ptr<VTKF>(new
                                         VTKF(zdgf,"fesol")));
  vtkwriter.write(ptree.get("output.filename","output"),
                  Dune::VTK::appendedraw);
}
