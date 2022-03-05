/********************************************************/
// Beware of line number changes, they may corrupt docu!
//! \brief Driver function to set up and solve the problem
/********************************************************/

template<typename Grid, int degree>
void driver (Grid& grid, Dune::ParameterTree& ptree)
{
  // get leaf gridview
  auto gv=grid.leafGridView();
  typedef decltype(gv) GV;

  // dimension and important types
  const int dim = GV::dimension;
  typedef typename GV::Grid::ctype DF; // type for coordinates
  typedef double RF;                   // type for computations

  // make user functions
  RF eta = ptree.get("problem.eta",(RF)1.0);
  Problem<RF> problem(eta);
  auto glambda =
    [&](const auto& e, const auto& x){return problem.g(e,x);};
  auto g = Dune::PDELab::makeGridFunctionFromCallable(gv,glambda);
  auto blambda =
    [&](const auto& i, const auto& x){return problem.b(i,x);};
  auto b = Dune::PDELab::makeBoundaryConditionFromCallable(gv,blambda);

  // Make grid function space
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,RF,degree> FEM;
  FEM fem(gv);
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  gfs.name("Vh");

  // Assemble constraints
  typedef typename GFS::template ConstraintsContainer<RF>::Type CC;
  CC cc;
  Dune::PDELab::constraints(b,gfs,cc); // assemble constraints

  // A coefficient vector
  using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
  Z z(gfs); // initial value

  // Fill the coefficient vector
  Dune::PDELab::interpolate(g,gfs,z);

  // adaptation loop
  // only grid function space and coefficient vector
  // live outside this loop
  int steps = ptree.get("fem.steps",(int)3);
  int uniformlevel = ptree.get("fem.uniformlevel",(int)2);
  for (int i=0; i<steps; i++)
  {
    std::stringstream s;
    s << i;
    std::string iter;
    s >> iter;
    std::cout << "Iteration: " << iter
      << "\thighest level in grid: " << grid.maxLevel()
      << std::endl;
    std::cout << "constrained dofs=" << cc.size()
      << " of " << gfs.globalSize() << std::endl;

    // Make local operator
    typedef NonlinearPoissonFEM<Problem<RF>,FEM> LOP;
    LOP lop(problem);

    // Make a global operator
    typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
    MBE mbe((int)pow(1+2*degree,dim));
    typedef Dune::PDELab::GridOperator<
      GFS,GFS,  /* ansatz and test space */
      LOP,      /* local operator */
      MBE,      /* matrix backend */
      RF,RF,RF, /* domain, range, jacobian field type */
      CC,CC     /* constraints for ansatz and test space */
      > GO;
    GO go(gfs,cc,gfs,cc,lop,mbe);

    // Select a linear solver backend
    typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO> LS;
    LS ls(100,0);

    // solve nonlinear problem
    Dune::PDELab::NewtonMethod<GO,LS> newton(go,ls);
    newton.setParameters(ptree.sub("newton"));
    newton.apply(z);

    // set up error estimator
    typedef Dune::PDELab::P0LocalFiniteElementMap<DF,RF,dim> P0FEM;
    P0FEM p0fem(Dune::GeometryTypes::simplex(dim));
    typedef Dune::PDELab::NoConstraints NCON;
    typedef Dune::PDELab::GridFunctionSpace<GV,P0FEM,NCON,VBE> P0GFS;
    P0GFS p0gfs(gv,p0fem);
    typedef NonlinearPoissonFEMEstimator<Problem<RF>,FEM> ESTLOP;
    ESTLOP estlop(problem);
    typedef Dune::PDELab::EmptyTransformation NCC;
    typedef Dune::PDELab::GridOperator<
      GFS,P0GFS,    /* one value per element */
      ESTLOP,       /* operator for error estimate */
      MBE,RF,RF,RF, /* same as before */
      NCC,NCC       /* no constraints */
      > ESTGO;
    ESTGO estgo(gfs,p0gfs,estlop,mbe);

    // compute local error contribution and global error
    using Z0 = Dune::PDELab::Backend::Vector<P0GFS,RF>;
    Z0 z0(p0gfs,0.0);
    estgo.residual(z,z0);
    auto estimated_error = sqrt(z0.one_norm());
    std::cout << "Estimated error in step " << i
      << " is " << estimated_error << std::endl;

    // vtk output
    std::cout << "VTK output" << std::endl;
    Dune::SubsamplingVTKWriter<GV>
      vtkwriter(gv,Dune::refinementIntervals(ptree.get("output.subsampling",(int)1)));
    typedef Dune::PDELab::DiscreteGridFunction<GFS,Z> ZDGF;
    ZDGF zdgf(gfs,z);
    typedef Dune::PDELab::VTKGridFunctionAdapter<ZDGF> VTKF;
    vtkwriter.addVertexData(
        std::shared_ptr<VTKF>(new VTKF(zdgf,"fesol")));
    typedef Dune::PDELab::DiscreteGridFunction<P0GFS,Z0> Z0DGF;
    Z0DGF z0dgf(p0gfs,z0);
    typedef Dune::PDELab::VTKGridFunctionAdapter<Z0DGF> VTKF0;
    vtkwriter.addCellData(
        std::shared_ptr<VTKF0>(new VTKF0(z0dgf,"error2")));
    vtkwriter.write(ptree.get("output.filename",
          (std::string)"output")+iter, Dune::VTK::appendedraw);

    // error control
    auto tol = ptree.get("fem.tol",(double)0.0);
    if (estimated_error<=tol) break;
    if (i==steps-1) break;

    // mark elements for refinement
    std::cout << "mark elements" << std::endl;
    auto fraction = ptree.get("fem.fraction",(double)0.5);
    RF eta_refine,eta_coarsen;
    Dune::PDELab::error_fraction(
        z0,fraction,0.0,eta_refine,eta_coarsen);
    if (fraction>=1.0 || i<uniformlevel) eta_refine=0.0;
    Dune::PDELab::mark_grid(grid,z0,eta_refine,0.0,2);

    // do refinement
    std::cout << "adapt grid and solution" << std::endl;
    Dune::PDELab::adapt_grid(grid,gfs,z,2*(degree+1));

    // recompute constraints
    std::cout << "constraints and stuff" << std::endl;
    Dune::PDELab::constraints(b,gfs,cc);

    // write correct boundary conditions in new vector
    Z znew(gfs);
    Dune::PDELab::interpolate(g,gfs,znew);

    // copy Dirichlet boundary to interpolated solution
    Dune::PDELab::copy_constrained_dofs(cc,znew,z);
  }
}
