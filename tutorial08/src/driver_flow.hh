//===============================================================
// a driver for solving the Navier-Stokes eqs.
//===============================================================

template<typename GV, typename Scheme, typename BCType, typename BConstraints, typename G>
void driver_flow (const GV& gv, const Scheme& scheme,
	     const BCType& bctype, const BConstraints& bconstraints, G& g,
	     Dune::ParameterTree& ptree)
{
  // constants and types
  const int dim = GV::dimension;
  typedef double RF;                   // type for computations
  const int degreeu = Scheme::degreeu;
  const int degreep = Scheme::degreep;
  const int nu = Scheme::nu;
  const int np = Scheme::np;
  const int mvol = Scheme::mvol;
  const int mbnd = Scheme::mbnd;
  const int faces = Scheme::faces;
  std::cout << "dim=" << dim << " degreeu=" << degreeu
	    << " degreep=" << degreep << " nu=" << nu
	    << " np=" << np << std::endl;  

  // grid function space for velocity
  typedef typename Scheme::FEMU FEMU;
  typedef Dune::PDELab::ISTL::VectorBackend<> VectorBackend;
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  typedef Dune::PDELab::VectorGridFunctionSpace
    <GV,FEMU,dim,VectorBackend,VectorBackend,CON> GFSU;
  GFSU gfsu(gv,scheme.femu);
  gfsu.name("velocity");

  // grid function space for pressure
  typedef typename Scheme::FEMP FEMP;
  typedef Dune::PDELab::GridFunctionSpace
    <GV,FEMP,CON,VectorBackend> GFSP;
  GFSP gfsp(gv,scheme.femp);
  gfsp.name("pressure");

  // combined Taylor-Hood grid function space
  typedef Dune::PDELab::CompositeGridFunctionSpace<
    VectorBackend,
    Dune::PDELab::LexicographicOrderingTag,
    GFSU,GFSP> GFSTH;
  GFSTH gfsth(gfsu,gfsp);
  gfsth.name("Taylor-Hood");
  // std::cout << "N=" << gfsth.globalSize() << std::endl; gibt einen Fehler hier ...
  //std::cout << "Nu=" << gfsu.globalSize() << " Np=" << gfsp.globalSize() << std::endl;
  
  // assemble the constraints using the given constraints function
  typedef typename GFSTH::template ConstraintsContainer<RF>::Type CC;
  CC cc;
  Dune::PDELab::constraints(bconstraints,gfsth,cc); // assemble constraints
  std::cout << "constrained dofs=" << cc.size() << " from " << gfsth.globalSize() << std::endl;

  // set up the coefficient vector
  using Z = Dune::PDELab::Backend::Vector<GFSTH,RF>;
  Z z(gfsth);

  // initialize time
  RF time = 0.0;

  // initialize coefficient vector
  g.setTime(time);
  Dune::PDELab::interpolate(g,gfsth,z);

  // local operator setup
  auto flambda = [](const auto& entity, const auto& position) { // right hand side function for momentum equation
    Dune::FieldVector<RF,GV::dimension> f(0.0); 
    return f;
  };
  const RF viscosity = ptree.get<RF>("problem.viscosity");
  const RF gamma = ptree.get<RF>("problem.gamma");
  const RF beta = ptree.get<RF>("problem.beta");
  using LOP = Dune::PDELab::ConformingNavierStokesLOP<decltype(flambda),BCType,FEMU,FEMP,degreeu,degreep,nu,np,mvol,mbnd,faces>;
  LOP lop(flambda,bctype,scheme.femu,scheme.femp,gv,viscosity,gamma,beta);
  using MLOP = Dune::PDELab::ConformingNavierStokesMassLOP<FEMU,degreeu,nu,mvol>;
  MLOP mlop(scheme.femu,gv);

  // grid operator setup
  typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MatrixBackend;
  typedef Dune::PDELab::GridOperator<GFSTH,GFSTH,LOP,MatrixBackend,RF,RF,RF,CC,CC> GOS;
  GOS gos(gfsth,cc,gfsth,cc,lop,MatrixBackend(66));
  typedef Dune::PDELab::GridOperator<GFSTH,GFSTH,MLOP,MatrixBackend,RF,RF,RF,CC,CC> MGOS;
  MGOS mgos(gfsth,cc,gfsth,cc,mlop,MatrixBackend(66));
  typedef Dune::PDELab::OneStepGridOperator<GOS,MGOS> IGO;
  IGO igo(gos,mgos);
  igo.divideMassTermByDeltaT();

  // typedef typename GOS::Jacobian M;
  // M m(gos);
  // std::cout << m.patternStatistics() << std::endl;
  // gos.jacobian(z,m);
  // return;

  // linear solver setup
  auto lineariterationsmax = ptree.get<int>("solver.lineariterationsmax");
  auto linearsolververbosity = ptree.get<int>("solver.linearsolververbosity");

  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_ILU0 LinearSolver;
  LinearSolver ls(lineariterationsmax,linearsolververbosity);

  // nonlinear solver setup
  typedef Dune::PDELab::NewtonMethod<IGO,LinearSolver> Newton;
  Newton newton(igo,ls);
  newton.setParameters(ptree.sub("newton"));

  // time integrator setup
  Dune::PDELab::OneStepThetaParameter<RF> crank_nicolson(0.5);
  Dune::PDELab::Alexander2Parameter<RF> alexander;
  Dune::PDELab::OneStepMethod<RF,IGO,Newton,Z,Z> osm(alexander,igo,newton);

  // prepare VTK writer
  int subsampling=ptree.get("output.subsampling",(int)0);
  using VTKWRITER = Dune::SubsamplingVTKWriter<GV>;
  VTKWRITER vtkwriter(gv,Dune::refinementLevels(subsampling));
  std::string filename=ptree.get("output.basename","output");
  struct stat st;
  if( stat( filename.c_str(), &st ) != 0 )
    {
      int stat = 0;
      stat = mkdir(filename.c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
      if( stat != 0 && stat != -1)
        std::cout << "Error: Cannot create directory "
                  << filename << std::endl;
    }
  using VTKSEQUENCEWRITER = Dune::VTKSequenceWriter<GV>;
  VTKSEQUENCEWRITER vtkSequenceWriter(
    std::make_shared<VTKWRITER>(vtkwriter),filename,filename,"");
  // add data field for all components of the space to the VTK writer
  Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter,gfsth,z,Dune::PDELab::vtk::DefaultFunctionNameGenerator("th"));
  // add divergence
  using Path0 = Dune::TypeTree::HybridTreePath<Dune::index_constant<0>>;
  using VelocitySubGFS = Dune::PDELab::GridFunctionSubSpace<GFSTH,Path0>;
  VelocitySubGFS velocitysubgfs(gfsth);
  using VDivDGF = Dune::PDELab::VectorDiscreteGridFunctionDiv<VelocitySubGFS,Z>;
  VDivDGF vdivdgf(velocitysubgfs,z);
  vtkSequenceWriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<VDivDGF> >(vdivdgf,"div(v)"));

  vtkSequenceWriter.write(time,Dune::VTK::appendedraw);

  // the time loop
  const RF dt = ptree.get<RF>("problem.dt");
  const RF T = ptree.get<RF>("problem.T");
  const auto every = ptree.get<int>("output.every");
  int step=1;
  while (time<T-1e-9)
    {
      // do time step
      Z znew(z);
      osm.apply(time,dt,z,g,znew);
      
      // advance time
      z = znew;
      time += dt;
      
      // store results
      if (step%every==0)
	vtkSequenceWriter.write(time,Dune::VTK::appendedraw);
      step++;
    }
}
