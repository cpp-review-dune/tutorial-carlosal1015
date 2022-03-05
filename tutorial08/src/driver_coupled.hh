//===============================================================
// Adapter for the transport problem
//===============================================================

template<typename GV, typename TBCType, typename TG, typename Velocity>
class GenericTransportProblem
{
  typedef typename TG::Traits::RangeFieldType RF;
  typedef Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type BCType;

  typename Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF>::PermTensorType I;
  RF time;
  const TBCType& Tbctype;
  TG Tg;
  const Velocity& velocity;
  
public:
  typedef Dune::PDELab::ConvectionDiffusionParameterTraits<GV,RF> Traits;

  GenericTransportProblem (Dune::ParameterTree& ptree, const TBCType& Tbctype_, TG Tg_, const Velocity& velocity_, double heatconductivity)
    : time(0.0), Tbctype(Tbctype_), Tg(Tg_), velocity(velocity_)
  {
    // diffusion tensor
    for (std::size_t i=0; i<Traits::dimDomain; i++)
      for (std::size_t j=0; j<Traits::dimDomain; j++)
        I[i][j] = (i==j) ? heatconductivity : 0;
  }

  //! yes it is constant
  static constexpr bool permeabilityIsConstantPerCell()
  {
    return true;
  }

  //! tensor diffusion coefficient
  typename Traits::PermTensorType
  A (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return I;
  }

  //! velocity field
  typename Traits::RangeType
  b (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Velocity::Traits::RangeType velo(0.0);
    velocity.evaluate(e,x,velo);
    return velo;
  }

  //! sink term
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! source term
  typename Traits::RangeFieldType
  f (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 0.0;
  }

  //! boundary condition type function
  /* return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Dirichlet for Dirichlet boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Neumann for flux boundary conditions
   * return Dune::PDELab::ConvectionDiffusionBoundaryConditions::Outflow for outflow boundary conditions
   */
  BCType
  bctype (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& xlocal) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(xlocal);
    return Tbctype(xglobal);
  }

  //! Dirichlet boundary condition value
  typename Traits::RangeFieldType
  g (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename TG::Traits::RangeType temperature;
    Tg.evaluate(e,x,temperature);
    return temperature;
  }

  //! flux boundary condition
  typename Traits::RangeFieldType
  j (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! outflow boundary condition
  typename Traits::RangeFieldType
  o (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x) const
  {
    return 0.0;
  }

  //! set time for subsequent evaluation
  void setTime (RF t)
  {
    time = t;
    Tg.setTime(t);
  }
};

//===============================================================
// a driver for solving the Navier-Stokes eqs.
//===============================================================

template<typename GV, typename Scheme, typename BCType, typename BConstraints, typename G, typename BCTypeT, typename GT>
void driver_coupled (const GV& gv, Scheme scheme,
		     const BCType& bctype, const BConstraints& bconstraints, G& g,BCTypeT& bctypeT, GT& gT, 
		     Dune::ParameterTree& ptree)
{
  // constants and types
  const int dim = GV::dimension;
  typedef typename GV::Grid::ctype DF;
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
  std::cout << "setting up flow" << std::endl;
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

  // make a vector valued grid function out of the velocity
  using Path0 = Dune::TypeTree::HybridTreePath<Dune::index_constant<0>>;
  using Path1 = Dune::TypeTree::HybridTreePath<Dune::index_constant<1>>;
  using VelocitySubGFS = Dune::PDELab::GridFunctionSubSpace<GFSTH,Path0>;
  VelocitySubGFS velocitysubgfs(gfsth); // subspace
  using VDGF = Dune::PDELab::VectorDiscreteGridFunction<VelocitySubGFS,Z>;
  VDGF vdgf(velocitysubgfs,z); // current velocity as a grid function

  //=============================
  // set up the transport problem
  //=============================

  // transport parameter class
  auto heatconductivity = ptree.get<RF>("problem.heatconductivity");
  typedef GenericTransportProblem<GV,BCTypeT,GT,VDGF> TP;
  TP tp(ptree,bctypeT,gT,vdgf,heatconductivity);

  // grid function space for transport problem
  typedef Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed,Scheme::nT> VBT;
  typedef Dune::PDELab::NoConstraints CONT;
  typedef Dune::PDELab::GridFunctionSpace<GV,typename Scheme::FEMT,CONT,VBT> GFST;
  GFST gfsT(gv,scheme.femT);
  gfsT.name("temperature");
  gfsT.ordering();
  std::cout << "temperature field has " << gfsT.globalSize() << " dofs" << std::endl;

  // coefficient vector
  using ZT = Dune::PDELab::Backend::Vector<GFST,RF>;
  ZT zT(gfsT);

  // initialize coefficient vector for transport with initial temperature field
  gT.setTime(time);
  Dune::PDELab::interpolate(gT,gfsT,zT);

  // make grid function
  typedef Dune::PDELab::DiscreteGridFunction<GFST,ZT> DGFT;
  DGFT dgfT(gfsT,zT);

  // compute force from temperature locally in the element
  auto rho_0 = ptree.get<RF>("problem.rho_0"); // rho = rho_0 - kappa*temperature
  auto alpha = ptree.get<RF>("problem.alpha");
  auto flambda = [&](const auto& entity, const auto& position) {
    Dune::FieldVector<RF,dim> f(0.0); // return value is the force vector
    Dune::FieldVector<DF,dim> x; for (int i=0; i<dim; i++) x[i] = position[i];
    typename DGFT::Traits::RangeType temperature;
    dgfT.evaluate(entity,x,temperature);
    f[dim-1] = -(rho_0 - alpha*temperature); // force is vertically down
    return f;
  };

  // =====================
  // continue with flow problem: local operator setup
  // =====================
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
  
  //=============================
  // continue with transport problem
  //=============================

  // assemblers for finite element problem
  typedef Dune::PDELab::ConvectionDiffusionDG<TP,typename Scheme::FEMT> LOPT;
  LOPT lopT(tp,Dune::PDELab::ConvectionDiffusionDGMethod::SIPG,Dune::PDELab::ConvectionDiffusionDGWeights::weightsOn,3.0);
  typedef typename GFST::template ConstraintsContainer<RF>::Type CCT;
  CCT ccT;
  
  typedef Dune::PDELab::GridOperator<GFST,GFST,LOPT,MatrixBackend,RF,RF,RF,CCT,CCT> GOST;
  GOST gosT(gfsT,ccT,gfsT,ccT,lopT,MatrixBackend(2*dim+1));
  typedef Dune::PDELab::L2 MLOPT;
  MLOPT mlopT(Scheme::degreeT+Scheme::degreeu+2);
  typedef Dune::PDELab::GridOperator<GFST,GFST,MLOPT,MatrixBackend,RF,RF,RF,CCT,CCT> MGOST;
  MGOST mgosT(gfsT,ccT,gfsT,ccT,mlopT,MatrixBackend(2*dim+1));
  typedef Dune::PDELab::OneStepGridOperator<GOST,MGOST> IGOT;
  IGOT igoT(gosT,mgosT);
  igoT.divideMassTermByDeltaT();

  // linear solver backends
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LinearSolverT;
  LinearSolverT lsT(lineariterationsmax,linearsolververbosity);

  // linear problem solver
  auto linearreduction = ptree.get<double>("solver.linearreduction");
  typedef Dune::PDELab::StationaryLinearProblemSolver<IGOT,LinearSolverT,ZT> PDESOLVERT;
  PDESOLVERT pdesolverT(igoT,lsT,linearreduction);

  // implicit time-stepper
  Dune::PDELab::Alexander2Parameter<RF> methodT;
  typedef Dune::PDELab::OneStepMethod<RF,IGOT,PDESOLVERT,ZT> OSMT;
  OSMT osmT(methodT,igoT,pdesolverT);
  osmT.setVerbosityLevel(2);

  //====================
  // prepare VTK writer
  //====================
  std::cout << "setting up output" << std::endl;
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
  typedef Dune::PDELab::VectorDiscreteGridFunctionDiv<VelocitySubGFS,Z> VDivDGF;
  VDivDGF vdivdgf(velocitysubgfs,z);
  vtkSequenceWriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<VDivDGF> >(vdivdgf,"div(v)"));
  // add temperature
  std::cout << "adding temperature to VTKWriter" << std::endl;
  vtkSequenceWriter.addVertexData(std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<DGFT> >(dgfT,"temperature"));
  std::cout << "write first file" << std::endl;
  vtkSequenceWriter.write(time,Dune::VTK::appendedraw);

  //////////////////////
  // subspace
  //////////////////////
  // using PressureSubGFS = Dune::PDELab::GridFunctionSubSpace<GFSTH,Path1>;
  // PressureSubGFS pressuresubgfs(gfsth);
  // auto pfunc = Dune::PDELab::DiscreteGridViewFunction<PressureSubGFS,Z>(pressuresubgfs,z);
  // auto plocal = localFunction(pfunc);
  // for (const auto& e : elements(gv))
  //   {
  //     plocal.bind(e);
  //     auto value = plocal(Dune::FieldVector<double,2>(1.0));
  //   }
  
  // the time loop
  const RF dt = ptree.get<RF>("problem.dt");
  const RF T = ptree.get<RF>("problem.T");
  const auto every = ptree.get<int>("output.every");
  int step=1;
  while (time<T-1e-9)
    {
      // do time step for flow
      Z znew(z);
      osm.apply(time,dt,z,g,znew);
      z = znew;
      
      // do time step for transport
      ZT zTnew(zT);
      osmT.apply(time,dt,zT,zTnew);
      zT = zTnew;

      // advance time
      time += dt;
      
      // store results
      if (step%every==0)
	vtkSequenceWriter.write(time,Dune::VTK::appendedraw);
      step++;
    }
}
