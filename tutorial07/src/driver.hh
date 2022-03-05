// -*- tab-width: 4; indent-tabs-mode: nil -*-

//===============================================================
// driver for general purpose hyperbolic solver
//===============================================================

template<typename GV, typename FEMDG, typename NUMFLUX>
void driver (const GV& gv, const FEMDG& femdg, NUMFLUX& numflux, Dune::ParameterTree& ptree)
{

  /// tex: extract
  // Choose domain and range field type
  using RF = typename NUMFLUX::RF; // type for computations
  static constexpr int dim = NUMFLUX::dim;
  static constexpr int m = NUMFLUX::m; //number of components
  /// tex: extract


  //initial condition
  auto u0lambda = [&](const auto& i, const auto& x)
    {return numflux.model().problem.u0(i,x);};
  auto u0 = Dune::PDELab::
    makeGridFunctionFromCallable(gv,u0lambda);

  typedef Dune::PDELab::NoConstraints CON;

  /// tex: gfs
  //block size for a component deducted by pdelab
  using VBE0 = Dune::PDELab::ISTL::VectorBackend<>;

  using VBE = Dune::PDELab::ISTL::VectorBackend<
    Dune::PDELab::ISTL::Blocking::fixed>;
  using OrderingTag = Dune::PDELab::EntityBlockedOrderingTag;

  using GFSDG = Dune::PDELab::GridFunctionSpace<GV,FEMDG,CON,VBE0>;
  GFSDG gfsdg(gv,femdg);

  using GFS = Dune::PDELab::PowerGridFunctionSpace<
    GFSDG,m,VBE,OrderingTag>;
  GFS gfs(gfsdg);

  typedef typename GFS::template ConstraintsContainer<RF>::Type C;
  C cg;
  gfs.update(); // initializing the gfs
  std::cout << "degrees of freedom: " << gfs.globalSize() << std::endl;
  /// tex: gfs

  /// tex: lop
  // Make instationary grid operator
  using LOP = Dune::PDELab::DGHyperbolicSpatialOperator<NUMFLUX,FEMDG>;
  LOP lop(numflux);
  using TLOP = Dune::PDELab::DGHyperbolicTemporalOperator<NUMFLUX,FEMDG>;
  TLOP tlop(numflux);
  /// tex: lop

  using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  MBE mbe(2*dim+1); // Maximal number of nonzeroes per row

  using GO0 = Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,C,C>;
  GO0 go0(gfs,cg,gfs,cg,lop,mbe);
  using GO1 = Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,RF,RF,RF,C,C>;
  GO1 go1(gfs,cg,gfs,cg,tlop,mbe);
  using IGO = Dune::PDELab::OneStepGridOperator<GO0,GO1,false>;
  IGO igo(go0,go1);

  /// tex: timestepping
  // select and prepare time-stepping scheme
  int torder = ptree.get("fem.torder",(int)1);

  Dune::PDELab::ExplicitEulerParameter<RF> method1;
  Dune::PDELab::HeunParameter<RF> method2;
  Dune::PDELab::Shu3Parameter<RF> method3;
  Dune::PDELab::RK4Parameter<RF> method4;
  Dune::PDELab::TimeSteppingParameterInterface<RF> *method;

  if (torder==1) {
    method=&method1;
    std::cout << "setting explicit Euler" << std::endl;}
  if (torder==2) {
    method=&method2;
    std::cout << "setting Heun" << std::endl;}
  if (torder==3) {
    method=&method3;
    std::cout << "setting Shu 3" << std::endl;}
  if (torder==4) {
    method=&method4;
    std::cout << "setting RK4" << std::endl;}
  if (torder<1||torder>4)
    std::cout << "torder should be in [1,4]" << std::endl;
  /// tex: timestepping

  igo.setMethod(*method);

  // set initial values
  typedef typename IGO::Traits::Domain V;
  V xold(gfs,0.0);

  Dune::PDELab::interpolate(u0,gfs,xold);

  // Make a linear solver backend
  typedef Dune::PDELab::ISTLBackend_OVLP_ExplicitDiagonal<GFS> LS;
  LS ls(gfs);
  /// tex: lsosm
  //typedef Dune::PDELab::ISTLBackend_SEQ_UMFPack LS;
  //LS ls(false);

  // time-stepping
  Dune::PDELab::ExplicitOneStepMethod<RF,IGO,LS,V,V> osm(*method,igo,ls);
  osm.setVerbosityLevel(2);
  /// tex: lsosm

  // prepare VTK writer and write first file
  int subsampling=ptree.get("output.subsampling",(int)0);
  using VTKWRITER = Dune::SubsamplingVTKWriter<GV>;
  //VTKWRITER vtkwriter(gv,subsampling);
  VTKWRITER vtkwriter(gv,Dune::refinementLevels(subsampling));

  std::string filename=ptree.get("output.filename","output");

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
  Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter,gfs,xold,Dune::PDELab::vtk::DefaultFunctionNameGenerator("u"));
  vtkSequenceWriter.write(0.0,Dune::VTK::appendedraw);

  // initialize simulation time
  RF time = 0.0;

  // time loop
  RF T = ptree.get("problem.T",(RF)1.0);
  RF dt = ptree.get("fem.dt",(RF)0.1);
  const int every = ptree.get<int>("output.every");
  int counter = 0;

  V x(gfs,0.0);

  while (time < T)
    {
      // do time step
      osm.apply(time,dt,xold,x);

      // accept time step
      xold = x;
      time += dt;
      counter++;

      //output to VTK file every n-th timestep
      if(counter % every == 0)
        vtkSequenceWriter.write(time,Dune::VTK::appendedraw);
    }
}
