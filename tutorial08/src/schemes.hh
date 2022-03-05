template<typename GridView>
class TaylorHood_21_Triangle
{
public:
  // define constants
  enum {dim=GridView::dimension}; // the dimension
  enum {degreeu=2}; // velocity polynomial degree
  enum {degreep=1}; // pressure polynomial degree
  enum {degreeT=2}; // temperature polynomial degree

  enum { mvol=12}; // number of quadrature points on element 
  enum { mbnd=3};  // number of quadrature points on face
  enum { faces=dim+1};  // number of faces

  enum { nu=Dune::PB::PkSize<degreeu,dim>::value};  // # velocity degrees of freedom
  enum { np=Dune::PB::PkSize<degreep,dim>::value};  // # pressure degrees of freedom

  // some types
  typedef GridView GV;        // export also the grid view type
  typedef typename GV::Grid::ctype DF; // type for coordinates
  typedef double RF;          // representation type for solution values

  // finite element map types
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,RF,degreeu> FEMU; // finite element map for velocity
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,RF,degreep> FEMP; // finite element map for pressure

  // finite element map objects
  FEMU femu;
  FEMP femp;

  // DG space for transport
#if HAVE_GMP
  typedef Dune::PDELab::OPBLocalFiniteElementMap<DF,RF,degreeT,dim,Dune::GeometryType::simplex,Dune::GMPField<512>,Dune::PB::BasisType::Pk> FEMT;
#else
  typedef Dune::PDELab::OPBLocalFiniteElementMap<DF,RF,degreeT,dim,Dune::GeometryType::simplex> FEMT;
#endif
  enum { nT=Dune::PB::PkSize<degreeT,dim>::value};  // # temperature degrees of freedom per element

  FEMT femT;

  // Constructor makes objects
  TaylorHood_21_Triangle (const GridView& gv)
    : femu(gv), femp(gv)
  {}
};


template<typename GridView>
class TaylorHood_32_Triangle
{
public:
  // define constants
  enum {dim=GridView::dimension}; // the dimension
  enum {degreeu=3}; // velocity polynomial degree
  enum {degreep=2}; // pressure polynomial degree
  enum {degreeT=3}; // temperature polynomial degree

  enum { mvol=25}; // number of quadrature points on element 
  enum { mbnd=5};  // number of quadrature points on face
  enum { faces=dim+1};  // number of faces

  enum { nu=Dune::PB::PkSize<degreeu,dim>::value};  // # velocity degrees of freedom
  enum { np=Dune::PB::PkSize<degreep,dim>::value};  // # pressure degrees of freedom

  // some types
  typedef GridView GV;        // export also the grid view type
  typedef typename GV::Grid::ctype DF; // type for coordinates
  typedef double RF;          // representation type for solution values

  // finite element map types
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,RF,degreeu> FEMU; // finite element map for velocity
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,RF,degreep> FEMP; // finite element map for pressure

  // finite element map objects
  FEMU femu;
  FEMP femp;

  // DG space for transport
#if HAVE_GMP
  typedef Dune::PDELab::OPBLocalFiniteElementMap<DF,RF,degreeT,dim,Dune::GeometryType::simplex,Dune::GMPField<512>,Dune::PB::BasisType::Pk> FEMT;
#else
  typedef Dune::PDELab::OPBLocalFiniteElementMap<DF,RF,degreeT,dim,Dune::GeometryType::simplex> FEMT;
#endif
  enum { nT=Dune::PB::PkSize<degreeT,dim>::value};  // # temperature degrees of freedom per element

  FEMT femT;

  // Constructor makes objects
  TaylorHood_32_Triangle (const GridView& gv)
    : femu(gv), femp(gv)
  {}
};


template<typename GridView>
class TaylorHood_21_Quadrilateral
{
public:
  // define constants
  enum {dim=GridView::dimension}; // the dimension
  enum {degreeu=2}; // velocity polynomial degree
  enum {degreep=1}; // pressure polynomial degree
  enum {degreeT=2}; // temperature polynomial degree

  enum { mvol=16}; // number of quadrature points on element 
  enum { mbnd=3};  // number of quadrature points on face
  enum { faces=2*dim};  // number of faces

  enum { nu=Dune::QkStuff::QkSize<degreeu,dim>::value};  // # velocity degrees of freedom
  enum { np=Dune::QkStuff::QkSize<degreep,dim>::value};  // # pressure degrees of freedom

  // some types
  typedef GridView GV;        // export also the grid view type
  typedef typename GV::Grid::ctype DF; // type for coordinates
  typedef double RF;          // representation type for solution values

  // finite element map types
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,RF,degreeu> FEMU; // finite element map for velocity
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,RF,degreep> FEMP; // finite element map for pressure

  // finite element map objects
  FEMU femu;
  FEMP femp;

  // DG space for transport
  typedef Dune::PDELab::QkDGLocalFiniteElementMap<DF,RF,degreeT,dim,Dune::PDELab::QkDGBasisPolynomial::legendre> FEMT;
  enum { nT=Dune::QkStuff::QkSize<degreeT,dim>::value};  // # temperature degrees of freedom per element

  FEMT femT;

  // Constructor makes objects
  TaylorHood_21_Quadrilateral (const GridView& gv)
    : femu(gv), femp(gv)
  {}
};

template<typename GridView>
class TaylorHood_21_Hexahedron
{
public:
  // define constants
  enum {dim=GridView::dimension}; // the dimension
  enum {degreeu=2}; // velocity polynomial degree
  enum {degreep=1}; // pressure polynomial degree
  enum {degreeT=2}; // temperature polynomial degree

  enum { mvol=64}; // number of quadrature points on element 
  enum { mbnd=9};  // number of quadrature points on face
  enum { faces=2*dim};  // number of faces

  enum { nu=Dune::QkStuff::QkSize<degreeu,dim>::value};  // # velocity degrees of freedom
  enum { np=Dune::QkStuff::QkSize<degreep,dim>::value};  // # pressure degrees of freedom

  // some types
  typedef GridView GV;        // export also the grid view type
  typedef typename GV::Grid::ctype DF; // type for coordinates
  typedef double RF;          // representation type for solution values

  // finite element map types
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,RF,degreeu> FEMU; // finite element map for velocity
  typedef Dune::PDELab::QkLocalFiniteElementMap<GV,DF,RF,degreep> FEMP; // finite element map for pressure

  // finite element map objects
  FEMU femu;
  FEMP femp;

  // DG space for transport
  typedef Dune::PDELab::QkDGLocalFiniteElementMap<DF,RF,degreeT,dim,Dune::PDELab::QkDGBasisPolynomial::legendre> FEMT;
  enum { nT=Dune::QkStuff::QkSize<degreeT,dim>::value};  // # temperature degrees of freedom per element

  FEMT femT;

  // Constructor makes objects
  TaylorHood_21_Hexahedron (const GridView& gv)
    : femu(gv), femp(gv)
  {}
};


template<typename GridView>
class TaylorHood_21_Tetrahedron
{
public:
  // define constants
  enum {dim=GridView::dimension}; // the dimension
  enum {degreeu=2}; // velocity polynomial degree
  enum {degreep=1}; // pressure polynomial degree
  enum {degreeT=2}; // temperature polynomial degree

  enum { mvol=60}; // number of quadrature points on element 
  enum { mbnd=7};  // number of quadrature points on face
  enum { faces=dim+1};  // number of faces

  enum { nu=Dune::PB::PkSize<degreeu,dim>::value};  // # velocity degrees of freedom
  enum { np=Dune::PB::PkSize<degreep,dim>::value};  // # pressure degrees of freedom

  // some types
  typedef GridView GV;        // export also the grid view type
  typedef typename GV::Grid::ctype DF; // type for coordinates
  typedef double RF;          // representation type for solution values

  // finite element map types
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,RF,degreeu> FEMU; // finite element map for velocity
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,RF,degreep> FEMP; // finite element map for pressure

  // finite element map objects
  FEMU femu;
  FEMP femp;

  // DG space for transport
#if HAVE_GMP
  typedef Dune::PDELab::OPBLocalFiniteElementMap<DF,RF,degreeT,dim,Dune::GeometryType::simplex,Dune::GMPField<512>,Dune::PB::BasisType::Pk> FEMT;
#else
  typedef Dune::PDELab::OPBLocalFiniteElementMap<DF,RF,degreeT,dim,Dune::GeometryType::simplex> FEMT;
#endif
  enum { nT=Dune::PB::PkSize<degreeT,dim>::value};  // # temperature degrees of freedom per element

  FEMT femT;

  // Constructor makes objects
  TaylorHood_21_Tetrahedron (const GridView& gv)
    : femu(gv), femp(gv)
  {}
};


template<typename GridView>
class TaylorHood_32_Tetrahedron
{
public:
  // define constants
  enum {dim=GridView::dimension}; // the dimension
  enum {degreeu=3}; // velocity polynomial degree
  enum {degreep=2}; // pressure polynomial degree
  enum {degreeT=3}; // temperature polynomial degree

  enum { mvol=175}; // number of quadrature points on element 
  enum { mbnd=16};  // number of quadrature points on face
  enum { faces=dim+1};  // number of faces

  enum { nu=Dune::PB::PkSize<degreeu,dim>::value};  // # velocity degrees of freedom
  enum { np=Dune::PB::PkSize<degreep,dim>::value};  // # pressure degrees of freedom

  // some types
  typedef GridView GV;        // export also the grid view type
  typedef typename GV::Grid::ctype DF; // type for coordinates
  typedef double RF;          // representation type for solution values

  // finite element map types
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,RF,degreeu> FEMU; // finite element map for velocity
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,DF,RF,degreep> FEMP; // finite element map for pressure

  // finite element map objects
  FEMU femu;
  FEMP femp;

  // DG space for transport
#if HAVE_GMP
  typedef Dune::PDELab::OPBLocalFiniteElementMap<DF,RF,degreeT,dim,Dune::GeometryType::simplex,Dune::GMPField<512>,Dune::PB::BasisType::Pk> FEMT;
#else
  typedef Dune::PDELab::OPBLocalFiniteElementMap<DF,RF,degreeT,dim,Dune::GeometryType::simplex> FEMT;
#endif
  enum { nT=Dune::PB::PkSize<degreeT,dim>::value};  // # temperature degrees of freedom per element

  FEMT femT;

  // Constructor makes objects
  TaylorHood_32_Tetrahedron (const GridView& gv)
    : femu(gv), femp(gv)
  {}
};
