// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_TAYLORHOODNAVIERSTOKES_HH
#define DUNE_PDELAB_LOCALOPERATOR_TAYLORHOODNAVIERSTOKES_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>

#include<dune/geometry/type.hh>
#include<dune/geometry/quadraturerules.hh>

#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/pdelab/gridfunctionspace/localvector.hh>

#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/l2.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    /* TODO List

       1 mass term: make a completely new operator that copies all the setup in the constructor
         but leaves out the gradient stuff and only implements the alpha volume term
         Use numerical Jacobians first. DONE 18/09/27.
       2 Test problem: start with a simple driven cavity as is usual!
         That can be done without alpha boundary, OK DONE on 18/09/27 It runs
       3 Implement the cylinder problem with do nothing.
       4 Implement directional do-nothing and Navier slip boundary conditions
         - boundary condition type enum
         - boundary condition type function (independent of time)
         - alpha boundary: both, Navier-slip and do nothing bc depend on solution
       5 Implement Jacobians. Do that only after testing the other stuff.
         Suggestion: start with the mass operator
       6 Quadrilateral meshes in Gmsh?
     */

    /** \brief Class to define the boundary condition types
     */
    struct NavierStokesBoundaryCondition
    {
      enum Type { noslip=1, slip=2, donothing=-1}; // BC requiring constraints must be >0 if
    };

    /** \brief A local operator for the Navier-Stokes equations.

        This class implements a local operator for conforming finite
        element discretizations of the Navier-Stokes equations with
        TaylorHood basis. It features i) directional do-nothing boundary
        conditions, ii) Navier-slip boundary condition (for axi-parallel boundary)
        and iii) grad-div stabilization.

        \f{align*}{
        \partial_t u + \nabla\cdot(u u^T)  &= \nabla \cdot \mathbb{S} + f &&\text{in $\Omega\times\Sigma$},\ \
        \mathbb{S} &= \nu \left( \nabla u + (\nabla u)^T \right) - pI,\ \
        \nabla\cdot u &= 0 .
        \f}

        with boundary conditions

        \f{align*}{
        u &= g                             &&\text{on $\Gamma_D$}, &&\text{(Dirichlet)}\ \
        \mathbb{S} n &= \frac12 (u\cdot n)_- u  &&\text{on $\Gamma_O$}, &&\text{(directional do-nothing)}\ \
        u\cdot n &= 0, \ n\times [\mathbb{S} n] \times n + \beta u = 0 &&\text{on $\Gamma_S$}. &&\text{(Navier slip)}
        \f}
    */

    template<typename FLOCAL, typename BCType, typename FEMU, typename FEMP, int degreeu, int degreep, int nu, int np, int mvol, int mbnd, int faces>
    class ConformingNavierStokesLOP :
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public Dune::PDELab::NumericalJacobianVolume<ConformingNavierStokesLOP<FLOCAL,BCType,FEMU,FEMP,degreeu,degreep,nu,np,mvol,mbnd,faces> >,
      public Dune::PDELab::NumericalJacobianApplyVolume<ConformingNavierStokesLOP<FLOCAL,BCType,FEMU,FEMP,degreeu,degreep,nu,np,mvol,mbnd,faces> >,
      public Dune::PDELab::NumericalJacobianBoundary<ConformingNavierStokesLOP<FLOCAL,BCType,FEMU,FEMP,degreeu,degreep,nu,np,mvol,mbnd,faces> >,
      public Dune::PDELab::NumericalJacobianApplyBoundary<ConformingNavierStokesLOP<FLOCAL,BCType,FEMU,FEMP,degreeu,degreep,nu,np,mvol,mbnd,faces> >,
      public InstationaryLocalOperatorDefaultMethods<typename FEMU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType>
    {
      // define useful types
      typedef typename FEMU::Traits::FiniteElementType FiniteElementTypeU;
      typedef typename FiniteElementTypeU::Traits::LocalBasisType LocalBasisTypeU;
      typedef typename FEMP::Traits::FiniteElementType FiniteElementTypeP;
      typedef typename FiniteElementTypeP::Traits::LocalBasisType LocalBasisTypeP;
      typedef typename LocalBasisTypeU::Traits::DomainType DomainType;
      typedef typename LocalBasisTypeU::Traits::DomainFieldType DF;
      typedef typename LocalBasisTypeU::Traits::RangeType RangeType;
      typedef typename LocalBasisTypeU::Traits::RangeFieldType RF;
      typedef typename LocalBasisTypeU::Traits::JacobianType JacobianType;

      enum {dim=LocalBasisTypeU::Traits::dimDomain};

      RF time = 0.0; // guess what

      // quadrature rules
      Dune::FieldVector<DF,dim> qpvol[mvol];   // quadrature points on volume
      RF wvol[mvol];                           // quadrature weight on refelem
      Dune::FieldVector<DF,dim-1> qpbnd[mbnd]; // quadrature points on boundary
      RF wbnd[mbnd];                           // quadrature weight on refelem

      // evaluations of basis functions on the reference element at quadrature points
      RF phihatvol[nu][mvol];          // velocity
      RF psihatvol[np][mvol];          // pressure
      RF Ghat[dim][nu][mvol];          // gradients of velocity basis functions (volume only)
      RF phihatbnd[faces][nu][mbnd];   // velocity

      // coefficient functions
      const FLOCAL& flocal;
      const BCType& bctype;

      // physical and numerical parameters
      RF viscosity; // kinematic viscosity $\mu/\rho$
      RF gamma;     // parameter in grad div stabilization
      RF beta;      // parameter in Navier slip condition

    public:

      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };
      enum { doAlphaBoundary = true };
      enum { doLambdaVolume = true };

      // constructor
      template<typename GV>
      ConformingNavierStokesLOP (const FLOCAL& flocal_, const BCType& bctype_, const FEMU& femu, const FEMP& femp, const GV& gv, RF viscosity_, RF gamma_, RF beta_)
        : flocal(flocal_), bctype(bctype_), viscosity(viscosity_), gamma(gamma_), beta(beta_)
      {
        // The idea is to do all computations on the reference element only once here in the constructor.
        // This implies we assume all elements are of the same type, e.g. simplices,
        // and we use the first element of the mesh as a template.
        // Moreover, we want to use compile-time known loop bounds. This is accomplished
        // by checking that the run-time given basius coomplis with the compile-time
        // given numbers.
        
        // get finite element basis using given rpresentative element
        auto felu = femu.find(*gv.template begin<0>());
        auto felp = femp.find(*gv.template begin<0>());

        // check size of the supplied basis
        if (felu.localBasis().size()!=nu) {
          std::cout << "Basis size mismatch for velocity!" << std::endl;
          exit(1);
        }
        if (felp.localBasis().size()!=np) {
          std::cout << "Basis size mismatch for pressure!" << std::endl;
          exit(1);
        }

        // find quadrature rules with the given number of points and
        // evaluate quadrature rules
        Dune::GeometryType gt = felu.type();
        //std::cout << "New Navier-Stokes LOP on " << gt << std::endl;
        int ordervol=-1; // volume rule
        // for (int order=1; order<=20; order++)
        //   {
        //     const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,order);
        //     std::cout << " order " << order << "qpoints " << rule.size() << std::endl;
        //   }
        // exit(1);
        for (int order=20; order>=1; order--)
          {
            const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,order);
            if (rule.size()==mvol)
              {
                ordervol = order;
                //std::cout << "order of volume quadrature with " << mvol << " points is " << order << std::endl;
                for (int i=0; i<mvol; i++) {
                  qpvol[i] = rule[i].position();
                  wvol[i] = rule[i].weight();
                }
                break;
              }
          }
        if (ordervol<0) {
          std::cout << "Could not find volume quadruture rule with that many points!" << std::endl;
          exit(1);
        }
        int orderbnd=-1; // boundary rule
        Dune::GeometryType gtbnd;
        if (gt.isSimplex()) gtbnd = Dune::GeometryTypes::simplex(dim-1);
        if (gt.isCube()) gtbnd = Dune::GeometryTypes::cube(dim-1);
        for (int order=20; order>=1; order--)
          {
            const Dune::QuadratureRule<RF,dim-1>& rule = Dune::QuadratureRules<RF,dim-1>::rule(gtbnd,order);
            if (rule.size()==mbnd)
              {
                orderbnd = order;
                //std::cout << "order of boundary quadrature with " << mbnd << " points is " << order << std::endl;
                for (int i=0; i<mbnd; i++) {
                  qpbnd[i] = rule[i].position();
                  wbnd[i] = rule[i].weight();
                }
                break;
              }
          }
        if (orderbnd<0) {
          std::cout << "Could not find boundary quadruture rule with that many points!" << std::endl;
          exit(1);
        }

        // evaluate basis functions on reference element for velocity and pressure
        for (int k=0; k<mvol; k++)
          {
            std::vector<RangeType> phi(nu);
            felu.localBasis().evaluateFunction(qpvol[k],phi);
            for (int i=0; i<nu; i++) phihatvol[i][k] = phi[i];
          }
        for (int k=0; k<mvol; k++)
          {
            std::vector<RangeType> psi(np);
            felp.localBasis().evaluateFunction(qpvol[k],psi);
            for (int i=0; i<np; i++) psihatvol[i][k] = psi[i];
          }
        for (const auto& is : intersections(gv,*gv.template begin<0>()))
          {
            auto fgeo_self = is.geometryInInside();
            int face = is.indexInInside();
            //std::cout << "face " << face << std::endl;
            for (int k=0; k<mbnd; k++) // loop over face quadrature points
              {
                auto qp = fgeo_self.global(qpbnd[k]);
                std::vector<RangeType> phi(nu);
                felu.localBasis().evaluateFunction(qp,phi);
                for (int i=0; i<nu; i++) phihatbnd[face][i][k] = phi[i];
              }
          }
        
        // evaluate gradients of basis functions on reference element
        std::vector<JacobianType> js(nu);
        for (int k=0; k<mvol; k++)
          {
            felu.localBasis().evaluateJacobian(qpvol[k],js);
            for (int i=0; i<nu; i++)
              for (int j=0; j<dim; j++)
                Ghat[j][i][k] = js[i][0][j];
          }
      }

      // set time in parameter class
      void setTime (RF t)
      {
        time = t;
      }

      //! residual contribution of volume source term
      template<typename EG, typename LFSTest, typename Residual>
      void lambda_volume (const EG& eg, const LFSTest& lfstest, Residual& residual) const
      {
        // define types
        using LFSV = typename LFSTest::template Child<0>::Type;   // the velocity node
        using LFSV0 = typename LFSV::template Child<0>::Type;     // first velocity component

        // extract local function spaces
        const auto& lfsv = lfstest.template child<0>();    // velocity node
        
        // evaluate geometry transformation at all quadrature points
        auto geo = eg.geometry();
        RF factor[mvol];       // quadrature weight times determinant
        if (geo.type().isSimplex())
          {
            auto integrationelement = geo.integrationElement(qpvol[0]);
            for (size_t alpha=0; alpha<mvol; ++alpha)
                factor[alpha] = integrationelement*wvol[alpha];
          }
        if (geo.type().isCube())
          for (size_t alpha=0; alpha<mvol; ++alpha)
            {
              auto integrationelement = geo.integrationElement(qpvol[alpha]);
              factor[alpha] = integrationelement*wvol[alpha];
            }

        // Now evaluate velocity
        // Here it turns out to be good idea to have the quadrature loop inside!
        
        // evaluate force at all quadrature points
        RF FT[dim][mvol] = {{0.0}};
        for (size_t alpha=0; alpha<mvol; ++alpha)
          {
            auto f = flocal(eg.entity(),qpvol[alpha]);
            for (size_t i=0; i<dim; ++i) FT[i][alpha] = f[i];
          }

        // Now come the residual contributions term by term.
        RF RT[dim][nu] = {{0.0}}; // residual contribution of this element to velocity test function [i,j]

        // mass residual
        for (size_t alpha=0; alpha<mvol; ++alpha) // reduction over alpha
          for (size_t i=0; i<dim; ++i) 
            for (size_t j=0; j<nu; ++j)
              RT[i][j] += FT[i][alpha]*phihatvol[j][alpha]*factor[alpha];

        // accumulate residual contributions for this element
        for (size_t i=0; i<dim; ++i)
          {
            const auto& lfsv_i = lfsv.child(i);
            for (size_t j=0; j<nu; ++j)
              residual.accumulate(lfsv_i,j,-RT[i][j]);
          }
      }
      
      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSTrial, typename X, typename LFSTest, typename Residual>
      void alpha_volume (const EG& eg, const LFSTrial& lfstrial, const X& x, const LFSTest& lfstest, Residual& residual) const
      {
        // extract local function spaces
        const auto& lfsu = lfstrial.template child<0>();    // velocity node
        const auto& lfsp = lfstrial.template child<1>();  // pressure

        // extract degrees of freedom
        RF ZT[dim][nu]; // velocity degrees of freedom
        for (size_t i=0; i<dim; ++i)
          {
            const auto& lfsu_i = lfsu.child(i);
            for (size_t j=0; j<nu; ++j) ZT[i][j] = x(lfsu_i,j);
          }
        RF YT[np]; // pressure degrees of freedom
        for (size_t j=0; j<np; ++j) YT[j] = x(lfsp,j);
        
        // evaluate geometry transformation at all quadrature points
        auto geo = eg.geometry();
        RF JT[dim][dim][mvol]; // Jacobian inverse transposed for gradient transformations
        RF factor[mvol];       // quadrature weight times determinant
        if (geo.type().isSimplex())
          {
            const auto S = geo.jacobianInverseTransposed(qpvol[0]);
            auto integrationelement = geo.integrationElement(qpvol[0]);
            for (size_t alpha=0; alpha<mvol; ++alpha)
                factor[alpha] = integrationelement*wvol[alpha];
            for (size_t i=0; i<dim; ++i)
              for (size_t j=0; j<dim; ++j)
                for (size_t alpha=0; alpha<mvol; ++alpha)
                  JT[i][j][alpha] = S[i][j];
          }
        if (geo.type().isCube())
          for (size_t alpha=0; alpha<mvol; ++alpha)
            {
              const auto S = geo.jacobianInverseTransposed(qpvol[alpha]);
              auto integrationelement = geo.integrationElement(qpvol[alpha]);
              factor[alpha] = integrationelement*wvol[alpha];
              for (size_t i=0; i<dim; ++i)
                for (size_t j=0; j<dim; ++j)
                  JT[i][j][alpha] = S[i][j];
            }

        // compute gradient of velocity basis functions at all quadrature points
        RF GT[dim][nu][mvol] = {{{0.0}}};
        for (size_t i=0; i<dim; ++i)
          for (size_t j=0; j<nu; ++j)
            for (size_t k=0; k<dim; ++k) // reduction over k
              for (size_t alpha=0; alpha<mvol; ++alpha)
                GT[i][j][alpha] += JT[i][k][alpha]*Ghat[k][j][alpha];

        // Now evaluate velocity, velocity gradient and pressure at quadrature points
        // Here it turns out to be good idea to have the quadrature loop inside!
        
        // compute velocity at all quadrature points
        RF UT[dim][mvol] = {{0.0}};
        for (size_t i=0; i<dim; ++i)
          for (size_t j=0; j<nu; ++j) // reduction over j
            for (size_t alpha=0; alpha<mvol; ++alpha)
              UT[i][alpha] += ZT[i][j]*phihatvol[j][alpha];

        // compute velocity gradient at all quadrature points
        RF GradUT[dim][dim][mvol] = {{{0.0}}};
        for (size_t i=0; i<dim; ++i)
          for (size_t j=0; j<dim; ++j)
            for (size_t k=0; k<nu; ++k) // reduction over k
              for (size_t alpha=0; alpha<mvol; ++alpha)
                GradUT[i][j][alpha] += ZT[i][k]*GT[j][k][alpha];

        // compute divergence of velocity at all quadrature points
        RF DivUT[mvol] = {0.0};
        for (size_t i=0; i<dim; ++i)  // double reduction over
          for (size_t j=0; j<nu; ++j) // i and j !
            for (size_t alpha=0; alpha<mvol; ++alpha)
              DivUT[alpha] += ZT[i][j]*GT[i][j][alpha];

        // compute pressure at all quadratue points
        RF PT[mvol] = {0.0};
        for (size_t j=0; j<np; ++j) // reduction over j
          for (size_t alpha=0; alpha<mvol; ++alpha)
            PT[alpha] += YT[j]*psihatvol[j][alpha];

        // Now come the residual contributions term by term.
        // These are ultimately reductions over the quadrature points but
        // possibly involve some intermediate results beforehand
        RF RUT[dim][nu] = {{0.0}}; // residual contribution of this element to velocity test function [i,j]
        RF RPT[np]  = {0.0};       // residual contribution of this element to pressure test function [j]

        // // Test purpose: mass on velocity
        // for (size_t alpha=0; alpha<mvol; ++alpha) // reduction over alpha
        //   for (size_t i=0; i<dim; ++i) 
        //     for (size_t j=0; j<nu; ++j)
        //       RUT[i][j] += UT[i][alpha]*phihatvol[j][alpha]*factor[alpha];
        // // Test purpose: Mass on pressure
        // for (size_t alpha=0; alpha<mvol; ++alpha) // reduction over alpha
        //   for (size_t j=0; j<np; ++j) // pressure test functions here!
        //     RPT[j] += PT[alpha]*psihatvol[j][alpha]*factor[alpha];

        // // Analysis:
        // std::cout << "velocity coefficients " << dim << " times " << nu << std::endl;
        // for (size_t i=0; i<dim; i++) {
        //   std::cout << i << ": ";
        //   for (size_t j=0; j<nu; ++j) std::cout << ZT[i][j] << " ";
        //   std::cout << std::endl;
        // }
        // std::cout << "pressure coefficients " << np << std::endl;
        // for (size_t j=0; j<np; ++j) std::cout << YT[j] << " ";
        // std::cout << std::endl;
        // std::cout << "velocity at quadrature points " << dim << " times " << mvol << std::endl;
        // for (size_t i=0; i<dim; i++) {
        //   std::cout << i << ": ";
        //   for (size_t alpha=0; alpha<mvol; ++alpha) std::cout << UT[i][alpha] << " ";
        //   std::cout << std::endl;
        // }
        // std::cout << "pressure at quadrature points " << mvol << std::endl;
        // for (size_t alpha=0; alpha<mvol; ++alpha) std::cout << PT[alpha] << " ";
        // std::cout << std::endl;
        // std::cout << "velocity RESIDUAL " << dim << " times " << nu << std::endl;
        // for (size_t i=0; i<dim; i++) {
        //   std::cout << i << ": ";
        //   for (size_t j=0; j<nu; ++j) std::cout << RUT[i][j] << " ";
        //   std::cout << std::endl;
        // }
        // std::cout << "pressure RESIDUAL " << np << std::endl;
        // for (size_t j=0; j<np; ++j) std::cout << RPT[j] << " ";
        // std::cout << std::endl;
        
        // Stress term
        {
          RF intermediate1[dim][nu][mvol] = {{{0.0}}};
          for (size_t i=0; i<dim; ++i) 
            for (size_t j=0; j<nu; ++j)
              for (size_t k=0; k<dim; ++k) // reduction over k
                for (size_t alpha=0; alpha<mvol; ++alpha)
                  intermediate1[i][j][alpha] += (GradUT[i][k][alpha]+GradUT[k][i][alpha])*GT[k][j][alpha];
          RF viscosity_factor[mvol];
          for (size_t alpha=0; alpha<mvol; ++alpha) viscosity_factor[alpha] = viscosity*factor[alpha];
          for (size_t i=0; i<dim; ++i) 
            for (size_t alpha=0; alpha<mvol; ++alpha) // reduction over alpha
              for (size_t j=0; j<nu; ++j)
                RUT[i][j] += intermediate1[i][j][alpha]*viscosity_factor[alpha];
        }
        
        // convection term
        {
          RF intermediate1[dim][mvol] = {{0.0}}; // grad u_h * u_h at quadrature points
          for (size_t i=0; i<dim; ++i) 
            for (size_t j=0; j<dim; ++j) // reduction over j
              for (size_t alpha=0; alpha<mvol; ++alpha)
                intermediate1[i][alpha] += GradUT[i][j][alpha]*UT[j][alpha]; //grad u_h * u_h
          for (size_t i=0; i<dim; ++i) 
            for (size_t alpha=0; alpha<mvol; ++alpha)
              intermediate1[i][alpha] += UT[i][alpha]*DivUT[alpha]; // add u_h * div u_h
          // now we dont want to store the outer product ...
          for (size_t alpha=0; alpha<mvol; ++alpha) // reduction over alpha
            for (size_t i=0; i<dim; ++i) 
              for (size_t j=0; j<nu; ++j)
                RUT[i][j] += intermediate1[i][alpha]*phihatvol[j][alpha]*factor[alpha];
        }

        // pressure term
        for (size_t alpha=0; alpha<mvol; ++alpha) // reduction over alpha
          for (size_t i=0; i<dim; ++i) 
            for (size_t j=0; j<nu; ++j)
              RUT[i][j] -= PT[alpha]*GT[i][j][alpha]*factor[alpha];

        // grad-div stabilization
        if (gamma>0.0)
          {
            for (size_t alpha=0; alpha<mvol; ++alpha) // reduction over alpha
              for (size_t i=0; i<dim; ++i) 
                for (size_t j=0; j<nu; ++j)
                  RUT[i][j] += gamma*DivUT[alpha]*GT[i][j][alpha]*factor[alpha];
          }

        // divergence term involving pressure test functions
        for (size_t alpha=0; alpha<mvol; ++alpha) // reduction over alpha
          for (size_t j=0; j<np; ++j) // pressure test functions here!
            RPT[j] -= DivUT[alpha]*psihatvol[j][alpha]*factor[alpha];

        // accumulate residual contributions for this element
        for (size_t i=0; i<dim; ++i)
          {
            const auto& lfsu_i = lfsu.child(i);
            for (size_t j=0; j<nu; ++j)
              residual.accumulate(lfsu_i,j, RUT[i][j]);
          }
        for (size_t j=0; j<np; ++j)
          residual.accumulate(lfsp,j, RPT[j]);
      }

      // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
      template<typename IG, typename LFSTrial, typename X,
               typename LFSTest, typename Residual>
      void alpha_boundary (const IG& ig,
                           const LFSTrial& lfstrial, const X& x,
                           const LFSTest& lfstest, Residual& residual) const
      {
        // check for Dirichlet boundary condition first
        auto bct = bctype(ig.geometry().center());
        if (bct == NavierStokesBoundaryCondition::noslip) return;

        // extract local function spaces
        const auto& lfsu = lfstrial.template child<0>();    // velocity node

        // extract degrees of freedom
        RF ZT[dim][nu]; // velocity degrees of freedom
        for (size_t i=0; i<dim; ++i)
          {
            const auto& lfsu_i = lfsu.child(i);
            for (size_t j=0; j<nu; ++j) ZT[i][j] = x(lfsu_i,j);
          }

        // evaluate geometry
        auto facegeo = ig.geometry(); // map from face refelem to face
        int face = ig.indexInInside(); // number of face

        // quadrature factor. Note: we integrate over face refefrence element
        RF factor[mbnd];       // quadrature weight times determinant
        if (facegeo.type().isSimplex())
          {
            auto integrationelement = facegeo.integrationElement(qpbnd[0]);
            for (size_t alpha=0; alpha<mbnd; ++alpha)
                factor[alpha] = integrationelement*wbnd[alpha];
          }
        if (facegeo.type().isCube())
          for (size_t alpha=0; alpha<mbnd; ++alpha)
            {
              auto integrationelement = facegeo.integrationElement(qpbnd[alpha]);
              factor[alpha] = integrationelement*wbnd[alpha];
            }

        // compute velocity at all quadrature points
        RF UT[dim][mbnd] = {{0.0}};
        for (size_t i=0; i<dim; ++i)
          for (size_t j=0; j<nu; ++j) // reduction over j
            for (size_t alpha=0; alpha<mbnd; ++alpha)
              UT[i][alpha] += ZT[i][j]*phihatbnd[face][j][alpha];
        
        // directional do nothing condition
        if (bct == NavierStokesBoundaryCondition::donothing)
          {
            // compute normal velocity over face at quadrature points
            // assume same normal direction
            RF normalvelocity[mbnd] = {0.0};
            auto normal = ig.centerUnitOuterNormal();
            for (size_t i=0; i<dim; ++i) normal[i] *= 0.5;
            for (size_t i=0; i<dim; ++i)
              for (size_t alpha=0; alpha<mbnd; ++alpha)
                normalvelocity[alpha] += UT[i][alpha]*normal[i];
            for (size_t alpha=0; alpha<mbnd; ++alpha)
              normalvelocity[alpha] = std::min(normalvelocity[alpha],0.0); // chop off positive entries

            // Now the residual contribution
            RF RT[dim][nu] = {{0.0}}; // residual contribution of this element to velocity test function [i,j]
            for (size_t alpha=0; alpha<mbnd; ++alpha) // reduction over alpha
              for (size_t i=0; i<dim; ++i) 
                for (size_t j=0; j<nu; ++j)
                  RT[i][j] += normalvelocity[alpha]*UT[i][alpha]*phihatbnd[face][j][alpha]*factor[alpha];

            // accumulate residual contributions for this element
            for (size_t i=0; i<dim; ++i)
              {
                const auto& lfsu_i = lfsu.child(i);
                for (size_t j=0; j<nu; ++j)
                  residual.accumulate(lfsu_i,j,-RT[i][j]);
              }
          }

        // Navier slip condition
        if (bct == NavierStokesBoundaryCondition::slip)
          {
            // check if we have beta>0
            if (beta<=0.0) return;
            
            // Now the residual contribution
            RF RT[dim][nu] = {{0.0}}; // residual contribution of this element to velocity test function [i,j]
            for (size_t alpha=0; alpha<mbnd; ++alpha) // reduction over alpha
              for (size_t i=0; i<dim; ++i) 
                for (size_t j=0; j<nu; ++j)
                  RT[i][j] += beta*UT[i][alpha]*phihatbnd[face][j][alpha]*factor[alpha];

            // accumulate residual contributions for this element
            for (size_t i=0; i<dim; ++i)
              {
                const auto& lfsu_i = lfsu.child(i);
                for (size_t j=0; j<nu; ++j)
                  residual.accumulate(lfsu_i,j,RT[i][j]);
              }
          }
      }

    };



    /** \brief A local operator for the mass term in the Navier-Stokes equations.

    */

    template<typename FEMU, int degreeu, int nu,int mvol>
    class ConformingNavierStokesMassLOP :
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      // public Dune::PDELab::NumericalJacobianVolume<ConformingNavierStokesMassLOP<FEMU,degreeu,nu,mvol> >,
      public Dune::PDELab::NumericalJacobianApplyVolume<ConformingNavierStokesMassLOP<FEMU,degreeu,nu,mvol> >,
      public InstationaryLocalOperatorDefaultMethods<typename FEMU::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType>
    {
      // define useful types
      typedef typename FEMU::Traits::FiniteElementType FiniteElementTypeU;
      typedef typename FiniteElementTypeU::Traits::LocalBasisType LocalBasisTypeU;
      typedef typename LocalBasisTypeU::Traits::DomainType DomainType;
      typedef typename LocalBasisTypeU::Traits::DomainFieldType DF;
      typedef typename LocalBasisTypeU::Traits::RangeType RangeType;
      typedef typename LocalBasisTypeU::Traits::RangeFieldType RF;
      typedef typename LocalBasisTypeU::Traits::JacobianType JacobianType;

      enum {dim=LocalBasisTypeU::Traits::dimDomain};

      RF time = 0.0; // guess what

      // quadrature rules
      Dune::FieldVector<DF,dim> qpvol[mvol];   // quadrature points on volume
      RF wvol[mvol];                           // quadrature weight on refelem

      // evaluations of basis functions on the reference element at quadrature points
      RF phihatvol[nu][mvol];          // velocity

    public:

      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      // constructor
      template<typename GV>
      ConformingNavierStokesMassLOP (const FEMU& femu, const GV& gv)
      {
        // The idea is to do all computations on the reference element only once here in the constructor.
        // This implies we assume all elements are of the same type, e.g. simplices,
        // and we use the first element of the mesh as a template.
        // Moreover, we want to use compile-time known loop bounds. This is accomplished
        // by checking that the run-time given basius coomplis with the compile-time
        // given numbers.
        
        // get finite element basis using given rpresentative element
        auto felu = femu.find(*gv.template begin<0>());

        // check size of the supplied basis
        if (felu.localBasis().size()!=nu) {
          std::cout << "Basis size mismatch for velocity!" << std::endl;
          exit(1);
        }

        // find quadrature rules with the given number of points and
        // evaluate quadrature rules
        Dune::GeometryType gt = felu.type();
        int ordervol=-1; // volume rule
        for (int order=20; order>=1; order--)
          {
            const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,order);
            if (rule.size()==mvol)
              {
                ordervol = order;
                std::cout << "order of volume quadrature with " << mvol << " points is " << order << std::endl;
                for (int i=0; i<mvol; i++) {
                  qpvol[i] = rule[i].position();
                  wvol[i] = rule[i].weight();
                }
                break;
              }
          }
        if (ordervol<0) {
          std::cout << "Could not find volume quadruture rule with that many points!" << std::endl;
          exit(1);
        }

        // evaluate basis functions on reference element for velocity and pressure
        for (int k=0; k<mvol; k++)
          {
            std::vector<RangeType> phi(nu);
            felu.localBasis().evaluateFunction(qpvol[k],phi);
            for (int i=0; i<nu; i++) phihatvol[i][k] = phi[i];
          }
      }

      // set time in parameter class
      void setTime (RF t)
      {
        time = t;
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSTrial, typename X, typename LFSTest, typename Residual>
      void alpha_volume (const EG& eg, const LFSTrial& lfstrial, const X& x, const LFSTest& lfstest, Residual& residual) const
      {
        // extract local function spaces
        const auto& lfsu = lfstrial.template child<0>();    // velocity node

        // extract degrees of freedom
        RF ZT[dim][nu]; // velocity degrees of freedom
        for (size_t i=0; i<dim; ++i)
          {
            const auto& lfsu_i = lfsu.child(i);
            for (size_t j=0; j<nu; ++j) ZT[i][j] = x(lfsu_i,j);
          }
        
        // evaluate geometry transformation at all quadrature points
        auto geo = eg.geometry();
        RF factor[mvol];       // quadrature weight times determinant
        if (geo.type().isSimplex())
          {
            auto integrationelement = geo.integrationElement(qpvol[0]);
            for (size_t alpha=0; alpha<mvol; ++alpha)
                factor[alpha] = integrationelement*wvol[alpha];
          }
        if (geo.type().isCube())
          for (size_t alpha=0; alpha<mvol; ++alpha)
            {
              auto integrationelement = geo.integrationElement(qpvol[alpha]);
              factor[alpha] = integrationelement*wvol[alpha];
            }

        // Now evaluate velocity
        // Here it turns out to be good idea to have the quadrature loop inside!
        
        // compute velocity at all quadrature points
        RF UT[dim][mvol] = {{0.0}};
        for (size_t i=0; i<dim; ++i)
          for (size_t j=0; j<nu; ++j) // reduction over j
            for (size_t alpha=0; alpha<mvol; ++alpha)
              UT[i][alpha] += ZT[i][j]*phihatvol[j][alpha];

        // Now come the residual contributions term by term.
        // These are ultimately reductions over the quadrature points but
        // possibly involve some intermediate results beforehand
        RF RUT[dim][nu] = {{0.0}}; // residual contribution of this element to velocity test function [i,j]

        // mass residual
        for (size_t alpha=0; alpha<mvol; ++alpha) // reduction over alpha
          for (size_t i=0; i<dim; ++i) 
            for (size_t j=0; j<nu; ++j)
              RUT[i][j] += UT[i][alpha]*phihatvol[j][alpha]*factor[alpha];

        // accumulate residual contributions for this element
        for (size_t i=0; i<dim; ++i)
          {
            const auto& lfsu_i = lfsu.child(i);
            for (size_t j=0; j<nu; ++j)
              residual.accumulate(lfsu_i,j, RUT[i][j]);
          }
      }

      template<typename EG, typename LFSTrial, typename X, typename LFSTest, typename M>
      void jacobian_volume (const EG& eg, const LFSTrial& lfstrial, const X& x, const LFSTest& lfstest,
                            M& mat) const
      {
        // extract local function spaces, assume Galerkin property: lfstrial = lfstest
        const auto& lfsu = lfstrial.template child<0>();    // velocity node

        // evaluate geometry transformation at all quadrature points
        auto geo = eg.geometry();
        RF factor[mvol];       // quadrature weight times determinant
        if (geo.type().isSimplex())
          {
            auto integrationelement = geo.integrationElement(qpvol[0]);
            for (size_t alpha=0; alpha<mvol; ++alpha)
                factor[alpha] = integrationelement*wvol[alpha];
          }
        if (geo.type().isCube())
          for (size_t alpha=0; alpha<mvol; ++alpha)
            {
              auto integrationelement = geo.integrationElement(qpvol[alpha]);
              factor[alpha] = integrationelement*wvol[alpha];
            }

        // Compute the local mass matrix for one component
        // This mass matrix block for one component is an outer product
        RF MT[nu][nu] = {{0.0}};
        for (size_t alpha=0; alpha<mvol; ++alpha) // reduction over alpha
          for (size_t i=0; i<nu; ++i) 
            for (size_t j=0; j<nu; ++j)
              MT[i][j] += phihatvol[i][alpha]*phihatvol[j][alpha]*factor[alpha];
        
        // accumulate jacobian contributions for this element
        for (size_t c=0; c<dim; ++c)
          {
            const auto& lfsu_c = lfsu.child(c);
            for (size_t i=0; i<nu; ++i)
              for (size_t j=0; j<nu; ++j)
                mat.accumulate(lfsu_c,i,lfsu_c,j,MT[i][j]);
          }
      }
        
    };


    //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif
