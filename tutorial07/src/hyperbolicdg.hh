// -*- tab-width: 2; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_LOCALOPERATOR_HYPERBOLICDG_HH
#define DUNE_PDELAB_LOCALOPERATOR_HYPERBOLICDG_HH

#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/idefault.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>


namespace Dune {
  namespace PDELab {
    /** Spatial local operator for discontinuous Galerkin method
        for the system of hyperbolic conservation laws:

        \nabla \cdot \{ F(u) \}  = 0 in \Omega

        Where u = (u1,...,um) is the solution with m components

        - Assumes that the local function space is a power space
        with m identical components.
        - Assumes Galerkin method, i.e. U=V

        \tparam P parameter class
        \tparam FEM Finite Element Map needed to select the cache
    */
    template<typename NUMFLUX, typename FEM>
    class DGHyperbolicSpatialOperator :
      public NumericalJacobianApplyVolume<DGHyperbolicSpatialOperator<NUMFLUX,FEM> >,
      public NumericalJacobianVolume<DGHyperbolicSpatialOperator<NUMFLUX,FEM> >,
      public NumericalJacobianApplySkeleton<DGHyperbolicSpatialOperator<NUMFLUX,FEM> >,
      public NumericalJacobianSkeleton<DGHyperbolicSpatialOperator<NUMFLUX,FEM> >,
      public NumericalJacobianApplyBoundary<DGHyperbolicSpatialOperator<NUMFLUX,FEM> >,
      public NumericalJacobianBoundary<DGHyperbolicSpatialOperator<NUMFLUX,FEM> >,
      public FullSkeletonPattern,
      public FullVolumePattern,
      public LocalOperatorDefaultFlags,
      public InstationaryLocalOperatorDefaultMethods<typename NUMFLUX::RF>
    {

        static constexpr int dim = NUMFLUX::dim;
        static constexpr int m = NUMFLUX::m; //number of components

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };
      enum { doPatternSkeleton = true };

      // residual assembly flags
      enum { doAlphaVolume  = true };
      enum { doAlphaSkeleton  = true };
      enum { doAlphaBoundary  = true };
      enum { doLambdaVolume  = true };

      // ! constructor
      DGHyperbolicSpatialOperator (NUMFLUX& numflux_, int overintegration_=0)
        : numflux(numflux_), overintegration(overintegration_), cache(20)
      {
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // Get types
        using namespace Indices;
        using RF = typename NUMFLUX::RF; // type for computations

        // get local function space that is identical for all components
        const auto& dgspace = child(lfsv,_0);

        // Get geometry
        auto geo = eg.geometry();

        // Transformation
        typename EG::Geometry::JacobianInverseTransposed jac;

        // Initialize vectors outside for loop
        Dune::FieldVector<RF,m> u(0.0);
        std::vector<Dune::FieldVector<RF,dim> > gradphi(dgspace.size());

        // loop over quadrature points
        const int order = dgspace.finiteElement().localBasis().order();
        const int intorder = overintegration+2*order;
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            const auto qp = ip.position();
            // evaluate basis functions
            auto& phi = cache[order].evaluateFunction(qp,dgspace.finiteElement().localBasis());

            // evaluate u
            u = 0.0;
            for (size_t k=0; k<m; k++) // for all components
              for (size_t j=0; j<dgspace.size(); j++) // for all basis functions
                u[k] += x(lfsv.child(k),j)*phi[j];
            // std::cout << "  u at " << qp << " : " << u << std::endl;

            // evaluate gradient of basis functions (we assume Galerkin method lfsu=lfsv)
            auto& js = cache[order].evaluateJacobian(qp,dgspace.finiteElement().localBasis());

            // compute global gradients
            jac = geo.jacobianInverseTransposed(qp);
            for (size_t i=0; i<dgspace.size(); i++)
              jac.mv(js[i][0],gradphi[i]);

            /// tex: fluxint
            Dune::FieldMatrix<RF,m,dim> F;
            numflux.model().flux(eg,qp,u,F);

            // integrate
            auto factor = ip.weight() * geo.integrationElement(qp);
            for (size_t k=0; k<dgspace.size(); k++)
              {
                // - F(u) \grad phi
                for (size_t i=0; i<m; i++)
                  for (size_t j=0; j<dim; j++)
                    r.accumulate(lfsv.child(i),k,-F[i][j]*gradphi[k][j]*factor);
              }
            /// tex: fluxint
          }
      }

      // skeleton integral depending on test and ansatz functions
      // each face is only visited ONCE!
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_skeleton (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                           R& r_s, R& r_n) const
      {
        // Get types
        using namespace Indices;
        using RF = typename NUMFLUX::RF;//range field

        // Get local function space that is identical for all components
        const auto& dgspace_s = child(lfsv_s,_0);
        const auto& dgspace_n = child(lfsv_n,_0);

        // References to inside and outside cells
        const auto& cell_inside = ig.inside();
        const auto& cell_outside = ig.outside();

        // Get geometries
        auto geo = ig.geometry();
        auto geo_in_inside = ig.geometryInInside();
        auto geo_in_outside = ig.geometryInOutside();

        // Initialize vectors outside for loop
        Dune::FieldVector<RF,m> u_s(0.0);
        Dune::FieldVector<RF,m> u_n(0.0);
        Dune::FieldVector<RF,m> f(0.0);

        // Loop over quadrature points
        const int order_s = dgspace_s.finiteElement().localBasis().order();
        const int order_n = dgspace_n.finiteElement().localBasis().order();
        const int intorder = overintegration+1+2*std::max(order_s,order_n);

        /// tex: skeleton
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            const auto qp = ip.position();
            // Position of quadrature point in local coordinates of elements
            auto iplocal_s = geo_in_inside.global(qp);
            auto iplocal_n = geo_in_outside.global(qp);

            // Evaluate basis functions
            auto& phi_s = cache[order_s].evaluateFunction(
              iplocal_s,dgspace_s.finiteElement().localBasis());
            auto& phi_n = cache[order_n].evaluateFunction(
              iplocal_n,dgspace_n.finiteElement().localBasis());

            // Evaluate u from inside and outside
            u_s = 0.0;
            for (size_t i=0; i<m; i++) // for all components
              for (size_t k=0; k<dgspace_s.size(); k++)
                u_s[i] += x_s(lfsv_s.child(i),k)*phi_s[k];
            u_n = 0.0;
            for (size_t i=0; i<m; i++) // for all components
              for (size_t k=0; k<dgspace_n.size(); k++)
                u_n[i] += x_n(lfsv_n.child(i),k)*phi_n[k];

            // Compute numerical flux at  the integration point
            numflux.numericalFlux(cell_inside, iplocal_s,
                                  cell_outside, iplocal_n,
                                  ig.unitOuterNormal(qp),u_s,u_n,f);

            // Integrate
            auto factor = ip.weight() * geo.integrationElement(qp);
            // loop over all vector-valued basis functions
            for (size_t k=0; k<dgspace_s.size(); k++)
              for (size_t i=0; i<m; i++) // loop over all components
                r_s.accumulate(lfsv_s.child(i),k, f[i]*phi_s[k]*factor);
            // loop over all vector-valued basis functions
            for (size_t k=0; k<dgspace_n.size(); k++)
              for (size_t i=0; i<m; i++) // loop over all components
                r_n.accumulate(lfsv_n.child(i),k, - f[i]*phi_n[k]*factor);
          }
        /// tex: skeleton
      }

      // Skeleton integral depending on test and ansatz functions
      template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_boundary (const IG& ig,
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s) const
      {
        // Get types
        using namespace Indices;
        using RF = typename NUMFLUX::RF; // type for computations

        // Get local function space that is identical for all components
        const auto& dgspace_s = child(lfsv_s,_0);

        // Reference to inside cell
        const auto& cell_inside = ig.inside();

        // Get geometries
        auto geo = ig.geometry();
        auto geo_in_inside = ig.geometryInInside();

        // Initialize vectors outside for loop
        Dune::FieldVector<RF,m> u_s(0.0);
        Dune::FieldVector<RF,m> f(0.0);

        // Loop over quadrature points
        const int order_s = dgspace_s.finiteElement().localBasis().order();
        const int intorder = overintegration+1+2*order_s;
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            const auto qp = ip.position();
            // Position of quadrature point in local coordinates of elements
            auto iplocal_s = geo_in_inside.global(qp);
            // Evaluate basis functions
            auto& phi_s = cache[order_s].evaluateFunction(iplocal_s,dgspace_s.finiteElement().localBasis());

            // Evaluate u from inside and outside
            u_s = 0.0;
            for (size_t i=0; i<m; i++) // for all components
              for (size_t k=0; k<dgspace_s.size(); k++) // for all basis functions
                u_s[i] += x_s(lfsv_s.child(i),k)*phi_s[k];
            // std::cout << "  u_s " << u_s << std::endl;

            /// tex: boundary
            // Evaluate boundary condition
            Dune::FieldVector<RF,m> u_n(
              numflux.model().problem.g(ig.intersection(),qp,u_s));

            // Compute numerical flux at integration point
            numflux.numericalFlux(cell_inside, iplocal_s,
                                  cell_inside, iplocal_s,
                                  ig.unitOuterNormal(qp),u_s,u_n,f);

            // Integrate
            auto factor = ip.weight() * geo.integrationElement(qp);
            // loop over all vector-valued (!) basis functions
            // (with identical components)
            for (size_t k=0; k<dgspace_s.size(); k++)
              for (size_t i=0; i<m; i++) // loop over all components
                r_s.accumulate(lfsv_s.child(i),k, f[i]*phi_s[k]*factor);

            /// tex: boundary
          }
        // std::cout << "  residual_s: ";
        // for (size_t i=0; i<r_s.size(); i++) std::cout << r_s[i] << " ";
        // std::cout << std::endl;
      }

      // Volume integral depending only on test functions
      template<typename EG, typename LFSV, typename R>
      void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
      {
        // Get types
        using namespace Indices;
        using DGSpace = TypeTree::Child<LFSV,0>;

        // Get local function space that is identical for all components
        const DGSpace& dgspace = child(lfsv,_0);

        // Reference to cell
        const auto& cell = eg.entity();

        // Get geometries
        auto geo = eg.geometry();

        // Loop over quadrature points
        const int order_s = dgspace.finiteElement().localBasis().order();
        const int intorder = overintegration+2*order_s;
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            const auto qp = ip.position();
            // Evaluate right hand side q
            auto q(numflux.model().problem.q(cell,qp));

            // Evaluate basis functions
            auto& phi = cache[order_s].evaluateFunction(qp,dgspace.finiteElement().localBasis());

            // Integrate
            auto factor = ip.weight() * geo.integrationElement(qp);
            for (size_t k=0; k<m; k++) // for all components
              for (size_t i=0; i<dgspace.size(); i++) // for all test functions of this component
                r.accumulate(lfsv.child(k),i, - q[k]*phi[i]*factor);
          }
      }

    private:
      NUMFLUX& numflux;
      int overintegration;
      using LocalBasisType = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
      using Cache = Dune::PDELab::LocalBasisCache<LocalBasisType>;
      std::vector<Cache> cache;
    };


    /** a local operator for the mass operator of a vector valued lfs (L_2 integral)
     *
     * \f{align*}{
     \int_\Omega uv dx
     * \f}
     */
    template<typename NUMFLUX, typename FEM>
    class DGHyperbolicTemporalOperator :
      public NumericalJacobianApplyVolume<DGHyperbolicTemporalOperator<NUMFLUX,FEM> >,
        public LocalOperatorDefaultFlags,
        public InstationaryLocalOperatorDefaultMethods<typename NUMFLUX::RF>
    {

      static constexpr int dim = NUMFLUX::dim;
      static constexpr int m = NUMFLUX::m;

    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      DGHyperbolicTemporalOperator (NUMFLUX& numflux_, int overintegration_=0)
        : numflux(numflux_), overintegration(overintegration_), cache(20)
      {}

      // define sparsity pattern of operator representation
      template<typename LFSU, typename LFSV, typename LocalPattern>
      void pattern_volume (const LFSU& lfsu, const LFSV& lfsv,
                           LocalPattern& pattern) const
      {
        for (size_t k=0; k<TypeTree::degree(lfsv); k++)
          for (size_t i=0; i<lfsv.child(k).size(); ++i)
            for (size_t j=0; j<lfsu.child(k).size(); ++j)
              pattern.addLink(lfsv.child(k),i,lfsu.child(k),j);
      }

      // volume integral depending on test and ansatz functions
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // get types
        using namespace Indices;
        using RF = typename NUMFLUX::RF; // type for computations

        // get local function space that is identical for all components
        const auto& dgspace = child(lfsv,_0);

        // get geometry
        auto geo = eg.geometry();

        // Initialize vectors outside for loop
        Dune::FieldVector<RF,m> u(0.0);

        // loop over quadrature points
        const int order = dgspace.finiteElement().localBasis().order();
        const int intorder = overintegration+2*order;
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // evaluate basis functions
            auto& phi = cache[order].evaluateFunction(ip.position(),dgspace.finiteElement().localBasis());

            // evaluate u
            u = 0.0;
            for (size_t k=0; k<m; k++) // for all components
              for (size_t j=0; j<dgspace.size(); j++) // for all basis functions
                u[k] += x(lfsv.child(k),j)*phi[j];

            // integrate
            auto factor = ip.weight() * geo.integrationElement(ip.position());
            for (size_t k=0; k<m; k++) // for all components
              for (size_t i=0; i<dgspace.size(); i++) // for all test functions of this component
                r.accumulate(lfsv.child(k),i, u[k]*phi[i]*factor);
          }
      }

      // jacobian of volume term
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                            M & mat) const
      {
        // get types
        using namespace Indices;

        // get local function space that is identical for all components
        const auto& dgspace = child(lfsv,_0);

        // get geometry
        auto geo = eg.geometry();

        // loop over quadrature points
        const int order = dgspace.finiteElement().localBasis().order();
        const int intorder = overintegration+2*order;
        for (const auto& ip : quadratureRule(geo,intorder))
          {
            // evaluate basis functions
            auto& phi = cache[order].evaluateFunction(ip.position(),dgspace.finiteElement().localBasis());

            // integrate
            auto factor = ip.weight() * geo.integrationElement(ip.position());
            for (size_t k=0; k<m; k++) // for all components
              for (size_t i=0; i<dgspace.size(); i++) // for all test functions of this component
                for (size_t j=0; j<dgspace.size(); j++) // for all ansatz functions of this component
                  mat.accumulate(lfsv.child(k),i,lfsu.child(k),j, phi[j]*phi[i]*factor);
          }
      }

    private:
      NUMFLUX& numflux;
      int overintegration;
      using LocalBasisType = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
      using Cache = Dune::PDELab::LocalBasisCache<LocalBasisType>;
      std::vector<Cache> cache;
    };

  }
}

#endif // DUNE_PDELAB_LOCALOPERATOR_HYPERBOLICDG_HH
