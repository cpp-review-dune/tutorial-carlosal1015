#ifndef NUMERICALFLUX
#define NUMERICALFLUX

//local Lax-Friedrichs Flux
template<typename MODEL>
class LLFflux
{

public:

  static constexpr int dim = MODEL::dim;
  static constexpr int m = MODEL::m;

  using Model = MODEL;
  using RF = typename MODEL::RangeField; // type for computations
  using DF = RF;

  LLFflux (const MODEL& model)
    : model_(model)
  {
  }

  /// tex: llf
  template<typename E, typename X>
  void numericalFlux(const E& inside, const X& x_inside,
                     const E& outside, const X& x_outside,
                     const Dune::FieldVector<DF,dim> n_F,
                     const Dune::FieldVector<RF,m>& u_s,
                     const Dune::FieldVector<RF,m>& u_n,
                     Dune::FieldVector<RF,m>& f) const
  {
    f = 0.0;
    // fetch flux
    Dune::FieldMatrix<RF,m,dim> Fs;
    Dune::FieldMatrix<RF,m,dim> Fn;

    //evaluate flux
    model().flux(inside,x_inside,u_s,Fs);
    model().flux(outside,x_outside,u_n,Fn);

    //Fs*n_F + Fn*n_F
    Fs.umv(n_F,f);
    Fn.umv(n_F,f);
    f *= 0.5;

    //max eigenvalue
    RF alpha(0.0);
    model().max_eigenvalue(inside,x_inside,outside,x_outside,u_s,u_n,n_F,alpha);

    //add diffusion
    for (size_t i =0 ; i<m;i++)
      f[i] = f[i] + 0.5*alpha*(u_s[i] - u_n[i]);
  }
  /// tex: llf

  const MODEL& model() const
  {
    return model_;
  }

private:

  const MODEL& model_;

};// LLF

//Flux Vector splitting
template<typename MODEL>
class FluxVectorSplitting
{

public:

  static constexpr int dim = MODEL::dim;
  static constexpr int m = MODEL::m;

  using Model = MODEL;
  using RF = typename MODEL::RangeField; // type for computations
  using DF = RF;

  FluxVectorSplitting (const MODEL& model)
    : model_(model)
  {
  }

  /// tex: fvs
  template<typename E, typename X>
  void numericalFlux(const E& inside, const X& x_inside,
                     const E& outside, const X& x_outside,
                     const Dune::FieldVector<DF,dim> n_F,
                     const Dune::FieldVector<RF,m>& u_s,
                     const Dune::FieldVector<RF,m>& u_n,
                     Dune::FieldVector<RF,m>& f) const
  {
    Dune::FieldMatrix<DF,m,m> D(0.0);
    // fetch eigenvalues
    model().diagonal(inside,x_inside,D);

    Dune::FieldMatrix<DF,m,m> Dplus(0.0);
    Dune::FieldMatrix<DF,m,m> Dminus(0.0);

    for (size_t i =0 ; i<m;i++)
      (D[i][i] > 0)
        ? Dplus[i][i] = D[i][i]
        : Dminus[i][i] = D[i][i];

    // fetch eigenvectors
    Dune::FieldMatrix<DF,m,m> Rot;
    model().eigenvectors(inside,x_inside,n_F,Rot);

    // compute B+ = RD+R^-1 and B- = RD-R^-1
    Dune::FieldMatrix<DF,m,m> Bplus(Rot);
    Dune::FieldMatrix<DF,m,m> Bminus(Rot);

    //multiply by D+-
    Bplus.rightmultiply(Dplus);
    Bminus.rightmultiply(Dminus);

    //multiply by R^-1
    Rot.invert();
    Bplus.rightmultiply(Rot);
    Bminus.rightmultiply(Rot);

    // Compute numerical flux at  the integration point
    f = 0.0;
    // f = Bplus*u_s + Bminus*u_n
    Bplus.umv(u_s,f);
    Bminus.umv(u_n,f);
  }
  /// tex: fvs

  const MODEL& model() const
  {
    return model_;
  }

private:

  const MODEL& model_;

};// FVS


//Flux Vector splitting for discontinuous coefficients
template<typename MODEL>
class VariableFluxVectorSplitting
{

public:

  static constexpr int dim = MODEL::dim;
  static constexpr int m = MODEL::m;
  static constexpr int mstar = MODEL::mstar;

  using Model = MODEL;
  using RF = typename MODEL::RangeField; // type for computations
  using DF = RF;

  VariableFluxVectorSplitting (const MODEL& model)
    : model_(model)
    , fluxVectorSplitting_(model)
  {
  }

  template<typename E, typename X>
  void numericalFlux(const E& inside, const X& x_inside,
                     const E& outside, const X& x_outside,
                     const Dune::FieldVector<DF,dim> n_F,
                     const Dune::FieldVector<RF,m>& u_s,
                     const Dune::FieldVector<RF,m>& u_n,Dune::FieldVector<RF,m>& f) const
  {
    // check for discontinuity
   if ( model().problem.material(inside,x_inside) == model().problem.material(outside,x_outside)  )
    {
      fluxVectorSplitting_.numericalFlux(inside, x_inside, outside, x_outside, n_F, u_s, u_n, f);
    }
    else // discontinuous coefficient case
    {
      // fetch eigenvalues
      Dune::FieldMatrix<DF,m,m> D_s(0.0), D_n(0.0);
      model().diagonal(inside,x_inside,D_s);
      model().diagonal(outside,x_outside,D_n);

      // split positive and negative eigenvalues
      Dune::FieldMatrix<DF,m,m> Dplus_s(0.0);
      Dune::FieldMatrix<DF,m,m> Dminus_n(0.0);
      for (size_t i=0 ; i<m; i++)
        (D_s[i][i] > 0) ? Dplus_s[i][i] = D_s[i][i] : Dminus_n[i][i] = D_n[i][i];

      // fetch eigenvectors
      Dune::FieldMatrix<DF,m,m> R_s, R_n;
      model().eigenvectors(inside,x_inside,n_F,R_s);
      model().eigenvectors(outside,x_outside,n_F,R_n);

      // compute B+ = RD+R^-1 and B- = RD-R^-1
      Dune::FieldMatrix<DF,m,m> Bplus_s(R_s);
      Dune::FieldMatrix<DF,m,m> Bminus_n(R_n);

      //multiply by D+-
      Bplus_s.rightmultiply(Dplus_s);
      Bminus_n.rightmultiply(Dminus_n);

      //multiply by R^-1
      Dune::FieldMatrix<DF,m,m> Rinv_s(R_s), Rinv_n(R_n);
      Rinv_s.invert(); Rinv_n.invert();
      Bplus_s.rightmultiply(Rinv_s);
      Bminus_n.rightmultiply(Rinv_n);

      // compute rectangular S matrix
      Dune::FieldMatrix<DF,m,mstar> S(0.0);
      for (int j=0; j<mstar; j++)
        if (D_s[j][j]>0.0)
          {
            for (int i=0; i<m; i++) S[i][j] = R_n[i][j]*D_n[j][j];
          }
        else
          {
            for (int i=0; i<m; i++) S[i][j] = -R_s[i][j]*D_s[j][j];
          }

      // compute square normal matrix
      Dune::FieldMatrix<DF,mstar,mstar> SS(0.0);
      for (int i=0; i<mstar; i++)
        for (int j=0; j<mstar; j++)
          for (int k=0; k<m; k++)
            SS[i][j] += S[k][i]*S[k][j];

      // compute right hand side of interface problem
      Dune::FieldVector<RF,m> fd(0.0);
      Bplus_s.umv(u_s,fd);
      Bminus_n.usmv(-1.0,u_n,fd);
      Dune::FieldVector<RF,mstar> rhs(0.0);
      for (int i=0; i<mstar; i++)
        for (int j=0; j<m; j++)
          rhs[i] += S[j][i]*fd[j];

      // Solve interface system
      Dune::FieldVector<RF,mstar> alpha(0.0);
      SS.solve(alpha,rhs);

      // extend alpha
      Dune::FieldVector<RF,m> alpha_ext(0.0);
      for (int i=0; i<mstar; i++)
        if (D_s[i][i]<0.0) alpha_ext[i] = D_s[i][i]*alpha[i];

      // compute flux
      f = 0.0;
      Bplus_s.umv(u_s,f);
      R_s.umv(alpha_ext,f);
    }
  }

  const MODEL& model() const
  {
    return model_;
  }

private:

  const MODEL& model_;
  FluxVectorSplitting<MODEL> fluxVectorSplitting_;

};// VFVS

#endif //NUMERICALFLUX
