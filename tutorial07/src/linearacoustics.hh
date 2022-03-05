#ifndef ACOUSTICS_MODEL
#define ACOUSTICS_MODEL
/*
 Linear Acoustics model class
*/
/** \brief provide matrix which contains rowwise the eigenvectors of linear acoustics problem
    \tparam dim the space dimension
    \param c speed of sound
    \param n unit outer normal vector
    \param RT matrix to be filled
*/

template<typename PROBLEM>
class Model
{
public:
  static constexpr int dim = PROBLEM::dim;   // space dimension
  static constexpr int m = dim+1; // system size
  static constexpr int mplus = 1; // number of positive eigenvalues
  static constexpr int mminus = 1; // number of negative eigenvalues
  static constexpr int mstar = mplus+mminus; // number of nonzero eigenvalues

  using RangeField = typename PROBLEM::RangeField;

  Model (PROBLEM& p)
  : problem(p)
  {
  }


  /// tex: eigenvectors
  template<typename E, typename X, typename T2, typename T3>
  void eigenvectors (const E& e, const X& x,
                     const Dune::FieldVector<T2,dim>& n,
                     Dune::FieldMatrix<T3,m,m>& RT) const
  {
    auto c = problem.c(e,x);

    //TODO find a way to write eigenvectors independently of dim
    if (dim == 1) {
      RT[0][0] =  1; RT[1][0] = c*n[0];
      RT[0][1] = -1; RT[1][1] = c*n[0];
    }
    if (dim == 2) {
      RT[0][0] =  1.0/c;  RT[0][1] = -1.0/c;  RT[0][2] = 0.0;
      RT[1][0] =  n[0]; RT[1][1] = n[0];  RT[1][2] = -n[1];
      RT[2][0] =  n[1]; RT[2][1] = n[1];  RT[2][2] = n[0];
    }
  }
  /// tex: eigenvectors


  /// tex: diagonal
  template<typename E, typename X, typename RF>
  void diagonal (const E& e, const X& x,
                 Dune::FieldMatrix<RF,m,m>& D) const
  {
    auto c = problem.c(e,x);

    for (size_t i=0; i<m; i++)
      for (size_t j=0; j<m; j++)
        D[i][j] = 0.0;
    D[0][0] = c;
    D[1][1] = -c ;
  }
  /// tex: diagonal


  template<typename E, typename X, typename RF>
  void max_eigenvalue (const E& inside, const X& x_inside,
                       const E& outside, const X& x_outside,
                       const Dune::FieldVector<RF,m>& u_s,
                       const Dune::FieldVector<RF,m>& u_n,
                       const Dune::FieldVector<RF,dim>& n_F,
                       RF& alpha) const
  {
    alpha = std::max( problem.c(inside,x_inside),
                      problem.c(outside,x_outside) );
  }

  /// tex: flux
  //Flux function
  template<typename E, typename X, typename RF>
  void flux (const E& e, const X& x,
             const Dune::FieldVector<RF,m>& u,
             Dune::FieldMatrix<RF,m,dim>& F) const
  {
    auto c = problem.c(e,x);

    for (size_t i=0; i<dim; i++) {
      F[0][i] = u[i+1];
      F[i+1][i] = c*c*u[0];
    }

  }
  /// tex: flux

  const PROBLEM& problem;

};// Model

#endif //ACOUSTICS_MODEL
