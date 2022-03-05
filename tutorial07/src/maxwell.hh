#ifndef MAXWELL_MODEL
#define MAXWELL_MODEL
/*
 Maxwell model class
*/
/** \brief provide matrix which contains rowwise the eigenvectors of linear acoustics problem
    \tparam dim the space dimension
    \param c speed of sound
    \param n unit outer normal vector
    \param RT matrix to be filled
*/

template<int dim, typename PROBLEM>
class Model ;


template<typename PROBLEM>
class Model<3,PROBLEM>
{

public:
  static constexpr int dim = 3;
  static constexpr int m = 6;

  using RangeField = typename PROBLEM::RangeField;

  Model (PROBLEM& p)
    : problem(p)
  {
  }

  template<typename E, typename X, typename T1, typename T2>
  void eigenvectors (const E& e, const X& x, const Dune::FieldVector<T1,dim>& n, Dune::FieldMatrix<T2,m,m>& R) const
  {
    int eps = 1.0;
    int mu  = 1.0;

    auto a=n[0], b=n[1], c=n[2];

    Dune::FieldVector<T2,dim> alpha, beta;

    if (std::abs(c)<0.5)
      {
        alpha[0]=a*c; alpha[1]=b*c; alpha[2]=c*c-1;
        beta[0]=-b;  beta[1]=a;   beta[2]=0;
      }
    else
      {
        alpha[0]=a*b; alpha[1]=b*b-1; alpha[2]=b*c;
        beta[0]=c;  beta[1]=0.0;   beta[2]=-a;
      }

    // \lambda_0,1 = s
    R[0][0] =  alpha[0];   R[0][1] =  -beta[0];
    R[1][0] =  alpha[1];   R[1][1] =  -beta[1];
    R[2][0] =  alpha[2];   R[2][1] =  -beta[2];
    R[3][0] =  beta[0];    R[3][1] =  alpha[0];
    R[4][0] =  beta[1];    R[4][1] =  alpha[1];
    R[5][0] =  beta[2];    R[5][1] =  alpha[2];

    // \lambda_2,3 = -s
    R[0][2] =  beta[0];   R[0][3] =  alpha[0];
    R[1][2] =  beta[1];   R[1][3] =  alpha[1];
    R[2][2] =  beta[2];   R[2][3] =  alpha[2];
    R[3][2] =  alpha[0];  R[3][3] =  -beta[0];
    R[4][2] =  alpha[1];  R[4][3] =  -beta[1];
    R[5][2] =  alpha[2];  R[5][3] =  -beta[2];

    // \lambda_4,5 = 0
    R[0][4] =   a;  R[0][5] =   0;
    R[1][4] =   b;  R[1][5] =   0;
    R[2][4] =   c;  R[2][5] =   0;
    R[3][4] =   0;  R[3][5] =   a;
    R[4][4] =   0;  R[4][5] =   b;
    R[5][4] =   0;  R[5][5] =   c;

    // apply scaling
    T1 weps=sqrt(eps);
    T1 wmu=sqrt(mu);
    for (std::size_t i=0; i<3; i++)
      for (std::size_t j=0; j<6; j++)
        R[i][j] *= weps;
    for (std::size_t i=3; i<6; i++)
      for (std::size_t j=0; j<6; j++)
        R[i][j] *= wmu;

    return;
  };
  //one can also provide eigenvectors inverse

  template<typename RF>
  static void coefficients (Dune::FieldMatrix<RF,m,m>& A)
  {
    RF mu(1.0);
    RF ep(1.0);

    A[0][0] = 0.0;   A[0][1] = 0.0;   A[0][2] = 0.0;   A[0][3] = 0.0;   A[0][4] = 1./mu; A[0][5] =-1./mu;
    A[1][0] = 0.0;   A[1][1] = 0.0;   A[1][2] = 0.0;   A[1][3] =-1./mu; A[1][4] = 0.0;   A[1][5] = 1./mu;
    A[2][0] = 0.0;   A[2][1] = 0.0;   A[2][2] = 0.0;   A[2][3] = 1./mu; A[2][4] =-1./mu; A[2][5] = 1.0;
    A[3][0] = 0.0;   A[3][1] =-1./mu; A[3][2] = 1./mu; A[3][3] = 0.0;   A[3][4] = 1.0;   A[3][5] = 1.0;
    A[4][0] = 1./mu; A[4][1] = 0.0;   A[4][2] =-1./mu; A[4][3] = 0.0;   A[4][4] = 1.0;   A[4][5] = 1.0;
    A[5][0] =-1./mu; A[5][1] = 1./mu; A[5][2] = 0.0;   A[5][3] = 0.0;   A[5][4] = 1.0;   A[5][5] = 1.0;
  }

  template<typename E, typename X, typename RF>
  void diagonal (const E& e, const X& x, Dune::FieldMatrix<RF,m,m>& D) const
  {
    RF c(1.0);

    D[0][0] = c;   D[0][1] = 0.0; D[0][2] = 0.0; D[0][3] = 0.0; D[0][4] = 0.0; D[0][5] = 0.0;
    D[1][0] = 0.0; D[1][1] = c;   D[1][2] = 0.0; D[1][3] = 0.0; D[1][4] = 0.0; D[1][5] = 0.0;
    D[2][0] = 0.0; D[2][1] = 0.0; D[2][2] = -c;  D[2][3] = 0.0; D[2][4] = 0.0; D[2][5] = 0.0;
    D[3][0] = 0.0; D[3][1] = 0.0; D[3][2] = 0.0; D[3][3] = -c;  D[3][4] = 0.0; D[3][5] = 0.0;
    D[4][0] = 0.0; D[4][1] = 0.0; D[4][2] = 0.0; D[4][3] = 0.0; D[4][4] = 0.0; D[4][5] = 0.0;
    D[5][0] = 0.0; D[5][1] = 0.0; D[5][2] = 0.0; D[5][3] = 0.0; D[5][4] = 0.0; D[5][5] = 0.0;
  }

  template<typename E, typename X, typename RF>
  void max_eigenvalue (const E& inside, const X& x_inside,
                       const E& outside, const X& x_outside,
                       const Dune::FieldVector<RF,m>& u_s,
                       const Dune::FieldVector<RF,m>& u_n,
                       const Dune::FieldVector<RF,dim>& n_F,
                       RF& alpha) const
  {
    using std::sqrt;
    using std::max;
    RangeField c_inside = 1./sqrt(problem.eps(inside,x_inside) * problem.mu(inside,x_inside));
    RangeField c_outside = 1./sqrt(problem.eps(outside,x_outside) * problem.mu(outside,x_outside));
    alpha = max(c_inside,c_outside);
  }

  //Flux function
  template<typename E, typename X, typename RF>
  void flux (const E& e, const X& x, const Dune::FieldVector<RF,m>& u, Dune::FieldMatrix<RF,m,dim>& F) const
  {
    RF mu(1.0);
    RF ep(1.0);

    F[0][0] = 0.0      ; F[0][1] =-1/mu*u[5]; F[0][2] = 1/mu*u[4];
    F[1][0] = 1/mu*u[5]; F[1][1] = 0.0;       F[1][2] =-1/mu*u[3];
    F[2][0] =-1/mu*u[4]; F[2][1] = 1/mu*u[3]; F[2][2] = 0.0;
    F[3][0] = 0.0      ; F[3][1] = 1/ep*u[2]; F[3][2] =-1/ep*u[1];
    F[4][0] =-1/ep*u[2]; F[4][1] = 0.0;       F[4][2] = 1/ep*u[0];
    F[5][0] = 1/ep*u[1]; F[5][1] =-1/ep*u[0]; F[5][2] = 0.0;
  }

  const PROBLEM& problem;

};
#endif //MAXWELL_MODEL
