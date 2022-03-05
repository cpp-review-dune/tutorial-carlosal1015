#ifndef SHALLOWWATER_MODEL
#define SHALLOWWATER_MODEL
/*
 Shallow Water model class
*/

#include <cmath>

template<int dim, typename PROBLEM>
class Model_ ;

//wrapper for dimension specialisation
template<typename PROBLEM>
using Model = Model_<PROBLEM::dim, PROBLEM>;

template<typename PROBLEM>
class Model_<1,PROBLEM>
{

public:
  static constexpr int dim = 1;
  static constexpr int m = 2;

  using RangeField = typename PROBLEM::RangeField;

  Model_ (PROBLEM& p)
  : problem(p)
  {
  }

  template<typename E, typename X, typename RF>
  void max_eigenvalue (const E& inside, const X& x_inside,
                       const E& outside, const X& x_outside,
                       const Dune::FieldVector<RF,m>& u_s,
                       const Dune::FieldVector<RF,m>& u_n,
                       const Dune::FieldVector<RF,dim>& n_F,
                       RF& alpha) const
  {
    auto g = problem.gravity(inside,x_inside);

    RF alpha_s(0.0);
    RF alpha_n(0.0);

    using std::abs;
    using std::max;
    using std::sqrt;

    alpha_s = abs(u_s[1]/u_s[0]) + sqrt(g*u_s[0]);
    alpha_n = abs(u_n[1]/u_n[0]) + sqrt(g*u_n[0]);

    alpha = max(alpha_s, alpha_n);
  }

  //Flux function
  template<typename E, typename X, typename RF>
  void flux (const E& e, const X& x,
             const Dune::FieldVector<RF,m>& u,
             Dune::FieldMatrix<RF,m,dim>& F) const
  {
    auto g = problem.gravity(e,x);
    F[0][0] = u[1];
    F[1][0] = u[1]*u[1]/u[0] + 0.5*g*u[0]*u[0];
  }

  const PROBLEM& problem;
};


template<typename PROBLEM>
class Model_<2,PROBLEM>
{

public:
  static constexpr int dim = 2;
  static constexpr int m = 3;

  using RangeField = typename PROBLEM::RangeField;

  Model_ (PROBLEM& p)
  : problem(p)
  {
  }

  template<typename E, typename X, typename RF>
  void max_eigenvalue (const E& inside, const X& x_inside,
                       const E& outside, const X& x_outside,
                       const Dune::FieldVector<RF,m>& u_s,
                       const Dune::FieldVector<RF,m>& u_n,
                       const Dune::FieldVector<RF,dim>& n_F,
                       RF& alpha) const
  {
    auto g = problem.gravity(inside,x_inside);

    RF alpha_s(0.0);
    RF alpha_n(0.0);

    for(size_t k=0;k<dim;++k)
    {
      alpha_s +=  n_F[k]*u_s[k+1];
      alpha_n += -n_F[k]*u_n[k+1];
    }

    using std::abs;
    using std::max;
    using std::sqrt;

    alpha_s = abs(alpha_s) / u_s[0] + sqrt(g*u_s[0]);
    alpha_n = abs(alpha_n) / u_n[0] + sqrt(g*u_n[0]);

    alpha = max(alpha_s, alpha_n);
  }


  //Flux function
  template<typename E, typename X, typename RF>
  void flux (const E& e, const X& x, const Dune::FieldVector<RF,m>& u, Dune::FieldMatrix<RF,m,dim>& F) const
  {
    auto g = problem.gravity(e,x);
    F[0][0] = u[1]                            ; F[0][1] = u[2];
    F[1][0] = u[1]*u[1]/u[0] + 0.5*g*u[0]*u[0]; F[1][1] = u[1]*u[2]/u[0];
    F[2][0] = u[1]*u[2]/u[0]                  ; F[2][1] = u[2]*u[2]/u[0] + 0.5*g*u[0]*u[0];
  }

  const PROBLEM& problem;
};

#endif //SHALLOWWATER_MODEL
