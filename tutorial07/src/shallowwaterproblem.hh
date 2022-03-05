#ifndef SHALLOWWATER_RIEMANNPROBLEM
#define SHALLOWWATER_RIEMANNPROBLEM
template<typename GV, typename NUMBER>
class Problem
{
public:

  using RangeField = NUMBER;

  //problem specification depends on dimension
  static constexpr int dim = GV::dimension;
  static constexpr int m = dim+1;

  using Range = Dune::FieldVector<NUMBER,m>;

  Problem (RangeField gravity = 1.0)
    : _time(0.0)
    , _gravity(gravity)
  {
  }

  //! Boundary condition value - reflecting bc
  template<typename I, typename X, typename R>
  Range g (const I& is, const X& x, const R& s) const
  {
    Range u(0.0);
    u[0] =  s[0];
    for (int i=0; i<dim; i++)
       u[i+1] = -s[i+1];

    return u;
  }

  //! right hand side
  template<typename E, typename X>
  Range q (const E& e, const X& x) const
  {
    Range rhs(0.0);
    return rhs;
  }

  //! initial value -> the same as tutorial04
  template<typename E, typename X>
  Range u0 (const E& e, const X& x) const
  {
    X xglobal = e.geometry().global(x);
    Range u(0.0);

    //hump
    using std::pow; using std::exp;
    auto tmp = 0.0;
    for (int i=0; i<dim; i++)
      tmp += pow((xglobal[i]-0.5)/0.2,2)/2;

    u[0] += exp(-tmp) + 0.5;

    return u;
  }

  //! gravity
  template<typename E, typename X>
  RangeField gravity (const E& e, const X& x) const
  {
    return _gravity;
  }

  //! set time for subsequent evaluation
  void setTime (NUMBER t)
  {
    _time = t;
  }

private:

  NUMBER _time;
  RangeField _gravity;

};
#endif //SHALLOWWATER_RIEMANNPROBLEM
