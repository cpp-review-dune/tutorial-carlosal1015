#ifndef ACOUSTICS_RIEMANNPROBLEM
#define ACOUSTICS_RIEMANNPROBLEM
template<typename GV, typename NUMBER>
class Problem
{
public:

  using RangeField = NUMBER;

  //problem specification depends on dimension
  static constexpr int dim = GV::dimension;
  static constexpr int m = dim+1;

  using Range = Dune::FieldVector<NUMBER,m>;

  Problem ()
    : time(0.0)
  {
  }

  /// tex: material
  //! material 
  // this function is used to decide if we work
  // with discontinous coefficient case  
  template<typename E, typename X>
  int material (const E& e, const X& x) const 
  {
    auto xglobal = e.geometry().center();
    if (xglobal[0]>1.625)
      return 1;
    else
      return 2;
  }
  /// tex: material

  /// tex: speedofsound
  //! speed of sound
  template<typename E, typename X>
  NUMBER c (const E& e, const X& x) const
  {
    auto xglobal = e.geometry().center();
    if (xglobal[0]>1.625)
      return 0.33333;
    else
      return 1.0;
  }
  /// tex: speedofsound

  /// tex: bc
  //! Boundary condition value - reflecting bc
  template<typename I, typename X, typename R>
  Range g (const I& is, const X& x, const R& s) const
  {
    Range u(0.0);
    u[0] =  s[0];
    
    for (size_t i=0; i<dim; i++) 
      u[i+1] = -s[i+1];
    
    return u;
  }
  /// tex: bc

  /// tex: rhs
  //! right hand side
  template<typename E, typename X>
  Range q (const E& e, const X& x) const
  {
    Range rhs(0.0);
    return rhs;
  }
  /// tex: rhs

  /// tex: init
  //! initial value -> the same as tutorial04
  template<typename E, typename X>
  Range u0 (const E& e, const X& x) const
  {
    X xglobal = e.geometry().global(x);
    Range u(0.0);
    for (int i=0; i<dim; i++)
      u[0] += (xglobal[i]-0.375)*(xglobal[i]-0.375);
    u[0] = std::max(0.0,1.0-8.0*sqrt(u[0]));

    return u;
  }
  /// tex: init

  //! set time for subsequent evaluation
  void setTime (NUMBER t)
  {
    time = t;
  }

private:

  NUMBER time;

};
#endif //ACOUSTICS_RIEMANNPROBLEM
