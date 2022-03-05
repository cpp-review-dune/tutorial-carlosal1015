#ifndef MAXWELL_RIEMANNPROBLEM
#define MAXWELL_RIEMANNPROBLEM
template<typename GV, typename NUMBER>
class Problem
{
public:

  using RangeField = NUMBER;

  //problem specification depends on dimension
  static constexpr int dim = 3;
  static constexpr int m = 6;

  using Range = Dune::FieldVector<NUMBER,m>;

  Problem ()
    : time(0.0)
  {
  }

  //! material 
  // this function is used to decide if we work
  // with discontinous coefficient case  
  template<typename E, typename X>
  int material (const E& e, const X& x) const 
  {
    return 1;
  }


  //! permittivity
  template<typename E, typename X>
  NUMBER eps (const E& e, const X& x) const
  {
    return 1.0;
  }

  //! permeability
  template<typename E, typename X>
  NUMBER mu (const E& e, const X& x) const
  {
    return 1.0;
  }

  //! speed of sound
  template<typename E, typename X>
  NUMBER c (const E& e, const X& x) const
  {
    return 1.0;
  }

  //! Boundary condition 
  template<typename I, typename X, typename R>
  Range g (const I& is, const X& x, const R& s) const
  {
    Range u(0.0);
    //reflecting bc
    // u[0] = -s[0];
    // u[1] = -s[1];
    // u[2] = -s[2];
    // u[3] = -s[3];
    // u[4] = -s[4];
    // u[5] = -s[5];
    return u;
  }

  //! right hand side
  template<typename E, typename X>
  Range q (const E& e, const X& x) const
  {
    Range rhs(0.0);
    return rhs;
  }

  //! initial value
  template<typename E, typename X>
  Range u0 (const E& e, const X& p) const
  {
    auto xglobal = e.geometry().global(p);


    Range u(0.0);

    auto x=xglobal[0];
    
    auto alpha = 1.0;

    auto s1 = std::sin(alpha*x);
    auto c1 = std::cos(alpha*x);
    auto st = std::sin(alpha*0.0);
    auto ct = std::cos(alpha*0.0);
    u[0] += 0;       // E_x
    u[1] += -c1*st;  // E_y
    u[2] += s1*ct;   // E_z
    u[3] += 0;       // H_x
    u[4] += c1*st;   // H_y
    u[5] += s1*ct;   // H_z

    return u;
  }

  //! set time for subsequent evaluation
  void setTime (NUMBER t)
  {
    time = t;
  }

private:

  NUMBER time;

};
#endif //MAXWELL_RIEMANNPROBLEM
