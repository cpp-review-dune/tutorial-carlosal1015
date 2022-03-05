template<typename Number>
class Problem
{
  Number eta;
  Number t;

public:
  typedef Number value_type;

  //! Constructor without arg sets nonlinear term to zero
  Problem () : eta(0.0), t(0.0) {}

  //! Constructor takes eta parameter
  Problem (const Number& eta_) : eta(eta_), t(0.0) {}

  //! nonlinearity
  Number q (Number u) const
  {
    return 0.0;
  }

  //! derivative of nonlinearity
  Number qprime (Number u) const
  {
    return 0.0;
  }

  //! right hand side
  template<typename E, typename X>
  Number f (const E& e, const X& x) const
  {
    return 0.0;
  }

  //! boundary condition type function (true = Dirichlet)
  template<typename I, typename X>
  bool b (const I& i, const X& x) const
  {
    return false;
  }

  //! Dirichlet extension
  template<typename E, typename X>
  Number g (const E& e, const X& x) const
  {
    auto global = e.geometry().global(x);
    Number s=1.0;
    using std::abs;
    using std::min;
    using std::max;
    for (std::size_t i=0; i<global.size(); i++) {
      Number f = 0.5-8*(abs(global[i]-0.5) - 0.25);
      s *= min(1.0,max(0.0,f));
    }
    return s;
  }

  //! Neumann boundary condition
  template<typename I, typename X>
  Number j (const I& i, const X& x) const
  {
    return 0.0;
  }

  //! Set time in instationary case
  void setTime (Number t_)
  {
    t = t_;
  }
};
