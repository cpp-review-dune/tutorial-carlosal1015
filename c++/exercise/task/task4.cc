#include <iostream>
#include <vector>


template<typename T>
T squared(T x)
{
  return x * x;
}


template<typename T>
struct Scaled
{

  template<typename U>
  U operator()(U x) const
  {
    return x * factor;
  }

  Scaled(T factor_)
    : factor(factor_)
  {}

  T factor;

};


// Call f for each entry in vec and print the result
template<typename Vector, typename Function>
void calculate(Vector& vec, Function f)
{
  for (auto& x : vec)
    std::cout << f(x) << " ";
  std::cout << std::endl;
}


int main()
{
  // create vector with doubles
  auto v = std::vector<double>{ 3, 5, 1, 9 };

  // print out contents of v
  for (auto x : v)
    std::cout << x << " ";
  std::cout << std::endl;

  // calculate squares of v - using a template function
  // note that we have to provide the number type
  calculate(v,squared<double>);

  double scale_factor = 2;
  // Scale values of v by scale_factor
  calculate(v,Scaled<double>(scale_factor));

  return 0;
}
