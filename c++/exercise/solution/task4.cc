#include <iostream>
#include <vector>


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

  // store lambda in a variable
  auto squared = [](auto x)
    {
      return x * x;
    };
  calculate(v,squared);

  double scale_factor = 2;
  // Scale values of v by scale_factor
  // create lambda function in place and capture the
  // variable scale factor
  calculate(v,[&](auto x)
    {
      return x * scale_factor;
    }
  );

  return 0;
}
