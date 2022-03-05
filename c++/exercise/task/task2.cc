#include <iostream>
#include <vector>

double cubed(double x)
{
  return x * x * x;
}


int main()
{
  // create vector with doubles
  auto v = std::vector<double>{ 3, 5, 1, 9 };

  // print out contents of v
  for (auto x : v)
    std::cout << x << " ";
  std::cout << std::endl;

  // calculate cubes of entries
  for (auto x : v)
    x = cubed(x);

  // print out cubes
  for (auto x : v)
    std::cout << x << " ";
  std::cout << std::endl;

  return 0;
}
