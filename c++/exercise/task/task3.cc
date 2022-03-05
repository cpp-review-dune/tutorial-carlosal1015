#include <iostream>
#include <vector>
#include <cmath>

// Replace this with a template function
int cubed(int x)
{
  return x * x * x;
}


int main()
{
  {
    // create vector with ints
    auto v = std::vector<int>{ 3, 5, 1, 9 };

    // print out contents of v
    for (auto x : v)
      std::cout << x << " ";
    std::cout << std::endl;

    // calculate cubes of entries
    for (auto& x : v)
      x = {cubed(x)};

    // print out cubes
    for (auto x : v)
      std::cout << x << " ";
    std::cout << std::endl;

  }


  {
    // create vector with doubles
    auto v = std::vector<double>{ 3.1, 5.2, 1.6, 8.4 };

    // print out contents of v
    for (auto x : v)
      std::cout << x << " ";
    std::cout << std::endl;

    // calculate cubes of entries
    for (auto& x : v)
      x = cubed(x);

    // print out cubes
    for (auto x : v)
      std::cout << x << " ";
    std::cout << std::endl;

  }


  return 0;
}
