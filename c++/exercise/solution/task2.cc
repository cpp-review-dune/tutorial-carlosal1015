#include <iostream>
#include <vector>

int cubed(int x)
{
  return x * x * x;
}


int main()
{
  // create vector with doubles
  auto v = std::vector<int>{ 3, 5, 1, 9 };

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

  return 0;
}
