#include <vector>

int f() {
  return 1;
}

int main()
{
  auto var1 = 5678;   // var1 has type int
  auto var2 = 'x';    // var2 has type char
  auto var3 = f();    // var3 has return type of f() i.e. int
  auto vector = std::vector<double>(5);
}
