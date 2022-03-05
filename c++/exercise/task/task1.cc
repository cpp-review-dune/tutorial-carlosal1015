#include <iostream>
#include <algorithm>

using namespace std;

template<typename T>
T min(T a, T b)
{
  if (a < b)
    return a;
  else
    return b;
}

int main()
{
  cout << min(4,8) << endl;
  return 0;
}
