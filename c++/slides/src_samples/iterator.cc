#include<iostream>
#include<vector>

int main(int argc, char** argv)
{
  std::vector<int> vec;
  for (int i=0; i<10 ;++i)
    vec.push_back(i);

  for (auto it=vec.begin(); it!=vec.end(); ++it)
    std::cout << *it << " ";
  std::cout << std::endl;
}
