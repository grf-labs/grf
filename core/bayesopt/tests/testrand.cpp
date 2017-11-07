#include "randgen.hpp"

int main()
{
  randEngine eng;
  randFloat sample(eng, realUniformDist(0,1));

  for (int i = 0; i<1000; i++)
    std::cout << sample() << std::endl;

  return 0;
}
