#include <boost/numeric/ublas/assignment.hpp>
#include "gridsampling.hpp"

int main()
{
  vectori divs(3); 
  divs <<= 3,4,5;

  vecOfvec grid;
  bayesopt::utils::buildGrid(divs,grid);

  for(size_t i=0; i<grid.size(); ++i)
    {
      for(size_t j=0; j<grid[i].size(); ++j)
	{
	  std::cout << grid[i](j) << ", ";	  
	}
      std::cout << std::endl;
    }
}
