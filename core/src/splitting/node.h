#ifndef NODE_H_
#define NODE_H_

namespace grf {

class node {
public:
  
  double weight_sum, F_squared, G_squared, cross, F_sum, G_sum; //information in the node
  double F_prop, G_prop; // amount to increase the children
  int left, right; // the node contains information about the interval [left, right]
  
  node() {
    weight_sum = F_squared = G_squared = cross = F_sum = G_sum = F_prop = G_prop = 0;
  }
};

}

#endif /* NODE_H_ */
