#ifndef SEGTREE_H_
#define SEGTREE_H_

#include <cstdio>
#include <vector>
#include <cmath>
#include "node.h"
//using namespace std;



namespace grf {

class segtree {
public:
  std::vector<node> T;
  std::vector<double> X;
  static double eps;
  
  int N, n; // n is the number of spaces between points of cdf (i.e. n+1 is the total number of discontinuities of F and G)
  // N is the first power of 2 that is >= n
  
  double n_F, n_G; //count the number of points in F and G
  
  //initialization
  segtree(std::vector <double> x) { // x is the vector of future x-values
    n_F = n_G = 0;
    
    // sort X and keep only unique values
    X = x;
    sort( X.begin(), X.end() );
    X.erase( unique( X.begin(), X.end() ), X.end() );
    
    n = X.size() - 1; // one leaf represents space between two values of X
    for(N = 1; N < n; N = N*2); // compute first power of 2 that is >= n
    
    // tree has in total 2*N nodes, 1 is the index of the root, node k has left and right children 2*k and 2*k+1 respectively
    T.resize(2*N);
    
    // initialization of the values in the leaves
    for(int i = 0; i < N; ++i){
      if(i < n) T[N + i].weight_sum = X[i + 1] - X[i]; // size of space between two values of X, other weights remain zero
      T[N + i].left = T[N + i].right = i;
    }
    
    // initialization of other nodes from their children, note that always the children have already been computed
    for(int i = N - 1; i > 0; --i){
      T[i].weight_sum = T[2*i].weight_sum + T[2*i + 1].weight_sum;
      T[i].left = T[2*i].left;
      T[i].right = T[2*i + 1].right;
    }
  }
  
  void propagate(int idx){ // propagate the change from the index idx to its children and reset the propagation counter
    if(abs(T[idx].F_prop) > eps){
      update_node_value(2*idx, T[idx].F_prop, true);
      update_node_value(2*idx + 1, T[idx].F_prop, true);
      T[idx].F_prop = 0;
    }
    if(abs(T[idx].G_prop) > eps){
      update_node_value(2*idx, T[idx].G_prop, false);
      update_node_value(2*idx + 1, T[idx].G_prop, false);
      T[idx].G_prop = 0;
    }
  }
  
  void compute_from_children(int idx){ // compute the values of the variables in node idx from the corresponding variables of its children
    T[idx].F_sum = T[2*idx].F_sum + T[2*idx + 1].F_sum;
    T[idx].G_sum = T[2*idx].G_sum + T[2*idx + 1].G_sum;
    
    T[idx].F_squared = T[2*idx].F_squared + T[2*idx + 1].F_squared;
    T[idx].G_squared = T[2*idx].G_squared + T[2*idx + 1].G_squared;
    
    T[idx].cross = T[2*idx].cross + T[2*idx + 1].cross;
  }
  
  void update_node_value(int idx, double amount, bool F){ // increase all values in the subinterval covered by index idx by 'amount', remember to propagate the changes, bool F denotes if we change F or G
    if(F){
      T[idx].F_prop += amount;
      T[idx].F_squared += 2 * amount * T[idx].F_sum + amount * amount * T[idx].weight_sum;
      T[idx].F_sum += T[idx].weight_sum * amount;
      T[idx].cross += amount * T[idx].G_sum;
    }
    else{
      T[idx].G_prop += amount;
      T[idx].G_squared += 2 * amount * T[idx].G_sum + amount * amount * T[idx].weight_sum;
      T[idx].G_sum += T[idx].weight_sum * amount;
      T[idx].cross += amount * T[idx].F_sum;
    }
  }
  
  void update(int idx, int l, int r, double amount, bool F){ // update recursively the values in the tree, idx is current position, we need to increase all values in the interval [l, r] by 'amount', bool F denotes if we change F or G
    if(T[idx].right < l || T[idx].left > r) return; // if the subinteraval covered by node idx is disjoint from [l, r]], go back
    
    if(l <= T[idx].left && T[idx].right <= r){ // if the subinterval covered by node idx is within [l, r], update node values, memorize that e need to propagate the change later and go back
      update_node_value(idx, amount, F);
      return;
    }
    
    // now we will call recursion for both children so need to propagate the changes that I have not done yet
    propagate(idx);
    update(2*idx, l, r, amount, F);
    update(2*idx + 1, l, r, amount, F);
    compute_from_children(idx); // compute values of node idx from its children
    return;
  }
  
  int find(double x){ // uses binary search to return the index of x in the array X, which is already sorted and has unique values
    int lo = 0;
    int hi = X.size() - 1;
    int mid;
    
    while(lo < hi){
      int mid = (lo + hi) / 2;
      if(X[mid] < x - eps){
        lo = mid + 1;
      }
      else{
        hi = mid;
      }
    }
    return lo;
  }
  
  void insert_point(double x, int which){ // add point x to F (if which = 0) or to G (if which = 1)
    if(which == 0){
      ++n_F;
      update(1, find(x), n-1, 1, true);
    }
    else{
      ++n_G;
      update(1, find(x), n-1, 1, false);
    }
  }
  
  void remove_point(double x, int which){ // remove point x from F (if which = 0) or from G (if which = 1)
    if(which == 0){
      --n_F;
      update(1, find(x), n-1, -1, true);
    }
    else{
      --n_G;
      update(1, find(x), n-1, -1, false);
    }
  }
  
  double energy(){ // return the energy statistic
    double ret = T[1].F_squared / (n_F * n_F) - 2 * T[1].cross / (n_F * n_G) + T[1].G_squared / (n_G * n_G);
    return  ((n_F * n_G) / (n_F + n_G)) * ret;
    //  *
  }
};

}

#endif /* SEGTREE_H_ */
